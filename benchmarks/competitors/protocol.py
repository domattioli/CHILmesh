"""Standard-library control protocol for an isolated measured region."""
import copy
import hashlib
import importlib.metadata
import os
import json
import math
import platform
import sysconfig
import sys
import time
from contextlib import contextmanager


class ProtocolError(ValueError):
    """A control message or transition is invalid."""


def clock_metadata():
    info = time.get_clock_info("perf_counter")
    return {"name": "perf_counter", "monotonic": info.monotonic,
            "implementation": info.implementation, "resolution_s": info.resolution,
            "adjustable": info.adjustable}


class ProtocolState:
    sequence = ("READY", "GO", "START", "END", "ACK")

    def __init__(self, repeat_id):
        if not isinstance(repeat_id, str) or not repeat_id.strip():
            raise ProtocolError("missing repeat identity")
        self.repeat_id = repeat_id
        self.events = []
        self.receipts = []

    @property
    def complete(self):
        return len(self.events) == len(self.sequence)

    def finish(self):
        if not self.complete:
            raise ProtocolError("incomplete protocol sequence")
        return self.events

    def accept(self, message, *, received_s=None):
        if not isinstance(message, dict) or type(message.get("protocol_version")) is not int or message.get("protocol_version") != 1 or message.get("repeat_id") != self.repeat_id or message.get("clock") != "perf_counter":
            raise ProtocolError("invalid protocol identity or clock")
        metadata = message.get("clock_metadata")
        if (not isinstance(metadata, dict) or metadata.get("monotonic") is not True
                or not isinstance(metadata.get("implementation"), str) or not metadata["implementation"]
                or type(metadata.get("resolution_s")) not in (float, int)
                or not math.isfinite(metadata["resolution_s"]) or metadata["resolution_s"] <= 0):
            raise ProtocolError("invalid clock metadata")
        if self.events and metadata != self.events[0]["clock_metadata"]:
            raise ProtocolError("clock metadata changed")
        stamp = message.get("monotonic_s")
        if isinstance(stamp, bool) or not isinstance(stamp, (float, int)) or not math.isfinite(stamp):
            raise ProtocolError("invalid timestamp")
        if self.complete or message.get("event") != self.sequence[len(self.events)]:
            raise ProtocolError("unexpected protocol event")
        if self.events and stamp < self.events[-1]["monotonic_s"]:
            raise ProtocolError("timestamps out of order")
        if received_s is not None and (type(received_s) not in (int, float) or not math.isfinite(received_s) or received_s < stamp):
            raise ProtocolError("invalid receipt timestamp")
        self.events.append(copy.deepcopy(message))
        self.receipts.append(received_s)


class ChildProtocol:
    def __init__(self, repeat_id, incoming, outgoing, result_stream, clock=time.perf_counter):
        self.repeat_id, self.incoming, self.outgoing, self.result_stream = repeat_id, incoming, outgoing, result_stream
        self.clock = clock
        self.state = ProtocolState(repeat_id)
        self.new_modules = None
        self.result_written = False
        self.runtime = child_runtime_identity()

    def _send(self, event):
        message = dict(protocol_version=1, repeat_id=self.repeat_id, event=event, monotonic_s=self.clock(), clock="perf_counter", clock_metadata=clock_metadata())
        self.state.accept(message)
        self.outgoing.write(json.dumps(message, allow_nan=False) + "\n")
        self.outgoing.flush()

    def _receive(self):
        try:
            message = json.loads(self.incoming.readline())
            self.state.accept(message)
        except (ValueError, TypeError) as exc:
            raise ProtocolError("invalid control message") from exc

    def ready(self):
        self._send("READY")
        self._receive()

    @contextmanager
    def region(self):
        modules = set(sys.modules)
        self._send("START")
        yield
        self.new_modules = sorted(set(sys.modules) - modules)
        self._send("END")
        self._receive()

    def result(self, data):
        if not self.state.complete or self.new_modules is None or self.result_written:
            raise ProtocolError("result before END acknowledgement")
        self.result_stream.write(json.dumps({**data, "new_modules": self.new_modules, "runtime": self.runtime}, allow_nan=False) + "\n")
        self.result_stream.flush()
        self.result_written = True

    def run(self, operation, *, materialize, verify):
        """Run prepared work, eagerly materialize inside the region, verify after ACK.

        Call ready() after imports and prerequisites. The materializer must return
        the complete actual output; the verifier must not reparse the input.
        """
        with self.region():
            output = materialize(operation())
        data = verify(output)
        self.result(data)
        return output


class MappingUnavailable(ValueError):
    """An adapter cannot establish the original identifiers."""


def canonical_output(output):
    """Normalize explicit original IDs or a declared row-offset map."""
    required = ('node_ids', 'element_ids', 'coordinates', 'connectivity', 'mapping_provenance')
    if not all(output.get(k) is not None for k in required) or not output['mapping_provenance']:
        raise MappingUnavailable('missing materialized records or original-ID mapping provenance')
    node_ids, elem_ids = output['node_ids'], output['element_ids']
    points, elements = output['coordinates'], output['connectivity']
    if output.get('n_nodes') is None or output.get('n_elems') is None:
        raise MappingUnavailable('missing output counts')
    if any(type(output[key]) is not int for key in ('n_nodes', 'n_elems')):
        raise ValueError('output counts must be integer counts')
    if not (len(points) == len(node_ids) == output['n_nodes'] > 0 and len(elements) == len(elem_ids) == output['n_elems'] > 0):
        raise ValueError('empty/count mismatch')
    for ids in (node_ids, elem_ids):
        if any(type(i) is not int or i <= 0 for i in ids) or len(set(ids)) != len(ids):
            raise ValueError('invalid/duplicate original IDs')
    nodes = []
    for ident, xyz in zip(node_ids, points):
        if len(xyz) != 3 or any(not math.isfinite(float(v)) for v in xyz):
            raise ValueError('invalid binary64 coordinates/depth')
        nodes.append((ident, ' '.join([str(ident), *(repr(float(v)) for v in xyz)])))
    indexing = output.get('connectivity_indexing', 'original-id')
    mapping = output.get('index_to_node_id')
    if indexing == 'row-offset' and (mapping is None or list(mapping) != list(node_ids)):
        raise MappingUnavailable('missing/ambiguous row-offset to original-ID map')
    if indexing not in ('row-offset', 'original-id'):
        raise MappingUnavailable('unknown connectivity indexing')
    valid_ids, elems = set(node_ids), []
    for ident, row in zip(elem_ids, elements):
        if indexing == 'row-offset':
            if any(type(i) is not int or i < 0 or i >= len(mapping) for i in row):
                raise ValueError('invalid storage offset')
            row = [mapping[i] for i in row]
        if len(row) not in (3, 4) or len(set(row)) != len(row) or any(type(i) is not int or i not in valid_ids for i in row):
            raise ValueError('invalid connectivity/arity/winding')
        elems.append((ident, ' '.join(map(str, [ident, len(row), *row]))))
    def digest(rows):
        return hashlib.sha256('\n'.join(line for _, line in sorted(rows)).encode()).hexdigest()
    return dict(n_nodes=len(nodes), n_elems=len(elems), node_sha256=digest(nodes), element_sha256=digest(elems))



def protocol_from_environment():
    """Open inherited channels; all three must be supplied together.

    No-channel use is a diagnostic invocation and cannot establish measurement
    validity. US2's controller supplies these file descriptors for campaigns.
    """
    keys = ('CHILMESH_CONTROL_IN_FD', 'CHILMESH_CONTROL_OUT_FD', 'CHILMESH_RESULT_FD')
    values = [os.environ.get(k) for k in keys]
    if not any(values):
        return None
    if not all(values) or not os.environ.get('CHILMESH_REPEAT_ID'):
        raise ProtocolError('incomplete inherited control channels')
    return ChildProtocol(os.environ['CHILMESH_REPEAT_ID'],
                         *(os.fdopen(int(value), mode, buffering=1) for value, mode in zip(values, ('r', 'w', 'w'))))


def warm_netcdf_reader(netCDF4, open_dataset, element_variable='connectivity', vertex_dimension='max_face_nodes'):
    """Open and load a three-node NetCDF file before READY.

    xarray loads optional packages (dask, numbagg, scipy, backend plugins) on the first
    open, not at import. This moves those imports outside the timed region for every
    NetCDF reader alike. The file is a temporary one-triangle mesh, never the measured
    mesh, so the measured file stays cold.
    """
    import tempfile  # netCDF4 is passed in by the adapter: this module imports only the standard library
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, 'warm.nc')
        with netCDF4.Dataset(path, 'w', format='NETCDF4') as ds:
            ds.createDimension('node', 3); ds.createDimension('nele', 1); ds.createDimension(vertex_dimension, 3)
            for name in ('x', 'y', 'depth'):
                ds.createVariable(name, 'f8', ('node',), fill_value=False)[:] = [0., 1., 2.]
            ds.createVariable('node_id', 'i8', ('node',), fill_value=False)[:] = [1, 2, 3]
            ds.createVariable('element_id', 'i8', ('nele',), fill_value=False)[:] = [1]
            ds.createVariable('arity', 'i8', ('nele',), fill_value=False)[:] = [3]
            ds.createVariable(element_variable, 'i8', ('nele', vertex_dimension), fill_value=False)[:] = [[1, 2, 3]]
        warm = open_dataset(path)
        warm.load()
        warm.close()


def eager_dataset(dataset):
    """Force a lazy dataset's full data into memory within the operation."""
    dataset.load()
    return dataset


def dataset_output(dataset, element_dimension):
    """Extract only D7's explicit original-ID variable contract, after END ACK."""
    sizes = dict(dataset.sizes)
    out = dict(n_nodes=sizes.get('node'), n_elems=sizes.get(element_dimension))
    try:
        def values(name):
            return dataset[name].values.tolist()
        ids, elements = values('node_id'), values('element_id')
        x, y, depth = values('x'), values('y'), values('depth')
        arity, connectivity = values('arity'), values('connectivity')
        if not (len(x) == len(y) == len(depth) == len(ids)):
            raise ValueError('coordinate lengths differ')
        if len(arity) != len(elements) or len(connectivity) != len(elements):
            raise ValueError('element lengths differ')
        if any(type(a) is not int or a not in (3, 4) or len(row) < a or any(v != 0 for v in row[a:]) for a, row in zip(arity, connectivity)):
            raise ValueError('invalid D7 arity/fill')
        out.update(node_ids=ids, element_ids=elements, coordinates=[list(p) for p in zip(x, y, depth)],
                   connectivity=[row[:a] for row, a in zip(connectivity, arity)],
                   connectivity_indexing='original-id',
                   mapping_provenance={'kind': 'fair-netcdf-v1 explicit node_id/element_id', 'element_dimension': element_dimension})
        out.update(canonical_output(out), verification_state='materialized')
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        out.update(verification_state='incomplete', reason=f'original-ID/content mapping unavailable: {exc}')
    return out


def thalassa_output(dataset):
    """Extract from Thalassa's NORMALIZED ADCIRC dataset (dims node, face, max_no_vertices)."""
    sizes = dict(dataset.sizes)
    out = dict(n_nodes=sizes.get('node'), n_elems=sizes.get('face'))
    try:
        def values(name):
            return dataset[name].values.tolist()
        ids, elements = values('node_id'), values('element_id')
        x, y, depth = values('lon'), values('lat'), values('depth')
        arity, faces = values('arity'), values('face_nodes')
        if not (len(x) == len(y) == len(depth) == len(ids)):
            raise ValueError('coordinate lengths differ')
        if len(arity) != len(elements) or len(faces) != len(elements):
            raise ValueError('element lengths differ')
        if any(type(v) is not int for v in ids) or any(type(v) is not int for v in elements):
            raise ValueError('non-integer ID')
        if any(type(a) is not int or a not in (3, 4) or len(row) < a or any(type(v) is not int for v in row)
               or any(v != -1 for v in row[a:]) for a, row in zip(arity, faces)):
            raise ValueError('invalid arity/padding')
        out.update(node_ids=ids, element_ids=elements, coordinates=[list(p) for p in zip(x, y, depth)],
                   connectivity=[row[:a] for row, a in zip(faces, arity)],
                   connectivity_indexing='row-offset', index_to_node_id=list(ids),
                   mapping_provenance={'kind': 'thalassa ADCIRC normalization: face_nodes = element - 1; row -> node_id variable',
                                       'element_dimension': 'face'})
        out.update(canonical_output(out), verification_state='materialized')
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        out.update(verification_state='incomplete', reason=f'original-ID/content mapping unavailable: {exc}')
    return out


def adcircpy_output(mesh):
    """Extract original IDs from adcircpy 1.2.7 (API checked in the pinned environment).

    ``mesh.nodes`` is a frame indexed by file node ID with columns x, y, value_1;
    ``mesh.elements.elements`` is a frame indexed by file element ID with columns
    n, node_1..node_n holding file node IDs. adcircpy stores value_1 = -depth, an
    exact sign flip, so depth is recovered bit for bit.
    """
    out = dict(n_nodes=None, n_elems=None)
    try:
        nodes, elems = mesh.nodes, mesh.elements.elements
        out.update(n_nodes=int(len(nodes)), n_elems=int(len(elems)))
        x, y, value = nodes['x'].tolist(), nodes['y'].tolist(), nodes['value_1'].tolist()
        arity = [int(n) for n in elems['n'].tolist()]
        cols = [elems[f'node_{k}'].tolist() for k in range(1, max(arity) + 1)]
        out.update(node_ids=[int(i) for i in nodes.index.tolist()], element_ids=[int(i) for i in elems.index.tolist()],
                   coordinates=[[a, b, -c] for a, b, c in zip(x, y, value)],
                   connectivity=[[int(cols[k][r]) for k in range(n)] for r, n in enumerate(arity)],
                   connectivity_indexing='original-id',
                   mapping_provenance={'kind': 'adcircpy 1.2.7 frame index = file ID; connectivity = file node IDs',
                                       'depth': 'negated value_1 (exact sign flip)'})
        out.update(canonical_output(out), verification_state='materialized')
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        out.update(verification_state='incomplete', reason=f'original-ID/content mapping unavailable: {exc}')
    return out


def run_adapter(tool, operation, *, materialize=lambda value: value, extract):
    """Metadata before READY; materialize inside region; hashes after END ACK."""
    version = importlib.metadata.version(tool)
    protocol = protocol_from_environment()
    def verify(output):
        content = extract(output)
        return dict(tool=tool, version=version, input_format='fort.14' if tool == 'adcircpy' else 'NetCDF',
                    materialization=True, output=content, n_nodes=content.get('n_nodes'), n_elems=content.get('n_elems'))
    if protocol is not None:
        protocol.ready()
        protocol.run(operation, materialize=materialize, verify=verify)
    else:
        modules = set(sys.modules)
        start = time.perf_counter()
        output = materialize(operation())
        end = time.perf_counter()
        new_modules = sorted(set(sys.modules) - modules)
        data = verify(output)
        data.update(parse_s=end - start, new_modules=new_modules, protocol_missing=True)
        print(json.dumps(data, allow_nan=False))
    return 0


def child_runtime_identity():
    """Collect build/digest/thread evidence before READY, without tool imports."""
    path = os.path.realpath(sys.executable)
    with open(path, 'rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return dict(version=sys.version, implementation=sys.implementation.name,
                implementation_fields={"name": sys.implementation.name, "cache_tag": sys.implementation.cache_tag,
                                       "version": list(sys.implementation.version), "hexversion": sys.hexversion},
                thread_environment={key: os.environ.get(key) for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')},
                platform=sysconfig.get_platform(), architecture=platform.machine(), executable=sys.executable, realpath=path, sha256=digest)

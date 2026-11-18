"""Independent fair-netcdf-v1 writer. Uses validated reference artifacts only.

Invoke outside repeat measurement. The written file is reopened with netCDF4;
content, IDs, arity and fill are verified before lineage can be accepted.
"""
import argparse
import bisect
import hashlib
import json
from pathlib import Path

REVISION = 'fair-netcdf-v1'


def _hash(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def _rows(reference, kind):
    path = Path(reference[kind + '_path'])
    expected = reference['canonical_' + ('node' if kind == 'nodes' else 'element') + '_sha256']
    if _hash(path) != expected:
        raise ValueError('reference artifact hash mismatch')
    with path.open() as stream:
        for line in stream:
            yield line.split()


def _row_index(ids, node_id):
    position = bisect.bisect_left(ids, node_id)
    if position == len(ids) or ids[position] != node_id:
        raise ValueError('element references unknown node ID')
    return position


def convert(reference, reader, output, *, netcdf=None):
    """Write and independently reopen the declared D7 content contract."""
    if netcdf is None:
        import netCDF4 as netcdf
    if reader not in ('xarray', 'thalassa') or reference['verification_state'] != 'verified':
        raise ValueError('unsupported reader or incomplete reference')
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    nn, ne = reference['header_node_count'], reference['header_element_count']
    max_arity = max(int(row[1]) for row in _rows(reference, 'elements'))
    base = {'node_id': 'i8', 'x': 'f8', 'y': 'f8', 'depth': 'f8', 'element_id': 'i8', 'arity': 'i8'}
    if reader == 'xarray':
        dimensions, indexing, conn = {'node': nn, 'nele': ne, 'max_face_nodes': max_arity}, 'original-file-ids', 'connectivity'
    else:
        dimensions, indexing, conn = {'node': nn, 'nele': ne, 'nvertex': max_arity}, 'one-based-row-index', 'element'
    pad_dim = 'max_face_nodes' if reader == 'xarray' else 'nvertex'
    # fill_value 0 is the literal pad value; no _FillValue attribute is declared.
    contract = dict(revision=REVISION, dimensions=dimensions, variables={**base, conn: 'i8'}, fill_value=0,
                    indexing=indexing, ordering='ascending original IDs')
    dimension = 'nele'
    ids = []  # the single node-id array (ascending); row index = position
    with netcdf.Dataset(str(output), 'w', format='NETCDF4') as ds:
        for name, size in contract['dimensions'].items(): ds.createDimension(name, size)
        for name, dtype in contract['variables'].items():
            dims = ('node',) if name in ('node_id', 'x', 'y', 'depth') else (dimension, pad_dim) if name == conn else (dimension,)
            ds.createVariable(name, dtype, dims, fill_value=False)
        for index, row in enumerate(_rows(reference, 'nodes')):
            ds['node_id'][index] = int(row[0]); ids.append(int(row[0]))
            for name, value in zip(('x', 'y', 'depth'), row[1:]): ds[name][index] = float(value)
        for index, row in enumerate(_rows(reference, 'elements')):
            ident, arity, *nodes = map(int, row)
            ds['element_id'][index], ds['arity'][index] = ident, arity
            if reader == 'thalassa':
                nodes = [_row_index(ids, v) + 1 for v in nodes]
            ds[conn][index, :] = nodes + [0] * (max_arity - arity)
    hashes = {}
    # Stream verification so STOFS does not require another mesh-sized object.
    with netcdf.Dataset(str(output), 'r') as ds:
        ds.set_auto_mask(False)
        if len(ds.dimensions['node']) != nn or len(ds.dimensions[dimension]) != ne:
            raise ValueError('converted dimension mismatch')
        node_ids = [int(v) for v in ds['node_id'][:]] if reader == 'thalassa' else None
        for name, expected_dtype in contract['variables'].items():
            if ds[name].dtype.kind != expected_dtype[0] or ds[name].dtype.itemsize != 8:
                raise ValueError('converted dtype mismatch')
        for kind, expected_count in (('nodes', nn), ('elements', ne)):
            digest, count = hashlib.sha256(), 0
            for index, source in enumerate(_rows(reference, kind)):
                if kind == 'nodes':
                    row = [str(int(ds['node_id'][index])), *(repr(float(ds[k][index])) for k in ('x', 'y', 'depth'))]
                else:
                    arity = int(ds['arity'][index]); cells = ds[conn][index].tolist()
                    if arity not in (3, 4) or any(v != 0 for v in cells[arity:]):
                        raise ValueError('converted arity/fill mismatch')
                    cells = cells[:arity]
                    if node_ids is not None:
                        if any(v < 1 or v > len(node_ids) for v in cells):
                            raise ValueError('converted element index out of range')
                        cells = [node_ids[v - 1] for v in cells]
                    row = list(map(str, [int(ds['element_id'][index]), arity, *cells]))
                if row != source: raise ValueError('converted content/ID mapping mismatch')
                digest.update((('\n' if count else '') + ' '.join(row)).encode()); count += 1
            if count != expected_count: raise ValueError('converted count mismatch')
            hashes['node_sha256' if kind == 'nodes' else 'element_sha256'] = digest.hexdigest()
    expected = {'node_sha256': reference['canonical_node_sha256'], 'element_sha256': reference['canonical_element_sha256']}
    if hashes != expected: raise ValueError('converted canonical hash mismatch')
    return dict(source_reference_id=reference['reference_id'], source_sha256=reference['raw_file_sha256'],
                converted_sha256=_hash(output), converted_path=str(output.resolve()), converted_format='NetCDF',
                reader=reader, contract=contract, converter_revision=REVISION,
                converter_sha256=_hash(__file__), netcdf4_version=netcdf.__version__,
                native_versions={'netcdf': netcdf.__netcdf4libversion__, 'hdf5': netcdf.__hdf5libversion__},
                source_to_row_maps={'nodes_path': reference['nodes_path'], 'elements_path': reference['elements_path'],
                                    'rule': 'row is ascending original ID rank', 'hashes': expected},
                source_to_converted_id_map={'rule': 'identity original IDs' if reader == 'xarray' else 'one-based row index via node_id variable', 'independently_reopened': True},
                source_content_verification=dict(verification_state='verified', n_nodes=nn, n_elems=ne, **hashes))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', required=True)
    parser.add_argument('--reader', required=True, choices=['xarray', 'thalassa'])
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    reference = json.loads(Path(args.reference).read_text())
    print(json.dumps(convert(reference, args.reader, args.output), allow_nan=False))


if __name__ == '__main__':
    main()

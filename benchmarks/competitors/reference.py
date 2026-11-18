"""Independent fort.14 content reference, bounded sorting and disk-backed ID checks.

Only the declared nodes/elements are canonicalized. Boundary sections remain
opaque but are included in the digest of the very same opened input stream.
"""
import hashlib
import heapq
import math
import re
import sqlite3
import tempfile
from contextlib import ExitStack
from pathlib import Path

CANONICALIZATION_REVISION = 'binary64-repr-v1'
PARSER_REVISION = 'fort14-reference-v1'


class ReferenceError(ValueError):
    """Source bytes cannot establish an independent content reference."""


def significant_digits(token):
    mantissa = re.split('[eE]', token.lstrip('+-'))[0].replace('.', '')
    return len(mantissa.lstrip('0'))


def _positive(token):
    if not re.fullmatch(r'[0-9]+', token) or int(token) <= 0:
        raise ReferenceError('IDs, counts and arities must be positive integers')
    return int(token)


def _number(token):
    if not re.fullmatch(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?', token):
        raise ReferenceError('malformed/nonfinite numeric token')
    value = float(token)
    if not math.isfinite(value):
        raise ReferenceError('nonfinite numeric token')
    return repr(value)


def _sort_chunks(lines, directory, prefix, size):
    chunks, batch = [], []
    for line in lines:
        batch.append(line)
        if len(batch) == size:
            chunks.append(_chunk(batch, directory, prefix, len(chunks)))
            batch = []
    if batch:
        chunks.append(_chunk(batch, directory, prefix, len(chunks)))
    return chunks


def _chunk(batch, directory, prefix, index):
    path = Path(directory) / f'{prefix}-{index}'
    path.write_text('\n'.join(sorted(batch, key=lambda line: int(line.split(' ', 1)[0]))) + '\n')
    return path


def _merge_group(paths, output):
    with ExitStack() as stack:
        files = [stack.enter_context(p.open()) for p in paths]
        stream = stack.enter_context(output.open('w'))
        for row in heapq.merge(*files, key=lambda line: int(line.split(' ', 1)[0])):
            stream.write(row)


def _canonical_file(chunks, output, directory):
    # Bound open file descriptors as well as the in-memory record chunk.
    generation = 0
    while len(chunks) > 32:
        merged = []
        for i in range(0, len(chunks), 32):
            target = Path(directory) / f'merge-{output.name}-{generation}-{i}'
            _merge_group(chunks[i:i + 32], target)
            merged.append(target)
        chunks = merged
        generation += 1
    digest, previous, count = hashlib.sha256(), None, 0
    with ExitStack() as stack:
        streams = [stack.enter_context(p.open()) for p in chunks]
        out = stack.enter_context(output.open('xb'))
        for line in heapq.merge(*streams, key=lambda row: int(row.split(' ', 1)[0])):
            ident = int(line.split(' ', 1)[0])
            if previous == ident:
                raise ReferenceError('duplicate ID')
            raw = (('\n' if count else '') + line.rstrip('\n')).encode('utf-8')
            digest.update(raw); out.write(raw)
            previous, count = ident, count + 1
    return digest.hexdigest(), count


def build_reference(path, output_dir, *, name=None, version=None, served_sha256=None,
                    valence_asset=None, chunk_size=100000):
    """Bind header and validated canonical content to one raw-byte SHA-256.

    Canonical artifacts are persistent and contain every validated record for
    independent conversion. Caller owns output_dir; existing files never change.
    """
    if type(chunk_size) is not int or chunk_size < 1:
        raise ReferenceError('chunk_size must be positive')
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    precision, total = 0, 0
    sequential = {'nodes': True, 'elements': True}
    digest = hashlib.sha256()
    with tempfile.TemporaryDirectory(dir=output_dir) as scratch, Path(path).open('rb') as source:
        db = sqlite3.connect(str(Path(scratch) / 'ids.sqlite'))
        try:
            db.execute('CREATE TABLE nodes (id INTEGER PRIMARY KEY)')
            def read():
                nonlocal total
                raw = source.readline(); digest.update(raw); total += len(raw)
                return raw.decode('utf-8').split()
            if not read():
                raise ReferenceError('missing title')
            header = read()
            if len(header) != 2:
                raise ReferenceError('invalid header')
            ne, nn = map(_positive, header)
            def nodes():
                nonlocal precision
                for _ in range(nn):
                    fields = read()
                    if len(fields) != 4:
                        raise ReferenceError('missing/malformed node record')
                    ident = _positive(fields[0])
                    if ident != _ + 1: sequential['nodes'] = False
                    db.execute('INSERT INTO nodes VALUES (?)', (ident,))
                    values = [_number(t) for t in fields[1:]]
                    precision += sum(significant_digits(t) > 17 for t in fields[1:])
                    yield ' '.join([str(ident), *values])
            node_chunks = _sort_chunks(nodes(), scratch, 'node', chunk_size)
            def elements():
                for _ in range(ne):
                    fields = read()
                    if len(fields) < 2:
                        raise ReferenceError('missing element record')
                    ident, arity = map(_positive, fields[:2])
                    if ident != _ + 1: sequential['elements'] = False
                    if arity not in (3, 4) or len(fields) != arity + 2:
                        raise ReferenceError('unsupported arity/malformed element')
                    ids = [_positive(t) for t in fields[2:]]
                    if len(set(ids)) != arity:
                        raise ReferenceError('duplicate vertex in element')
                    if any(db.execute('SELECT 1 FROM nodes WHERE id=?', (i,)).fetchone() is None for i in ids):
                        raise ReferenceError('missing node ID in connectivity')
                    yield ' '.join(map(str, [ident, arity, *ids]))
            element_chunks = _sort_chunks(elements(), scratch, 'element', chunk_size)
            boundary_digest = hashlib.sha256()
            while True:
                raw = source.read(1 << 20)
                if not raw: break
                digest.update(raw); boundary_digest.update(raw); total += len(raw)
            raw_hash = digest.hexdigest()
            if served_sha256 is not None and served_sha256 != raw_hash:
                raise ReferenceError('served digest mismatch')
            nodes_path = output_dir / f'{raw_hash}.nodes'
            elements_path = output_dir / f'{raw_hash}.elements'
            nh, nc = _canonical_file(node_chunks, nodes_path, scratch)
            eh, ec = _canonical_file(element_chunks, elements_path, scratch)
        except (UnicodeError, sqlite3.IntegrityError, OverflowError) as exc:
            raise ReferenceError(f'invalid source/duplicate ID: {exc}') from exc
        finally:
            db.close()
    return dict(reference_id=raw_hash, name=name, exact_version=version, aliases=[],
                valence_asset=valence_asset, provenance='Valence registry' if valence_asset else 'caller supplied',
                served_sha256=served_sha256, raw_file_sha256=raw_hash, raw_size_bytes=total,
                header_node_count=nn, header_element_count=ne, n_nodes=nc, n_elems=ec,
                canonical_node_sha256=nh, canonical_element_sha256=eh,
                canonicalization_revision=CANONICALIZATION_REVISION, parser_revision=PARSER_REVISION,
                precision_token_count_gt17=precision, malformed=[], nonfinite=[], duplicate=[], missing=[],
                verification_state='incomplete' if precision else 'verified',
                nodes_path=str(nodes_path.resolve()), elements_path=str(elements_path.resolve()),
                evidence_hashes={'nodes': nh, 'elements': eh}, rights_ruling_id=None,
                boundary_section_sha256=boundary_digest.hexdigest(),
                file_order_ids_sequential=dict(sequential),
                sort_chunks=len(node_chunks) + len(element_chunks))

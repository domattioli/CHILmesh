"""Verify materialized tool output; never open the measured input."""
from .reference import CANONICALIZATION_REVISION


from .protocol import MappingUnavailable, canonical_output


def row_order_mapping(output, reference):
    """Map row order to original IDs only when the reference proves IDs are 1..N in file order.

    Applies to an output that declares record_order='row-order' and carries no IDs.
    The canonical hashes still have to match, so a reader that reorders or drops rows fails.
    """
    if output.get('record_order') != 'row-order' or output.get('node_ids') is not None or output.get('element_ids') is not None:
        return output
    proof = reference.get('file_order_ids_sequential') or {}
    if proof.get('nodes') is not True or proof.get('elements') is not True:
        return output
    n, m = output.get('n_nodes'), output.get('n_elems')
    if type(n) is not int or type(m) is not int:
        return output
    ids = list(range(1, n + 1))
    return dict(output, node_ids=ids, element_ids=list(range(1, m + 1)), index_to_node_id=ids,
                mapping_provenance={'kind': 'row order; reference proves file IDs are 1..N in file order',
                                    'reference_id': reference.get('reference_id'), 'triangle_padding': 'repeated first vertex'})


def verify_output(output, reference):
    expected = dict(n_nodes=reference['header_node_count'], n_elems=reference['header_element_count'],
                    node_sha256=reference['canonical_node_sha256'], element_sha256=reference['canonical_element_sha256'])
    evidence = dict(n_nodes=output.get('n_nodes'), n_elems=output.get('n_elems'), node_sha256=None,
                    element_sha256=None, expected=expected, mapping_provenance=output.get('mapping_provenance'),
                    normalization_revision=CANONICALIZATION_REVISION, reasons=[], missing={})
    try:
        output = row_order_mapping(output, reference)
        evidence['mapping_provenance'] = output.get('mapping_provenance')
        evidence.update(canonical_output(output))
        evidence['reasons'] = [f'{k} mismatch' for k, value in expected.items() if evidence[k] != value]
        state = 'failed' if evidence['reasons'] else 'verified'
        if reference['verification_state'] != 'verified':
            state = 'incomplete'; evidence['reasons'].append('reference precision/content incomplete')
    except MappingUnavailable as exc:
        state = 'incomplete'; evidence['reasons'].append(str(exc))
    except (ValueError, TypeError, KeyError, IndexError, OverflowError) as exc:
        state = 'failed'; evidence['reasons'].append(str(exc))
    evidence['verification_state'] = state
    evidence['missing'] = {k: '; '.join(evidence['reasons']) for k, v in evidence.items() if v is None}
    return evidence


def verify_written_output(path, reference):
    """Read the actual written artifact with the independent parser after END.

    Boundary bytes are conservatively required to match too. Formatting changes
    in boundaries remain incomplete until an independent semantic boundary
    contract exists; matching nodes alone never certifies a full roundtrip.
    """
    import tempfile
    from pathlib import Path
    from .reference import build_reference, ReferenceError
    if path is None or not Path(path).is_file():
        return {'verification_state': 'incomplete', 'reason': 'written artifact unavailable'}
    try:
        with tempfile.TemporaryDirectory() as scratch:
            actual = build_reference(path, scratch)
            keys = ('header_node_count', 'header_element_count', 'canonical_node_sha256', 'canonical_element_sha256')
            mismatches = [key for key in keys if actual[key] != reference[key]]
            state = 'failed' if mismatches else 'verified'
            if actual['verification_state'] != 'verified' or reference['verification_state'] != 'verified':
                state = 'incomplete'; mismatches.append('precision incomplete')
            if actual.get('boundary_section_sha256') != reference.get('boundary_section_sha256'):
                state = 'incomplete'; mismatches.append('boundary bytes differ; semantic boundary proof unavailable')
            return dict(verification_state=state, reason='; '.join(mismatches), path=str(path),
                        sha256=actual['raw_file_sha256'], n_nodes=actual['n_nodes'], n_elems=actual['n_elems'],
                        node_sha256=actual['canonical_node_sha256'], element_sha256=actual['canonical_element_sha256'])
    except (OSError, ValueError, ReferenceError) as exc:
        return {'verification_state': 'failed', 'reason': str(exc), 'path': str(path)}

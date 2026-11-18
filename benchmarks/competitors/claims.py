"""Deterministic claim gates over recomputed, rights-cleared series cells.

Only the bounded 'less' template is certified. Arbitrary prose is never echoed
as an accepted sentence. Cell summaries come from the offline report, not from
user-supplied claim arithmetic.
"""
import math
from .aggregate import METRICS
from .campaign import _hash, _pin


def evaluate_claim(proposal, cells):
    failed = []
    def check(condition, predicate):
        if not condition and predicate not in failed:
            failed.append(predicate)
    check(isinstance(proposal, dict), 'proposal')
    proposal = proposal if isinstance(proposal, dict) else {}
    allowed = {'claim_id', 'template', 'numerator_cell_id', 'denominator_cell_id', 'metric', 'unit',
               'direction', 'mesh', 'exact_version', 'operation', 'formula', 'resource_condition_reference',
               'denominator_metric'}
    check(proposal.get('template') == 'less' and not set(proposal) - allowed, 'template')
    check(bool(proposal.get('claim_id')), 'identity')
    check(proposal.get('formula') == '100*(B-A)/B', 'formula')
    check(proposal.get('direction') == 'lower', 'direction')
    metric = proposal.get('metric')
    check(metric in METRICS and proposal.get('unit') == METRICS.get(metric)
          and proposal.get('denominator_metric', metric) == metric, 'metric')
    a = cells.get(proposal.get('numerator_cell_id'))
    b = cells.get(proposal.get('denominator_cell_id'))
    check(bool(a) and bool(b) and a is not b, 'cells')
    value, hashes, sentence = None, [], None
    if a and b:
        pair = (a, b)
        check(all(c.get('equivalence') == 'like-for-like' for c in pair), 'equivalence')
        check(all(c.get('verification') == 'verified' for c in pair), 'verification')
        check(all(c.get('complete') is True for c in pair), 'complete')
        check(all(c.get('stable') is True and all(c.get('metrics', {}).get(k, {}).get('unstable') is False
                  for k in METRICS) for c in pair), 'stable')
        check(all(c.get('conditions_valid') is True for c in pair), 'conditions')
        check(all(c.get('rights_cleared') is True for c in pair), 'rights')
        for predicate, keys in [('input', ('input_format', 'input_sha256', 'reference_id')),
                                ('contract', ('contract_id', 'backend', 'options')),
                                ('resources', ('manifest_sha256', 'series_id', 'host_id', 'host_allocation', 'thread_environment'))]:
            check(all(a.get(k) is not None and a.get(k) == b.get(k) for k in keys), predicate)
        # Converted NetCDF reads never establish native fort.14 advantages.
        check(all(c.get('input_format') == 'fort.14' and not c.get('conversion') for c in pair), 'input')
        check(all(c.get('mesh') == proposal.get('mesh') and c.get('exact_version') == proposal.get('exact_version')
                  and c.get('mesh') and c.get('exact_version') for c in pair), 'mesh')
        check(all(c.get('operation') == proposal.get('operation') and c.get('operation') for c in pair), 'operation')
        check(bool(proposal.get('resource_condition_reference')) and all(
            c.get('resource_condition_reference') == proposal.get('resource_condition_reference') for c in pair), 'resources')
        check(all(_pin(c.get('tool_version')) for c in pair), 'versions')
        check(all(c.get('source_record_hashes') and all(_hash(h) for h in c['source_record_hashes']) and c.get('tool_version') and c.get('attempt_ids')
                  and c.get('series_id') for c in pair), 'evidence')
        summaries = [c.get('metrics', {}).get(metric, {}) for c in pair]
        check(all(s.get('unit') == proposal.get('unit') and s.get('n', 0) >= c.get('required_repeat_count', 1)
                  and s.get('raw_values') and s.get('relative_spread_pct') is not None
                  for c, s in zip(pair, summaries)), 'metric')
        av, bv = (s.get('median') for s in summaries)
        valid_numbers = all(type(v) in (float, int) and math.isfinite(v) and v >= 0 for v in (av, bv))
        check(valid_numbers and bv > 0, 'denominator')
        if valid_numbers and bv > 0:
            value = 100 * (bv - av) / bv
            check(value >= 0, 'direction')
        hashes = sorted({h for c in pair for h in c.get('source_record_hashes', [])})
        if not failed:
            details = '; '.join(f"{c['tool']} {c['tool_version']}: cell {c['cell_id']}, series {c['series_id']}, "
                f"n={s['n']}, median={s['median']}, spread={s['relative_spread_pct']}%, "
                f"repeats={','.join(c['attempt_ids'])}" for c, s in zip(pair, summaries))
            sentence = (f"On {a['mesh']} {a['exact_version']}, {a['operation']}, {a['tool']} uses "
                f"{value:.12g}% less {metric} ({proposal['unit']}) than {b['tool']}, with {b['tool']} "
                f"as denominator; {details}; scope {proposal['resource_condition_reference']}; "
                f"record sha256 {','.join(hashes)}.")
    return dict(claim_id=proposal.get('claim_id'), accepted=not failed, failed_predicates=failed,
                value=value, sentence=sentence, source_record_hashes=hashes)

"""Deterministic, descriptive summaries of one immutable attempt series."""
import math
import statistics
from .campaign import REPEATS, validate_replacement

METRICS = {'operation_elapsed_s': 'seconds', 'rss_whole_process_peak_bytes': 'bytes', 'rss_growth_bytes': 'bytes'}
IDENTITY = ('campaign_id', 'manifest_sha256', 'cell_id', 'series_id', 'reference_id', 'tool', 'environment_id',
            'tool_version', 'contract_id', 'backend', 'host_id', 'host_allocation', 'thread_environment', 'options',
            'input_sha256', 'input_format')


def summarize_metric(pairs, *, unit):
    pairs = list(pairs)
    if any(type(v) not in (int, float) or not math.isfinite(v) or v < 0 for _, v in pairs):
        raise ValueError('metric values must be finite nonnegative numbers')
    values = [v for _, v in pairs]
    median = statistics.median(values) if values else None
    low, high = (min(values), max(values)) if values else (None, None)
    spread = high-low if values else None
    infinite = median == 0 and spread > 0
    relative = None if not values or infinite else (100*spread/median if median else 0)
    return dict(unit=unit, raw_values=[dict(repeat_id=i, value=v) for i,v in pairs], n=len(values),
                median=median, min=low, max=high, absolute_spread=spread, relative_spread_pct=relative,
                spread_state='infinite' if infinite else ('finite' if values else 'missing'),
                unstable=infinite or relative is None or relative > 10)


def aggregate_series(attempts, *, mesh_class, previous_series_id=None, replacement_rationale=None):
    attempts = list(attempts)
    if not attempts: raise ValueError('empty series')
    if previous_series_id is not None or replacement_rationale is not None:
        validate_replacement(previous_series_id, replacement_rationale)
        if previous_series_id == attempts[0]['series_id']: raise ValueError('cannot replace self')
    if any(any(r.get(k) != attempts[0].get(k) for k in IDENTITY) for r in attempts):
        raise ValueError('cannot pool different cell identities')
    def conversion_identity(r):
        conversion=r.get('conversion')
        if not isinstance(conversion, dict): return conversion
        return {k:v for k,v in conversion.items() if k not in ('read_s','end_to_end_s')}
    if any(conversion_identity(r) != conversion_identity(attempts[0]) for r in attempts):
        raise ValueError('cannot pool different conversion input provenance')
    actual = (attempts[0].get('adapter_result') or {}).get('actual_backend')
    if any((r.get('adapter_result') or {}).get('actual_backend') != actual for r in attempts):
        raise ValueError('cannot pool different actual backend/work evidence')
    ids = [r['repeat_id'] for r in attempts]
    indices = [r['repeat_index'] for r in attempts]
    if len(set(ids)) != len(ids) or len(set(indices)) != len(indices): raise ValueError('duplicate attempt')
    included, excluded = [], []
    for r in attempts:
        valid = r['outcome'] == 'ok' and r.get('verification') == 'verified' and all(r.get('measurement_validity', {}).get(k, {}).get('valid') is True for k in ('time','memory','cache','load','runtime'))
        valid = valid and all(type(r['measurement'].get(k)) in (int,float) and math.isfinite(r['measurement'][k]) and r['measurement'][k] >= 0 for k in METRICS)
        if valid: included.append(r)
        else: excluded.append(dict(repeat_id=r['repeat_id'], reason=r.get('reason') or 'invalid or missing measurement evidence'))
    metrics = {k: summarize_metric([(r['repeat_id'],r['measurement'][k]) for r in included],unit=unit) for k,unit in METRICS.items()}
    required = REPEATS[mesh_class]
    complete = not excluded and len(included) >= required and sorted(indices) == list(range(len(indices)))
    return dict(series_id=attempts[0]['series_id'], cell_id=attempts[0]['cell_id'], required_repeat_count=required,
                attempt_ids=ids, included_attempt_ids=[r['repeat_id'] for r in included], excluded_attempts=excluded,
                complete=complete, stable=all(not s['unstable'] for s in metrics.values()), metrics=metrics,
                eligible=complete and all(not s['unstable'] for s in metrics.values()),
                previous_series_id=previous_series_id, replacement_rationale=replacement_rationale)

"""Evidence for a Linux cold reset, performed only on explicit controller calls."""
import os
import re
import time
from datetime import datetime, timezone
from pathlib import Path


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def _drop(value):
    Path('/proc/sys/vm/drop_caches').write_text(value)
    return 0


def reset_cache(*, sync=None, drop=None, read_meminfo=None, clock=time.perf_counter):
    """No input file access. Inject operations for unprivileged offline testing."""
    sync = sync or (lambda: (os.sync(), 0)[1])
    drop = drop or _drop
    read_meminfo = read_meminfo or (lambda: Path('/proc/meminfo').read_text())
    result = dict(operations=[], cached_before_kib=None, cached_after_kib=None,
                  cached_before_bytes=None, cached_after_bytes=None, reset_utc=None, reset_s=None,
                  launch_s=None, success=False, cache='unknown', reason='', missing={})
    def operation(name, fn):
        entry = dict(operation=name, monotonic_s=clock(), utc=utc_now(), exit_code=None)
        if name == 'drop_caches': entry['value'] = '3'
        result['operations'].append(entry)
        try: entry['exit_code'] = fn()
        except OSError as exc: entry['error'] = str(exc); raise
        if entry['exit_code'] != 0: raise ValueError(f'{name} exit status {entry["exit_code"]}')
    def cached():
        match = re.search(r'^Cached:\s+(\d+) kB\s*$', read_meminfo(), re.M)
        if match is None: raise ValueError('Cached unavailable')
        return int(match[1])
    try:
        operation('sync', sync)
        result['cached_before_kib'] = cached()
        operation('drop_caches', lambda: drop('3'))
        result['cached_after_kib'] = cached()
        if result['cached_after_kib'] >= result['cached_before_kib']: raise ValueError('Cached did not decrease')
        result.update(success=True, cache='cold')
    except (OSError, ValueError) as exc: result['reason'] = str(exc)
    result.update(reset_s=clock(), reset_utc=utc_now())
    for prefix in ('before','after'):
        value = result[f'cached_{prefix}_kib']
        result[f'cached_{prefix}_bytes'] = None if value is None else value*1024
    result['missing'] = {k: result['reason'] or 'child not launched yet' for k,v in result.items() if v is None}
    return result

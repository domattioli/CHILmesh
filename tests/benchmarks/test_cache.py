import pytest


@pytest.mark.parametrize('mode',['ok','sync-fail','permission','no-decrease','missing'])
def test_reset_evidence_order_and_failure(mode):
    from benchmarks.competitors.cache import reset_cache
    calls=[]; reads=iter(['Cached: 100 kB\n', 'Cached: %s kB\n'% ('100' if mode=='no-decrease' else '20')])
    def read():
        calls.append('meminfo')
        if mode=='missing': raise OSError('missing')
        return next(reads)
    def sync(): calls.append('sync'); return 1 if mode=='sync-fail' else 0
    def drop(value):
        calls.append(('drop',value))
        if mode=='permission': raise PermissionError('denied')
        return 0
    ticks=iter(range(20))
    r=reset_cache(sync=sync,drop=drop,read_meminfo=read,clock=lambda:next(ticks))
    assert r['success'] == (mode=='ok')
    assert r['cache'] == ('cold' if mode=='ok' else 'unknown')
    if mode=='ok':
        assert calls==['sync','meminfo',('drop','3'),'meminfo']
        assert r['cached_before_bytes']==102400 and r['cached_after_bytes']==20480
        assert [x['operation'] for x in r['operations']]==['sync','drop_caches']
        assert r['reset_s'] >= r['operations'][-1]['monotonic_s']
    else: assert r['reason']

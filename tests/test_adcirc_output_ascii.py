"""ASCII ADCIRC output readers (spec 005 US4): fort.61/62/63/64, maxele.63, maxvel.63."""
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent / "fixtures"))
import synthetic_adcirc as syn  # noqa: E402

from chilmesh import examples, read_adcirc_output  # noqa: E402

N = 7
TIMES = [10.0, 20.0, 30.0]


def _scalar(nt=3, n=N):
    return np.arange(nt * n, dtype=float).reshape(nt, n) / 10.0 - 0.5


def test_full_scalar_values(tmp_path):
    vals = _scalar()
    p = syn.write_output_ascii(tmp_path / "fort.63", TIMES, vals)
    s = read_adcirc_output(p)
    assert s.kind == "scalar" and s.location == "nodes" and s.dropped_records == 0
    np.testing.assert_allclose(s.times, TIMES)
    np.testing.assert_array_equal(s.iterations, [1, 2, 3])
    np.testing.assert_allclose(s.values, vals, atol=1e-9)
    assert s.values.shape == (3, N)


def test_full_vector_values(tmp_path):
    vals = np.arange(3 * N * 2, dtype=float).reshape(3, N, 2) / 7.0
    s = read_adcirc_output(syn.write_output_ascii(tmp_path / "fort.64", TIMES, vals))
    assert s.kind == "vector" and s.values.shape == (3, N, 2)
    np.testing.assert_allclose(s.values, vals, atol=1e-9)


def test_sparse_matches_full_with_nan(tmp_path):
    vals = _scalar()
    vals[0, 2] = np.nan
    vals[2, 5] = np.nan
    full = read_adcirc_output(syn.write_output_ascii(tmp_path / "a.63", TIMES, vals))
    sparse = read_adcirc_output(syn.write_output_ascii(tmp_path / "b.63", TIMES, vals, sparse=True))
    assert np.isnan(sparse.values[0, 2]) and np.isnan(sparse.values[2, 5])
    np.testing.assert_allclose(sparse.values, full.values, equal_nan=True)


def test_sparse_unlisted_nodes_default_to_nan(tmp_path):
    vals = np.full((1, N), np.nan)
    vals[0, 3] = 1.25
    s = read_adcirc_output(syn.write_output_ascii(tmp_path / "fort.63", [5.0], vals, sparse=True))
    assert s.values[0, 3] == pytest.approx(1.25)
    assert np.isnan(np.delete(s.values[0], 3)).all()


def test_fill_value_becomes_nan(tmp_path):
    vals = _scalar()
    vals[1, 0] = -99999.0
    s = read_adcirc_output(syn.write_output_ascii(tmp_path / "fort.63", TIMES, vals))
    assert np.isnan(s.values[1, 0])
    assert not np.isnan(s.values[0, 0])


def test_station_files(tmp_path):
    vals = _scalar(nt=4, n=3)
    s = read_adcirc_output(syn.write_output_ascii(tmp_path / "fort.61", [1, 2, 3, 4], vals))
    assert s.location == "stations" and s.values.shape == (4, 3)
    v = np.ones((2, 3, 2))
    s2 = read_adcirc_output(syn.write_output_ascii(tmp_path / "fort.62", [1, 2], v))
    assert s2.location == "stations" and s2.kind == "vector"


def test_node_count_mismatch_names_both_counts(tmp_path):
    mesh = examples.annulus.__wrapped__()
    p = syn.write_output_ascii(tmp_path / "fort.63", TIMES, _scalar())
    with pytest.raises(ValueError) as ei:
        read_adcirc_output(p, mesh=mesh)
    assert str(N) in str(ei.value) and str(mesh.n_verts) in str(ei.value)


def test_node_count_match_ok(tmp_path):
    mesh = examples.annulus.__wrapped__()
    vals = np.linspace(0, 1, mesh.n_verts)[None, :]
    s = read_adcirc_output(syn.write_output_ascii(tmp_path / "fort.63", [1.0], vals), mesh=mesh)
    assert s.values.shape == (1, mesh.n_verts)


def test_single_step(tmp_path):
    vals = _scalar()
    p = syn.write_output_ascii(tmp_path / "fort.63", TIMES, vals)
    for k in (0, 1, 2, -1):
        s = read_adcirc_output(p, step=k)
        assert s.values.shape == (1, N)
        np.testing.assert_allclose(s.values[0], vals[k], atol=1e-9)
        np.testing.assert_allclose(s.times, [TIMES[k]])
    sp = syn.write_output_ascii(tmp_path / "sp.63", TIMES, vals, sparse=True)
    np.testing.assert_allclose(read_adcirc_output(sp, step=1).values[0], vals[1], atol=1e-9)
    with pytest.raises(IndexError):
        read_adcirc_output(p, step=3)


def test_single_step_file_with_one_record(tmp_path):
    vals = _scalar(nt=1)
    s = read_adcirc_output(syn.write_output_ascii(tmp_path / "fort.63", [9.0], vals))
    assert s.values.shape == (1, N)


def test_truncated_warns_with_dropped_count(tmp_path):
    vals = _scalar(nt=4)
    p = syn.write_output_ascii(tmp_path / "fort.63", [1, 2, 3, 4], vals)
    syn.truncate_lines(p, 3)  # cuts into set 4
    with pytest.warns(UserWarning, match=r"1 record"):
        s = read_adcirc_output(p)
    assert s.dropped_records == 1 and s.values.shape == (3, N)
    np.testing.assert_allclose(s.values, vals[:3], atol=1e-9)


def test_truncated_multiple_sets_dropped(tmp_path):
    vals = _scalar(nt=4)
    p = syn.write_output_ascii(tmp_path / "fort.63", [1, 2, 3, 4], vals)
    syn.truncate_lines(p, N + 2)  # drops set 4 entirely and cuts set 3
    with pytest.warns(UserWarning, match=r"2 record"):
        s = read_adcirc_output(p)
    assert s.dropped_records == 2 and s.values.shape == (2, N)


def test_not_truncated_does_not_warn(tmp_path):
    p = syn.write_output_ascii(tmp_path / "fort.63", TIMES, _scalar())
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        read_adcirc_output(p)


@pytest.mark.parametrize("name", ["maxele.63", "maxvel.63"])
def test_max_files_second_record_is_time_of_max(tmp_path, name):
    mx = np.linspace(0.1, 2.0, N)
    tom = np.linspace(100.0, 900.0, N)
    s = read_adcirc_output(syn.write_max_ascii(tmp_path / name, mx, tom))
    assert s.values.shape == (2, N)
    np.testing.assert_allclose(s.values[0], mx, atol=1e-9)
    np.testing.assert_allclose(s.values[1], tom, atol=1e-6)


def test_unsupported_irtype_raises(tmp_path):
    p = tmp_path / "fort.63"
    p.write_text("h\n 1 2 1.0 1 3\n 1.0 1\n 1 0.5\n 2 0.5\n")
    with pytest.raises(ValueError, match="IRTYPE"):
        read_adcirc_output(p)


def test_single_step_memory_1m_nodes(tmp_path):
    """SC-005: peak allocation of a single-step read <= 3.0 x returned array nbytes."""
    import tracemalloc
    n = 1_000_000
    v = np.linspace(-1, 1, n)[None, :]
    p = syn.write_output_ascii(tmp_path / "fort.63", [1.0], v)
    tracemalloc.start()
    s = read_adcirc_output(p, step=0)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert s.values.shape == (1, n)
    assert peak <= 3.0 * s.values.nbytes, (peak, s.values.nbytes)

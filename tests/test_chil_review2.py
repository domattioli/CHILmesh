"""Second independent review, 2026-11-19 (specs/010-chil-container-mvp T051p to T051s).

Written before the fixes.
"""
from __future__ import annotations

import importlib
import os
import tracemalloc
import types

import pytest

import _chil_builders as B

chil_io = importlib.import_module("chilmesh.chil_io")
Err = chil_io.ChilFormatError


def _refuses(check, fn, *a, **k):
    with pytest.raises(Err) as ei:
        fn(*a, **k)
    assert ei.value.check == check, ei.value
    return ei.value


def _head(nodes=4):
    lines = ["tiny", f"1 {nodes}"]
    lines += [f"{i + 1} {float(i)} 0.0 0.0" for i in range(nodes)]
    lines += ["1 3 1 2 3"]
    return ("\n".join(lines) + "\n").encode()


# ---------------------------------------------------------------- T051p

_PAD_LINE = b"# trailing comment line padding\n"


def _padded_f14() -> bytes:
    pad = _PAD_LINE * ((65 * 1024 * 1024) // len(_PAD_LINE) + 1)
    return _head() + b"1\n2\n2 0\n1\n0\n0\n0\n" + pad


def test_huge_boundary_tail_node_id_zero_refused_validate_and_pack(tmp_path):
    f14 = _padded_f14()
    p = tmp_path / "a.chil"
    p.write_bytes(B.build_archive({"fort.14": f14}))
    _refuses("node-id-zero", chil_io.validate_chil, str(p))
    src = tmp_path / "fort.14"
    src.write_bytes(f14)
    del f14
    dst = tmp_path / "o.chil"
    _refuses("node-id-zero", chil_io.pack, dst, src)
    assert not dst.exists()


def test_no_tail_cap_symbols():
    assert not hasattr(chil_io, "_TAIL_CAP")


def test_streamed_boundary_scan_matches_in_any_chunking():
    tail = b"1\n2\n2 0\n1\n0\n"
    data = _head() + tail
    for size in (1, 3, 7, len(data)):
        s = chil_io._Fort14Scanner()
        with pytest.raises(Err) as ei:
            for i in range(0, len(data), size):
                s.feed(data[i:i + size])
            s.finish()
        assert ei.value.check == "node-id-zero"


# ---------------------------------------------------------------- T051q

_SHORT = ["tiny", "1 3", "1 0.0 0.0 0.0", "2 1.0 0.0 0.0", "3 0.0 1.0 0.0", "1 4 1 2 3"]
_SHORT_F14 = ("\n".join(_SHORT) + "\n").encode()


def test_short_element_row_refused_validate(tmp_path):
    p = tmp_path / "a.chil"
    p.write_bytes(B.build_archive({"fort.14": _SHORT_F14}))
    e = _refuses("fort14-parse", chil_io.validate_chil, str(p))
    assert "1" in str(e)


def test_short_element_row_refused_pack(tmp_path):
    src = tmp_path / "fort.14"
    src.write_bytes(_SHORT_F14)
    dst = tmp_path / "o.chil"
    _refuses("fort14-parse", chil_io.pack, dst, src)
    assert not dst.exists()


def test_short_element_row_refused_load(tmp_path):
    p = tmp_path / "a.chil"
    p.write_bytes(B.build_archive({"fort.14": _SHORT_F14}))
    _refuses("fort14-parse", chil_io.load_chil, str(p))


# ---------------------------------------------------------------- T051r

def test_pack_member_limit_checked_before_read(tmp_path):
    src = tmp_path / "fort.14"
    with open(src, "wb") as fh:
        fh.truncate(200 * 1024 * 1024)
    dst = tmp_path / "o.chil"
    tracemalloc.start()
    try:
        with pytest.raises(Err) as ei:
            chil_io.pack(dst, src, limits=chil_io.ChilLimits(max_member_bytes=1000))
        _cur, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert ei.value.check == "limit-member-bytes"
    assert peak < 16 * 1024 * 1024
    assert not dst.exists()


# ---------------------------------------------------------------- T051s

def _src(tmp_path):
    p = tmp_path / "s.chil"
    p.write_bytes(B.build_archive({"fort.14": B.good_fort14_bytes()}))
    return str(p)


_BAD = [
    ("ext", [["a"]]),
    ("ext", 5),
    ("remove_parts", 5),
    ("remove_parts", "fort13"),
    ("set_parts", [["a"]]),
    ("set_parts", 5),
    ("ext", {"demo": {"m": types.MappingProxyType({"k": 1})}}),
    ("ext", {"demo": {"big": 10 ** 5000}}),
]


@pytest.mark.parametrize("name,value", _BAD)
def test_rewrite_bad_request_types_refused(tmp_path, name, value):
    src = _src(tmp_path)
    dst = tmp_path / "o.chil"
    _refuses("request", chil_io.rewrite_chil, src, dst, **{name: value})
    assert not dst.exists()
    assert [f for f in os.listdir(tmp_path) if f not in ("s.chil",)] == []


def _nest(depth):
    obj = "x"
    for _ in range(depth):
        obj = [obj]
    return {"k": obj}


def test_rewrite_nesting_never_fails_reparse(tmp_path):
    src = _src(tmp_path)
    seen = set()
    for d in range(55, 70):
        dst = tmp_path / f"o{d}.chil"
        try:
            chil_io.rewrite_chil(src, dst, ext={"demo": _nest(d)})
            seen.add("ok")
        except Err as exc:
            assert exc.check == "request", (d, exc)
            assert not dst.exists()
            seen.add("refused")
    assert seen == {"ok", "refused"}

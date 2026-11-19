"""Tests for the executive rulings of 2026-11-19 (specs/010-chil-container-mvp T051e to T051n).

Written before the fixes. Refusals assert ``ChilFormatError`` and the exact slug.
"""
from __future__ import annotations

import importlib
import json
import os
import struct
import time
import tracemalloc
import zipfile

import pytest

import _chil_builders as B

chil_io = importlib.import_module("chilmesh.chil_io")
Err = chil_io.ChilFormatError


def _refuses(check, fn, *a, **k):
    with pytest.raises(Err) as ei:
        fn(*a, **k)
    assert ei.value.check == check, ei.value
    return ei.value


def _tiny_lines(nodes=3):
    lines = ["tiny", f"1 {nodes}"]
    lines += [f"{i + 1} {float(i)} 0.0 0.0" for i in range(nodes)]
    lines += ["1 3 1 2 3"]
    return lines


def _f14(tail):
    return ("\n".join(_tiny_lines() + list(tail)) + "\n").encode()


def _archive(tmp_path, f14, name="a.chil"):
    p = tmp_path / name
    p.write_bytes(B.build_archive({"fort.14": f14}))
    return str(p)


def _good_chil(tmp_path, name="good.chil", **kw):
    p = tmp_path / name
    p.write_bytes(B.build_archive({"fort.14": B.good_fort14_bytes()}, **kw))
    return str(p)


def _members(path):
    with zipfile.ZipFile(path) as zf:
        return {i.filename: zf.read(i) for i in zf.infolist()}


def _manifest(path):
    return json.loads(_members(path)["manifest.json"])


# ------------------------------------------------------------- M1 (T051e)

_PARTIAL_OPEN = ["1", "2", "2 0"]  # NOPE=1, total 2, segment of 2 nodes


@pytest.mark.parametrize("rows,slug", [
    (["0", "garbage"], "node-id-zero"),
    (["9", "garbage"], "node-id-range"),
    (["1", "0", "x"], "node-id-zero"),
])
def test_boundary_bad_id_before_unparseable_row_refused(tmp_path, rows, slug):
    data = _f14(_PARTIAL_OPEN + rows)
    _refuses(slug, chil_io.validate_chil, _archive(tmp_path, data))


def test_boundary_bad_id_after_unparseable_row_not_refused(tmp_path):
    data = _f14(_PARTIAL_OPEN + ["garbage", "0"])
    chil_io.validate_chil(_archive(tmp_path, data))


def test_partial_block_complete_open_missing_land_node_zero(tmp_path):
    data = _f14(["1", "2", "2 0", "1", "0"])  # land block absent
    exc = _refuses("node-id-zero", chil_io.validate_chil, _archive(tmp_path, data))
    assert "open" in str(exc)


def test_pack_boundary_bad_id_before_unparseable_row_refused(tmp_path):
    dst = tmp_path / "o.chil"
    for i, tail in enumerate((_PARTIAL_OPEN + ["0", "garbage"],
                              ["1", "2", "2 0", "1", "0"])):
        src = tmp_path / f"x{i}"
        src.write_bytes(_f14(tail))
        _refuses("node-id-zero", chil_io.pack, dst, src)
        assert not dst.exists()


def test_pack_boundary_bad_id_after_unparseable_row_not_refused(tmp_path):
    src = tmp_path / "x"
    src.write_bytes(_f14(_PARTIAL_OPEN + ["garbage", "0"]))
    chil_io.pack(tmp_path / "o.chil", src)


# ------------------------------------------------------------- M3 (T051f)

def _feed_time(make_chunks):
    t = time.perf_counter()
    s = chil_io._Fort14Scanner()
    for c in make_chunks:
        s.feed(c)
    return time.perf_counter() - t, s


def test_scanner_long_single_line_is_linear():
    size, step = 16 * 1024 * 1024, 4096
    n = size // step
    head = [b"name\n0 0\n"]  # everything after the header is the held tail
    long_line = head + [b"x" * step] * n + [b"\n"]
    normal = head + [(b"y" * (step - 1) + b"\n")] * n
    t_normal, _ = _feed_time(normal)
    t_long, _ = _feed_time(long_line)
    assert t_long <= 10 * t_normal + 1.0, (t_long, t_normal)


def test_long_first_line_32_mib_validates_in_time(tmp_path):
    size = 32 * 1024 * 1024
    body = B.good_fort14_bytes()
    long_first = b"n" * size + b"\n" + body.split(b"\n", 1)[1]
    normal = body + b"\n".join([b"filler line"] * (size // 12)) + b"\n"
    p_long = _archive(tmp_path, long_first, "long.chil")
    p_norm = _archive(tmp_path, normal, "norm.chil")
    t = time.perf_counter()
    chil_io.validate_chil(p_norm)
    t_norm = time.perf_counter() - t
    t = time.perf_counter()
    chil_io.validate_chil(p_long)
    t_long = time.perf_counter() - t
    assert t_long <= 10 * t_norm + 1.0, (t_long, t_norm)


def test_scanner_lines_split_on_cr_lf_and_crlf_across_chunks():
    data = b"tiny\r\n1 3\r\n1 0 0 0\r2 1 0 0\n3 2 0 0\r\n1 3 1 2 3\r0\r0\r0\r0\r"
    s = chil_io._Fort14Scanner()
    for i in range(0, len(data), 3):
        s.feed(data[i:i + 3])
    s.finish()
    assert s._state == 4


# ------------------------------------------------------------- L5 (T051g)

def _odd_rows(n):
    """A plain fort.14 whose single element row (id 7) has ``n`` nodes."""
    cnt = max(n, 3)
    lines = ["x", f"1 {cnt}"] + [f"{i + 1} {float(i)} {float(i % 2)} 0.0" for i in range(cnt)]
    lines += [f"7 {n} " + " ".join(str(i + 1) for i in range(n)), "0", "0", "0", "0"]
    return ("\n".join(lines) + "\n").encode()


@pytest.mark.parametrize("n", [5, 2])
def test_non_tri_quad_element_row_refused_everywhere(tmp_path, n):
    from chilmesh import CHILmesh

    f14 = _odd_rows(n)
    src = tmp_path / "t.14"
    src.write_bytes(f14)
    e = _refuses("fort14-parse", chil_io.pack, tmp_path / "o.chil", src)
    assert "7" in str(e)
    arc = _archive(tmp_path, f14)
    e = _refuses("fort14-parse", chil_io.validate_chil, arc)
    assert "7" in str(e)
    e = _refuses("fort14-parse", CHILmesh.load, arc)
    assert "7" in str(e)
    e = _refuses("fort14-parse", chil_io.summary_chil, arc, deep=True)
    assert "7" in str(e)


def test_plain_fort14_reader_still_accepts_uniform_five_node_file(tmp_path):
    from chilmesh import CHILmesh

    p = tmp_path / "t.14"
    p.write_bytes(_odd_rows(5))
    mesh = CHILmesh.read_from_fort14(str(p), compute_layers=False)
    assert mesh.connectivity_list.shape == (1, 5)


# ------------------------------------------------------------- M5 (T051h)

def _sparse(path, size):
    with open(path, "wb") as f:
        f.truncate(size)
    return path


def test_rewrite_refuses_oversize_member_by_stat_without_reading(tmp_path):
    src = _good_chil(tmp_path)
    big = _sparse(tmp_path / "big.bin", 64 * 1024 * 1024)
    dst = tmp_path / "out.chil"
    lim = chil_io.ChilLimits(max_member_bytes=1024 * 1024)
    tracemalloc.start()
    try:
        _refuses("limit-member-bytes", chil_io.rewrite_chil, src, dst,
                 set_parts={"ext/demo/big.bin": big}, ext={"demo": {}}, limits=lim)
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert peak < 16 * 1024 * 1024, peak
    assert not dst.exists()


def test_rewrite_refuses_total_by_stat_without_reading(tmp_path):
    src = _good_chil(tmp_path)
    big = _sparse(tmp_path / "big.bin", 8 * 1024 * 1024)
    lim = chil_io.ChilLimits(max_total_bytes=1024 * 1024)
    tracemalloc.start()
    try:
        _refuses("limit-total-bytes", chil_io.rewrite_chil, src, tmp_path / "o.chil",
                 set_parts={"fort13": big}, limits=lim)
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert peak < 16 * 1024 * 1024, peak


# ------------------------------------------------------------- M6 (T051i)

def test_ext_numbers_preserved_by_value_not_text(tmp_path):
    raw = (json.dumps({"format_version": 0, "crs": "unknown",
                       "parts": B.default_manifest({"fort.14": B.good_fort14_bytes()})["parts"],
                       "ext": {"demo": {}}}, sort_keys=True)).replace('"demo": {}', '"demo": {"v": 1E2}')
    src = tmp_path / "in.chil"
    src.write_bytes(B.build_archive({"fort.14": B.good_fort14_bytes()}, manifest=raw.encode()))
    out = tmp_path / "out.chil"
    chil_io.rewrite_chil(str(src), str(out))
    v = _manifest(str(out))["ext"]["demo"]["v"]
    assert v == 100.0 and isinstance(v, float)


def test_ext_member_bytes_kept_byte_for_byte(tmp_path):
    blob = b"1E2\r\n\x00\xff   \r\n"
    out = tmp_path / "out.chil"
    chil_io.rewrite_chil(_good_chil(tmp_path), str(out),
                         ext={"demo": {"n": 1}}, set_parts={"ext/demo/raw.txt": blob})
    assert _members(str(out))["ext/demo/raw.txt"] == blob


# ------------------------------------------------------------- L1 (T051j)

def test_grid_name_cr_only_fort14(tmp_path):
    data = _f14([]).replace(b"\n", b"\r")
    data = data[:-1] + b"\r0\r0\r0\r0\r"
    src = tmp_path / "x"
    src.write_bytes(data)
    rep = chil_io.pack(tmp_path / "o.chil", src)
    assert rep.summary["grid_name"] == "tiny"


def test_grid_name_crlf_fort14(tmp_path):
    src = tmp_path / "x"
    src.write_bytes(_f14(["0", "0", "0", "0"]).replace(b"\n", b"\r\n"))
    assert chil_io.pack(tmp_path / "o.chil", src).summary["grid_name"] == "tiny"


def test_grid_name_very_long_first_line_truncated_at_char_boundary(tmp_path):
    name = "é" * 3000
    lines = _tiny_lines()
    lines[0] = name
    src = tmp_path / "x"
    src.write_bytes(("\n".join(lines + ["0", "0", "0", "0"]) + "\n").encode())
    rep = chil_io.pack(tmp_path / "o.chil", src)
    gn = rep.summary["grid_name"]
    assert len(gn.encode("utf-8")) <= 256 and len(gn) == 128
    assert name.startswith(gn)
    chil_io.validate_chil(str(tmp_path / "o.chil"))


# ------------------------------------------------------------- L2 (T051k)

def _reorder_central_directory(data: bytes) -> bytes:
    eocd = data.rfind(b"PK\x05\x06")
    n, size, off = struct.unpack_from("<HII", data, eocd + 10)
    recs, pos = [], off
    for _ in range(n):
        fn, ex, cm = struct.unpack_from("<HHH", data, pos + 28)
        end = pos + 46 + fn + ex + cm
        recs.append(data[pos:end])
        pos = end
    return data[:off] + b"".join(reversed(recs)) + data[pos:]


def test_manifest_listed_first_in_central_directory_but_not_first_locally(tmp_path):
    f14 = B.good_fort14_bytes()
    raw = B.build_archive({"fort.14": f14}, order=["fort.14", "manifest.json"])
    raw = _reorder_central_directory(raw)
    p = tmp_path / "x.chil"
    p.write_bytes(raw)
    with zipfile.ZipFile(p) as zf:
        assert zf.infolist()[0].filename == "manifest.json"
        assert zf.infolist()[0].header_offset > 0
    _refuses("manifest-not-first", chil_io.validate_chil, str(p))


# ------------------------------------------------------------- L4 (T051l)

def _db_pack(tmp_path, base, name="o.chil"):
    src = tmp_path / "f14"
    src.write_bytes(B.good_fort14_bytes())
    db = tmp_path / base
    db.write_bytes(b"poly")
    return chil_io.pack(tmp_path / name, src, domain_boundary=db)


def test_pack_domain_boundary_basename_200_ok_201_refused(tmp_path):
    _db_pack(tmp_path, "x" * 200)
    dst = tmp_path / "o2.chil"
    _refuses("member-name", _db_pack, tmp_path, "y" * 201, "o2.chil")
    assert not dst.exists()


def test_pack_domain_boundary_basename_201_bytes_multibyte(tmp_path):
    _refuses("member-name", _db_pack, tmp_path, "é" * 101, "o2.chil")


def test_rewrite_domain_boundary_basename_201_refused(tmp_path):
    src = _good_chil(tmp_path)
    dst = tmp_path / "o.chil"
    spec = chil_io.MemberSpec(b"poly", name="z" * 201)
    _refuses("member-name", chil_io.rewrite_chil, src, dst, set_parts={"domain_boundary": spec})
    assert not dst.exists()
    chil_io.rewrite_chil(src, dst, set_parts={"domain_boundary": chil_io.MemberSpec(b"p", name="z" * 200)})


def test_rewrite_domain_boundary_path_basename_201_refused(tmp_path):
    src = _good_chil(tmp_path)
    long = tmp_path / ("w" * 201)
    long.write_bytes(b"poly")
    _refuses("member-name", chil_io.rewrite_chil, src, tmp_path / "o.chil",
             set_parts={"domain_boundary": long})


def test_bytes_paths_refused_with_request(tmp_path):
    f14 = tmp_path / "f14"
    f14.write_bytes(B.good_fort14_bytes())
    good = _good_chil(tmp_path)
    _refuses("request", chil_io.pack, tmp_path / "o.chil", os.fsencode(f14))
    _refuses("request", chil_io.pack, os.fsencode(tmp_path / "o.chil"), f14)
    _refuses("request", chil_io.pack, tmp_path / "o.chil", f14, domain_boundary=os.fsencode(f14))
    _refuses("request", chil_io.rewrite_chil, os.fsencode(good), tmp_path / "o.chil")
    _refuses("request", chil_io.rewrite_chil, good, os.fsencode(tmp_path / "o.chil"))
    _refuses("request", chil_io.validate_chil, os.fsencode(good))
    assert not (tmp_path / "o.chil").exists()


# ------------------------------------------------------------- L6 (T051m)

BAD_CRS = ["", "a" * 257, "bad\x00", "tab\there", 5]


@pytest.mark.parametrize("crs", BAD_CRS)
def test_bad_crs_on_save_is_request(tmp_path, crs):
    from chilmesh import examples

    dst = tmp_path / "o.chil"
    _refuses("request", chil_io.save_chil, examples.annulus(), dst, crs=crs)
    assert not dst.exists()


@pytest.mark.parametrize("crs", BAD_CRS)
def test_bad_crs_on_rewrite_is_request(tmp_path, crs):
    dst = tmp_path / "o.chil"
    _refuses("request", chil_io.rewrite_chil, _good_chil(tmp_path), dst, crs=crs)
    assert not dst.exists()


@pytest.mark.parametrize("crs", BAD_CRS + [None])
def test_bad_crs_on_pack_is_request(tmp_path, crs):
    src = tmp_path / "f14"
    src.write_bytes(B.good_fort14_bytes())
    dst = tmp_path / "o.chil"
    _refuses("request", chil_io.pack, dst, src, crs=crs)
    assert not dst.exists()


@pytest.mark.parametrize("value", [5, ["a"], None])
def test_non_string_domain_boundary_key_is_domain_boundary(tmp_path, value):
    f14 = B.good_fort14_bytes()
    m = B.default_manifest({"fort.14": f14})
    m["domain_boundary"] = value
    p = tmp_path / "x.chil"
    p.write_bytes(B.build_archive({"fort.14": f14}, manifest=m))
    _refuses("domain-boundary", chil_io.validate_chil, str(p))


# ------------------------------------------------------- contract additions

def test_remove_missing_part_or_ext_member_is_request(tmp_path):
    src = _good_chil(tmp_path)
    for key in ("fort13", "domain_boundary", "ext/demo/none.bin"):
        _refuses("request", chil_io.rewrite_chil, src, tmp_path / "o.chil", remove_parts=[key])
    assert not (tmp_path / "o.chil").exists()


def test_ext_key_dropped_when_no_namespace_remains(tmp_path):
    a, b = tmp_path / "a.chil", tmp_path / "b.chil"
    chil_io.rewrite_chil(_good_chil(tmp_path), a, ext={"demo": {"k": 1}})
    assert "ext" in _manifest(str(a))
    chil_io.rewrite_chil(a, b, ext={"demo": None})
    assert "ext" not in _manifest(str(b))

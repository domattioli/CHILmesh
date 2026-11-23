"""NOAA auxiliary fetch records URL, SHA-256 and fetch date (plan D1)."""
from __future__ import annotations

import hashlib
import io

from benchmarks.competitors import noaa


class _Resp(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def test_fetch_aux_records_provenance(tmp_path):
    payload = b"netcdf-bytes" * 100
    res = noaa.fetch_aux("https://noaa.example/points.cwl.nc", tmp_path / "points.cwl.nc",
                         opener=lambda req, timeout=None: _Resp(payload),
                         today="2026-11-08", chunk_bytes=64)
    assert res["url"] == "https://noaa.example/points.cwl.nc"
    assert res["sha256"] == hashlib.sha256(payload).hexdigest()
    assert res["fetch_date"] == "2026-11-08"
    assert (tmp_path / "points.cwl.nc").read_bytes() == payload


def test_fetch_aux_failure_is_recorded_not_raised(tmp_path):
    def opener(req, timeout=None):
        raise OSError("offline")

    res = noaa.fetch_aux("https://noaa.example/x", tmp_path / "x", opener=opener, today="2026-11-08")
    assert res["status"] == "failed"
    assert "offline" in res["reason"]


def test_fetch_aux_rejects_non_https(tmp_path):
    res = noaa.fetch_aux("file:///etc/passwd", tmp_path / "x", opener=lambda *a, **k: None, today="2026-11-08")
    assert res["status"] == "failed"

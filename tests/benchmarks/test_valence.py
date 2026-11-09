"""Valence wrapper: streaming download, not-served handling, rights filter, no fallback."""
from __future__ import annotations

import hashlib
import io
import json
import types

import pytest

from benchmarks.competitors import valence as V


class _FakeResp(io.BytesIO):
    def __init__(self, data):
        super().__init__(data)
        self.max_read = 0

    def read(self, n=-1):
        self.max_read = max(self.max_read, n if n and n > 0 else 10**12)
        return super().read(n)

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class _NotFound(Exception):
    pass


class _Unavailable(Exception):
    pass


class _Auth(Exception):
    pass


def _fake_registry(digest="", *, hash_exc=None):
    def fetch_hash(name):
        if hash_exc:
            raise hash_exc
        return digest

    return types.SimpleNamespace(
        RegistryNotFoundError=_NotFound,
        RegistryUnavailableError=_Unavailable,
        RegistryAuthError=_Auth,
        RegistryError=Exception,
        _config=lambda: ("https://valence.example", "tok"),
        _check_name=lambda n: n,
        fetch_hash=fetch_hash,
        list_meshes=lambda: ["Shinnecock"],
    )


def test_streaming_download_is_chunked_and_verified(tmp_path):
    payload = b"x" * (3 * 1024 + 7)
    reg = _fake_registry(hashlib.sha256(payload).hexdigest())
    resp = _FakeResp(payload)
    calls = []

    def opener(req, timeout=None):
        calls.append(req.full_url)
        return resp

    res = V.fetch_mesh_file("Shinnecock", tmp_path, registry=reg, opener=opener, chunk_bytes=1024)
    assert res["status"] == "ok"
    assert res["sha256"] == hashlib.sha256(payload).hexdigest()
    assert (tmp_path / "Shinnecock.14").read_bytes() == payload
    assert resp.max_read <= 1024  # never reads the whole payload at once
    assert calls == ["https://valence.example/meshes/Shinnecock"]


def test_hash_mismatch_is_failed_and_file_removed(tmp_path):
    reg = _fake_registry("0" * 64)
    res = V.fetch_mesh_file("Shinnecock", tmp_path, registry=reg,
                            opener=lambda req, timeout=None: _FakeResp(b"abc"))
    assert res["status"] == "failed"
    assert "hash" in res["reason"]
    assert not (tmp_path / "Shinnecock.14").exists()


@pytest.mark.parametrize("exc", [_NotFound("nope"), _Unavailable("down")])
def test_not_served_or_unavailable_is_not_run_without_fallback(tmp_path, exc):
    reg = _fake_registry(hash_exc=exc)
    opened = []

    def opener(req, timeout=None):  # must never be reached
        opened.append(req)
        raise AssertionError("fallback download attempted")

    res = V.fetch_mesh_file("STOFS", tmp_path, registry=reg, opener=opener)
    assert res["status"] == "not-run"
    assert res["reason"] == "mesh not yet served by Valence"
    assert opened == []
    assert list(tmp_path.iterdir()) == []


def test_missing_token_is_not_run_with_distinct_reason(tmp_path):
    reg = _fake_registry(hash_exc=_Auth("no token"))
    res = V.fetch_mesh_file("Shinnecock", tmp_path, registry=reg, opener=None)
    assert res["status"] == "not-run"
    assert "VALENCE_TOKEN" in res["reason"]


def test_module_has_no_noaa_or_local_clone_fallback():
    src = open(V.__file__, encoding="utf-8").read().lower()
    assert "noaa" not in src
    assert "valence_data_dir" not in src


def test_rights_manual_file_excludes_unpermitted_meshes(tmp_path):
    rights = tmp_path / "rights.json"
    rights.write_text(json.dumps({
        "Shinnecock": {"publish_results": True, "note": "ok"},
        "WNAT_Hagen": {"publish_results": False, "note": "no"},
    }))
    allowed, excluded = V.apply_rights(["Shinnecock", "WNAT_Hagen", "ENPAC2003"], V.load_rights(rights))
    assert allowed == ["Shinnecock"]
    assert {e["mesh"] for e in excluded} == {"WNAT_Hagen", "ENPAC2003"}
    assert all(e["reason"] for e in excluded)


def test_rights_file_missing_excludes_everything(tmp_path):
    allowed, excluded = V.apply_rights(["Shinnecock"], V.load_rights(tmp_path / "absent.json"))
    assert allowed == []
    assert excluded[0]["mesh"] == "Shinnecock"

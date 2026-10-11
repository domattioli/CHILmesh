"""Tests for chilmesh.registry. No network: urlopen is mocked."""
import hashlib
import io
import json
import urllib.error
from pathlib import Path
from unittest import mock

import pytest

from chilmesh import registry as reg

TOK = "s3cr3t-token-xyz"
FIXTURE = Path(reg.__file__).parent / "data" / "quad_2x2.fort.14"


class Resp:
    def __init__(self, body):
        self._b = body

    def read(self):
        return self._b

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def http_err(code, headers=None):
    return urllib.error.HTTPError("http://x", code, "m", headers or {}, io.BytesIO(b""))


@pytest.fixture(autouse=True)
def env(monkeypatch):
    monkeypatch.setenv("VALENCE_TOKEN", TOK)
    monkeypatch.setenv("VALENCE_API_URL", "https://api.example.test/")


def patch(side):
    if isinstance(side, list):
        return mock.patch("urllib.request.urlopen", side_effect=side)
    return mock.patch("urllib.request.urlopen", side_effect=side)


def J(o):
    return Resp(json.dumps(o).encode())


def test_list_ok():
    with patch([J(["a.14", "b.14"])]) as m:
        assert reg.list_meshes() == ["a.14", "b.14"]
    req = m.call_args[0][0]
    assert req.full_url == "https://api.example.test/meshes"
    assert req.get_header("Authorization") == f"Bearer {TOK}"
    assert TOK not in req.full_url


@pytest.mark.parametrize("body", [{"a": 1}, [1, 2], ["a", 2]])
def test_list_malformed(body):
    with patch([J(body)]):
        with pytest.raises(reg.RegistryIntegrityError):
            reg.list_meshes()


def test_list_non_json():
    with patch([Resp(b"<html>")]):
        with pytest.raises(reg.RegistryIntegrityError):
            reg.list_meshes()


def test_hash_ok():
    h = "a" * 64
    with patch([J({"sha256": h})]):
        assert reg.fetch_hash("m.14") == h


@pytest.mark.parametrize(
    "body",
    [{"sha256": "abc"}, {"sha256": "g" * 64}, {"sha256": "a" * 64 + "\n"}, {"sha256": "A" * 64}, [1], {"sha256": "a" * 64, "size": 3}, {}],
)
def test_hash_malformed(body):
    with patch([J(body)]):
        with pytest.raises(reg.RegistryIntegrityError):
            reg.fetch_hash("m.14")


def test_bytes_ok():
    data = b"hello"
    with patch([J({"sha256": hashlib.sha256(data).hexdigest()}), Resp(data)]):
        assert reg.fetch_mesh_bytes("m.14") == data


def test_integrity_mismatch():
    with patch([J({"sha256": "0" * 64}), Resp(b"hello")]):
        with pytest.raises(reg.RegistryIntegrityError):
            reg.fetch_mesh_bytes("m.14")


@pytest.mark.parametrize(
    "code,exc",
    [
        (401, reg.RegistryAuthError),
        (404, reg.RegistryNotFoundError),
        (429, reg.RegistryUnavailableError),
        (503, reg.RegistryUnavailableError),
    ],
)
def test_status_mapping(code, exc):
    with patch(http_err(code)):
        with pytest.raises(exc) as ei:
            reg.list_meshes()
    assert TOK not in str(ei.value)


def test_405_is_hard_error():
    # Contract: 401/404/405 are hard errors; 405 must not look transient.
    with patch(http_err(405)):
        with pytest.raises(reg.RegistryError) as ei:
            reg.list_meshes()
    assert type(ei.value) is reg.RegistryError


def test_429_retry_after():
    with patch(http_err(429, {"Retry-After": "42"})):
        with pytest.raises(reg.RegistryUnavailableError, match="42"):
            reg.list_meshes()


def test_not_found_is_value_error():
    assert issubclass(reg.RegistryNotFoundError, ValueError)
    assert issubclass(reg.RegistryNotFoundError, reg.RegistryError)


@pytest.mark.parametrize(
    "side", [urllib.error.URLError("boom"), TimeoutError(), OSError("x")]
)
def test_network_errors(side):
    with patch(side):
        with pytest.raises(reg.RegistryUnavailableError):
            reg.list_meshes()


def test_missing_token(monkeypatch):
    monkeypatch.delenv("VALENCE_TOKEN")
    with patch([]) as m:
        with pytest.raises(reg.RegistryAuthError, match="VALENCE_TOKEN"):
            reg.list_meshes()
    assert m.call_count == 0


def test_default_url_used_when_unset(monkeypatch):
    monkeypatch.delenv("VALENCE_API_URL")
    with patch([J([])]) as m:
        reg.list_meshes()
    assert m.call_args[0][0].full_url.startswith("https://valence-api.iqh766c2.workers.dev/")


def test_url_override_wins():
    with patch([J([])]) as m:
        reg.list_meshes()
    assert m.call_args[0][0].full_url.startswith("https://api.example.test/")


def test_user_agent_not_urllib_default():
    with patch([J([])]) as m:
        reg.list_meshes()
    ua = m.call_args[0][0].get_header("User-agent")
    assert ua == "chilmesh-registry/1"
    assert not ua.startswith("Python-urllib")


@pytest.mark.parametrize("bad", ["../x", "a/b", "", "..", ".", "a b", "a%2Fb"])
def test_bad_names(bad):
    with patch([]) as m:
        for fn in (reg.fetch_hash, reg.fetch_mesh_bytes, reg.fetch_mesh):
            with pytest.raises(reg.RegistryNotFoundError):
                fn(bad)
        with pytest.raises(reg.RegistryNotFoundError):
            reg.download_mesh(bad, "/nonexistent/x")
    assert m.call_count == 0


def test_headers_on_all_calls():
    data = b"zz"
    with patch([J({"sha256": hashlib.sha256(data).hexdigest()}), Resp(data)]) as m:
        reg.fetch_mesh_bytes("m.14")
    for c in m.call_args_list:
        r = c[0][0]
        assert r.get_header("Authorization") == f"Bearer {TOK}"
        assert TOK not in r.full_url
        assert c[1]["timeout"] == 30.0


def test_download_mesh(tmp_path):
    data = b"abc"
    with patch([J({"sha256": hashlib.sha256(data).hexdigest()}), Resp(data)]):
        p = reg.download_mesh("m.14", tmp_path / "o.14")
    assert p.read_bytes() == data


def test_fetch_mesh_real(tmp_path):
    from chilmesh import CHILmesh

    data = FIXTURE.read_bytes()
    with patch([J({"sha256": hashlib.sha256(data).hexdigest()}), Resp(data)]):
        mesh = reg.fetch_mesh("quad.14", compute_layers=False)
    assert isinstance(mesh, CHILmesh)
    assert mesh.n_elems > 0

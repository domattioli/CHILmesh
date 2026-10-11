"""Drive every GET transcript case through the CHILmesh registry client."""
import hashlib
import io
import json
import urllib.error
from email.message import Message
from pathlib import Path

import pytest

from chilmesh import registry

CONTRACT = Path(__file__).parent / "data" / "registry-contract.v1.json"
CONTRACT_SHA256 = "702ecccf7d5ac5a2e043581479babf5652918cf028fab49428d14ce590f37f58"
CASES = [entry for entry in json.loads(CONTRACT.read_text())["transcript"] if entry["request"]["method"] == "GET"]


class _Response:
    def __init__(self, data): self.data = data
    def read(self): return self.data
    def __enter__(self): return self
    def __exit__(self, *args): return False


def _body(response):
    return response.get("body_text", "").encode() if "body_text" in response else json.dumps(response["body"]).encode()


@pytest.mark.parametrize("entry", CASES)
def test_get_transcript(monkeypatch, entry):
    assert hashlib.sha256(CONTRACT.read_bytes()).hexdigest() == CONTRACT_SHA256
    contract = json.loads(CONTRACT.read_text())
    assert contract["errors"]["405"] == {"error": "method_not_allowed"}
    monkeypatch.setenv("VALENCE_TOKEN", "token")
    response, calls = entry["response"], []
    hash_entry = next(item for item in CASES if item["request"]["path"].endswith("/hash") and item["response"]["status"] == 200)

    def urlopen(request, timeout):
        calls.append((request, timeout))
        current = hash_entry["response"] if request.full_url.endswith("/WNAT.14/hash") and entry["request"]["path"] == "/meshes/WNAT.14" else response
        status = current["status"]
        if status == 200:
            return _Response(_body(current))
        headers = Message()
        if status == 429:
            headers["Retry-After"] = "60"
        raise urllib.error.HTTPError(request.full_url, status, "contract", headers, io.BytesIO(_body(current)))

    monkeypatch.setattr(registry.urllib.request, "urlopen", urlopen)
    path, status = entry["request"]["path"], response["status"]
    adapter = registry.list_meshes if path == "/meshes" else (lambda: registry.fetch_hash("WNAT.14") if path.endswith("/hash") else registry.fetch_mesh_bytes("WNAT.14"))
    expected_error = {401: registry.RegistryAuthError, 404: registry.RegistryNotFoundError, 429: registry.RegistryUnavailableError, 503: registry.RegistryUnavailableError}.get(status)
    if expected_error:
        with pytest.raises(expected_error): adapter()
    else:
        got = adapter()
        expected = response.get("body_text", response.get("body"))
        if path == "/meshes/WNAT.14":
            assert got == expected.encode()
            assert hashlib.sha256(got).hexdigest() == hash_entry["response"]["body"]["sha256"]
        else:
            assert got == (expected if path == "/meshes" else expected["sha256"])
    assert calls
    assert all(request.get_header("Authorization") == "Bearer token" and not str(request.get_header("User-agent") or "Python-urllib").startswith("Python-urllib") and timeout == 30.0 for request, timeout in calls)


def test_malformed_hash_is_integrity_error(monkeypatch):
    monkeypatch.setenv("VALENCE_TOKEN", "token")
    monkeypatch.setattr(registry.urllib.request, "urlopen", lambda request, timeout: _Response(b'{"sha256":"x","extra":1}'))
    with pytest.raises(registry.RegistryIntegrityError):
        registry.fetch_hash("WNAT.14")

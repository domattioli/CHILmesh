import io
import pytest
from benchmarks.competitors import protocol as P


def test_state_sequence_and_validation():
    state = P.ProtocolState("r")
    for i, event in enumerate(("READY", "GO", "START", "END", "ACK")):
        state.accept(dict(protocol_version=1, repeat_id="r", event=event, monotonic_s=float(i), clock="perf_counter", clock_metadata=P.clock_metadata()))
    assert state.complete
    with pytest.raises(P.ProtocolError): state.accept(dict(protocol_version=1, repeat_id="r", event="END", monotonic_s=6., clock="perf_counter", clock_metadata=P.clock_metadata()))


def test_bad_identity_and_order():
    state = P.ProtocolState("r")
    with pytest.raises(P.ProtocolError): state.accept(dict(protocol_version=1, repeat_id="other", event="READY", monotonic_s=0., clock="perf_counter", clock_metadata=P.clock_metadata()))
    with pytest.raises(P.ProtocolError): state.accept(dict(protocol_version=1, repeat_id="r", event="START", monotonic_s=0., clock="perf_counter", clock_metadata=P.clock_metadata()))
    assert not state.complete


def test_child_ack_gates_verification_and_result():
    import json
    incoming = io.StringIO(json.dumps(message("GO", 1)) + "\n" + json.dumps(message("ACK", 3)) + "\n")
    outgoing = io.StringIO(); results = io.StringIO()
    child = P.ChildProtocol("r", incoming, outgoing, results, clock=iter([0., 2., 2.5]).__next__)
    child.ready()
    with child.region():
        value = [1, 2, 3]
    child.result({"value": value, "payload": "x" * 100000})
    assert len(results.getvalue()) > 100000
    assert [__import__('json').loads(x)["event"] for x in outgoing.getvalue().splitlines()] == ["READY", "START", "END"]
    assert "new_modules" in __import__('json').loads(results.getvalue())


def message(event, stamp=0, **changes):
    value = dict(protocol_version=1, repeat_id="r", event=event, monotonic_s=stamp, clock="perf_counter", clock_metadata=P.clock_metadata())
    value.update(changes)
    return value


@pytest.mark.parametrize("change", [dict(repeat_id="wrong"), dict(clock="other"),
                                    dict(protocol_version=2), dict(monotonic_s=float("nan")),
                                    dict(monotonic_s=True)])
def test_marker_identity_and_timestamp_rejected(change):
    with pytest.raises(P.ProtocolError):
        P.ProtocolState("r").accept(message("READY", **change))


def test_missing_markers_are_rejected_at_finish():
    state = P.ProtocolState("r")
    state.accept(message("READY"))
    with pytest.raises(P.ProtocolError):
        state.finish()


def test_duplicate_and_reversed_timestamps_rejected():
    state = P.ProtocolState("r")
    state.accept(message("READY", 2))
    with pytest.raises(P.ProtocolError):
        state.accept(message("READY", 3))
    with pytest.raises(P.ProtocolError):
        state.accept(message("GO", 1))


def test_malformed_control_and_missing_ack_block_verification():
    incoming = io.StringIO('{"event": "GO"}\n')
    child = P.ChildProtocol("r", incoming, io.StringIO(), io.StringIO(), clock=lambda: 0)
    with pytest.raises(P.ProtocolError):
        child.ready()
    with pytest.raises(P.ProtocolError):
        child.result({})


def test_region_tracks_new_modules_and_materializes_before_end(monkeypatch):
    import json
    import sys
    incoming = io.StringIO(json.dumps(message("GO", 1)) + "\n" + json.dumps(message("ACK", 4)) + "\n")
    output, results = io.StringIO(), io.StringIO()
    child = P.ChildProtocol("r", incoming, output, results, clock=iter([0, 2, 3]).__next__)
    child.ready()
    with child.region():
        value = list(range(100))
        monkeypatch.setitem(sys.modules, "fair_test_deferred_module", object())
        with pytest.raises(P.ProtocolError):
            child.result({"value": value})
    child.result({"value": value})
    result = json.loads(results.getvalue())
    assert result["value"] == list(range(100))
    assert result["new_modules"] == ["fair_test_deferred_module"]


def test_finish_accepts_complete_sequence():
    state = P.ProtocolState("r")
    for stamp, event in enumerate(state.sequence):
        state.accept(message(event, stamp))
    assert len(state.finish()) == 5


def test_invalid_json_control_is_protocol_error():
    child = P.ChildProtocol("r", io.StringIO("not json\n"), io.StringIO(), io.StringIO(), clock=lambda: 0)
    with pytest.raises(P.ProtocolError):
        child.ready()


def test_missing_ack_prevents_verification():
    import json
    incoming = io.StringIO(json.dumps(message("GO", 1)) + "\n")
    child = P.ChildProtocol("r", incoming, io.StringIO(), io.StringIO(), clock=iter([0, 2, 3]).__next__)
    child.ready()
    verified = []
    with pytest.raises(P.ProtocolError):
        with child.region():
            pass
        verified.append(True)
    assert not verified
    with pytest.raises(P.ProtocolError):
        child.result({})


@pytest.mark.parametrize("value", [True, 1.0, "1"])
def test_protocol_version_is_integer(value):
    with pytest.raises(P.ProtocolError): P.ProtocolState("r").accept(message("READY", protocol_version=value))


def test_state_retains_copy_and_parent_receipt_timestamp():
    marker = message("READY", 1)
    state = P.ProtocolState("r")
    state.accept(marker, received_s=2)
    marker["event"] = "END"
    assert state.events[0]["event"] == "READY"
    assert state.receipts == [2]


def test_clock_metadata_identifies_monotonic_shared_clock():
    metadata = P.clock_metadata()
    assert metadata["monotonic"] is True
    assert metadata["implementation"]
    assert metadata["resolution_s"] > 0


def test_region_run_materializes_then_ack_then_verifies():
    import json
    incoming = io.StringIO(json.dumps(message("GO", 1)) + "\n" + json.dumps(message("ACK", 4)) + "\n")
    outgoing, results = io.StringIO(), io.StringIO()
    child = P.ChildProtocol("r", incoming, outgoing, results, clock=iter([0, 2, 3]).__next__)
    order = []
    child.ready()
    def materialize(value):
        assert '"event": "END"' not in outgoing.getvalue()
        order.append("load")
        return list(value)
    def verify(value):
        assert child.state.complete
        order.append("verify")
        return {"value": value}
    child.run(lambda: iter([1, 2]), materialize=materialize, verify=verify)
    assert order == ["load", "verify"]
    assert json.loads(results.getvalue())["value"] == [1, 2]


def test_result_can_only_be_written_once():
    import json
    incoming = io.StringIO(json.dumps(message("GO", 1)) + "\n" + json.dumps(message("ACK", 4)) + "\n")
    child = P.ChildProtocol("r", incoming, io.StringIO(), io.StringIO(), clock=iter([0, 2, 3]).__next__)
    child.ready()
    with child.region(): pass
    child.result({})
    with pytest.raises(P.ProtocolError): child.result({})


def test_protocol_imports_only_stdlib():
    import ast
    from pathlib import Path
    import sys
    tree = ast.parse(Path(P.__file__).read_text())
    imports = [n.module.split('.')[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)]
    imports += [a.name.split('.')[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names]
    assert set(imports) <= sys.stdlib_module_names


@pytest.mark.parametrize("metadata", [None, {}, {"monotonic": False, "implementation": "clock", "resolution_s": 1e-9},
                                      {"monotonic": True, "implementation": "clock", "resolution_s": 0}])
def test_markers_require_complete_clock_metadata(metadata):
    marker = message("READY"); marker["clock_metadata"] = metadata
    with pytest.raises(P.ProtocolError): P.ProtocolState("r").accept(marker)


def test_clock_build_cannot_change_between_markers():
    state = P.ProtocolState("r"); state.accept(message("READY", 0))
    marker = message("GO", 1); marker["clock_metadata"] = dict(P.clock_metadata(), implementation="another clock")
    with pytest.raises(P.ProtocolError): state.accept(marker)


def test_region_exception_never_verifies_or_writes_result():
    import json
    incoming = io.StringIO(json.dumps(message("GO", 1)) + "\n")
    results = io.StringIO()
    child = P.ChildProtocol("r", incoming, io.StringIO(), results, clock=iter([0, 2]).__next__)
    child.ready()
    with pytest.raises(RuntimeError):
        child.run(lambda: (_ for _ in ()).throw(RuntimeError("operation failed")), materialize=list, verify=lambda value: {"value": value})
    assert not child.state.complete and results.getvalue() == ""


def test_child_can_import_from_repo_root_without_chilmesh_or_numpy(tmp_path):
    import subprocess
    import sys
    from pathlib import Path
    root = Path(P.__file__).resolve().parents[2]
    script = "import sys; sys.path.insert(0, sys.argv[1]); from benchmarks.competitors.protocol import ChildProtocol; assert 'numpy' not in sys.modules; assert 'chilmesh' not in sys.modules"
    result = subprocess.run([sys.executable, "-I", "-c", script, str(root)], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr

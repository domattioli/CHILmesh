import math
import pytest
from benchmarks.competitors import records as R


def sample():
    from test_records import fair_data
    return fair_data()


def test_v2_roundtrip_immutable_and_hash(tmp_path):
    rec = R.make_fair_record(**sample())
    path = tmp_path / "r.json"
    digest = R.write_fair_record(rec, path)
    assert R.read_record(path, expected_sha256=digest) == rec
    with pytest.raises(FileExistsError): R.write_fair_record(rec, path)
    with pytest.raises(R.RecordError): R.read_record(path, expected_sha256="0" * 64)


def test_v2_rejects_missing_and_nonfinite():
    data = sample(); del data["series_id"]
    with pytest.raises(R.RecordError): R.make_fair_record(**data)
    data = sample(); data["measurement"] = {"value": math.inf}
    with pytest.raises(R.RecordError): R.make_fair_record(**data)
    data = sample(); data["measurement"] = {"operation_elapsed_s": None}
    with pytest.raises(R.RecordError): R.make_fair_record(**data)


def test_v2_coverage_and_legacy_exclusion():
    rec = R.make_fair_record(**sample())
    R.check_fair_coverage([rec], [("c", "cell", "s", "parse", "python")])
    with pytest.raises(R.MissingCellError): R.check_fair_coverage([rec], [("c", "cell", "s", "parse", "cpp")])
    assert not R.fair_eligible({"schema_version": 1})

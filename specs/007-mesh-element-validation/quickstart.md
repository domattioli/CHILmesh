# Quickstart: Mesh Element Validity

## Run the suite against a built-in fixture

```python
import chilmesh
from tests._validity import validate_mesh_elements

mesh = chilmesh.examples.annulus()
report = validate_mesh_elements(mesh)

print(f"ok={report.ok}  n_elems={report.n_elems_checked}  runtime={report.runtime_s:.3f}s")
print(f"violations: {len(report.violations)}")
print(f"notes: {len(report.notes)}")

for v in report.violations[:5]:
    print(f"  [{v.category}] elems={v.element_ids} edges={v.edge_ids}: {v.detail}")
```

Expected output on a healthy `annulus`:

```text
ok=True  n_elems=104  runtime=0.012s
violations: 0
notes: 0
```

## Run the pytest suite

```bash
# All four built-in fixtures + the synthetic negative fixtures:
pytest -v tests/test_mesh_element_validity.py

# Skip slow block_o:
pytest -v -k 'not block_o' tests/test_mesh_element_validity.py

# Run just the negative-fixture sanity checks:
pytest -v -k 'synthetic' tests/test_mesh_element_validity.py
```

## Override the tolerance

```python
# For lat/lon meshes with very small bbox (~1 degree), the default 1e-12 * bbox_diag
# may be too tight. Override:
report = validate_mesh_elements(mesh, tol=1e-9)
```

## Inspect violations programmatically

```python
report = validate_mesh_elements(mesh)
bowties = [v for v in report.violations if v.category == "SELF_INTERSECTING_QUAD"]
overlaps = [v for v in report.violations if v.category == "INTERIOR_OVERLAP"]
```

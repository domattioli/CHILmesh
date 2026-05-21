# tests/_validity — internal helpers (spec 007)

Test-suite-only mesh element validator. **NOT public chilmesh API.**

See `specs/007-mesh-element-validation/spec.md` for requirements.
Promotion to `chilmesh.validate` is tracked at issue #142.

## Layout

- `types.py` — `Violation`, `InformationalNote`, `MeshValidityReport`
- `predicates.py` — orientation, segment-crossing, winding, classify_element
- `broadphase.py` — `UniformGridIndex` for >5k-element meshes
- `validator.py` — `validate_mesh_elements` entry point
- `fixtures.py` — synthetic negative fixtures (bowtie, interior tri, pentagon, overlap, edge-cross)

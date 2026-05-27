# Draft: Gmsh `.msh` I/O (ported-from-ADMESH agent output)

**Status:** preservation only — not wired into CHILmesh. Parked here so the
work isn't lost while the repo-home decision (ADMESH vs CHILmesh) is open.

## What this is

A complete first-pass implementation of Gmsh ASCII `.msh` read/write,
generated against **ADMESH spec 008** (`specs/008-gmsh-io-integration/` in
`domattioli/ADMESH`). It is **admesh-shaped**: `gmsh.py` imports
`admesh.api` (`Mesh`, `BoundarySegment`) and `admesh.boundary_types`
(`BoundaryType`), so it does not run inside CHILmesh as-is.

Reported by the generating agent: **82 passed, 2 skipped** (the 2 skips are
`gmsh -check` tests gated by `shutil.which("gmsh")`); no regression in the
existing 370-test ADMESH suite.

## Contents

| File | Role |
|---|---|
| `gmsh.py` | `read_msh`, `write_msh`, `GmshParseError`, `PHYSICAL_GROUP_MAP` — v2.2 + v4.1 parsers/writers, hand-rolled (no `gmsh` PyPI dep) |
| `test_gmsh_read.py` / `test_gmsh_write.py` / `test_gmsh_roundtrip.py` | test suite |
| `fixtures/*.msh` | hand-crafted Gmsh fixtures |
| `admesh-integration.patch` | the `admesh/api.py` (`Mesh.to_msh`) + `admesh/__init__.py` (export) edits, as a diff |

## If this lands in ADMESH

```
# from the ADMESH repo root:
cp <this-dir>/gmsh.py admesh/gmsh.py
cp <this-dir>/test_gmsh_*.py tests/
mkdir -p tests/fixtures/gmsh && cp <this-dir>/fixtures/*.msh tests/fixtures/gmsh/
git apply <this-dir>/admesh-integration.patch
```

## If this is re-homed in CHILmesh

`gmsh.py` must be re-targeted onto CHILmesh's mesh class (`CHILmesh`) and its
existing fort.14 I/O instead of `admesh.api.Mesh` — i.e. a port, not a copy.
Use this as the reference for grammar handling, the physical-group ↔ boundary
mapping, and the test matrix.

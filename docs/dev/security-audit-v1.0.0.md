<!--
INTERNAL ENGINEERING ARTIFACT — NOT the project's security policy.
This is a point-in-time audit report, not a vulnerability-disclosure policy.
The disclosure policy lives at repo-root `SECURITY.md` (GitHub Security tab).
Findings A1 and A5 were "added but unverified" at audit time; do not treat
their mitigations as confirmed until the C++ backend is rebuilt and tested in
CI. NOTE: this repo is public, so this path (and git history) remain publicly
readable — if A1/A5 stay unverified, prefer a private GHSA advisory over any
tracked file. Moved out of `docs/SECURITY.md` 2026-06-25 so it stops populating
the GitHub Security tab as if it were the disclosure policy.
-->

# CHILmesh v1.0.0 Security Audit

**Date:** 2026-05-22  
**Auditor:** Adversarial code review (gsd-security-auditor agent)  
**Status:** 5 threats identified; 3 mitigations complete, 2 pending C++ rebuild verification

---

## Threat Summary

| ID | Category | Severity | Component | Disposition |
|----|----------|----------|-----------|-------------|
| A1 | Tampering / DoS | **CRITICAL** | C++ `HalfEdgeMesh::build()` | Vertex bounds check (code added, unverified) |
| A2 | Denial of Service | **HIGH** | Fort14 parser | File-size + count cap (✅ implemented) |
| A3 | Denial of Service | **HIGH** | 2dm parser | File-size + record limit (✅ implemented) |
| A5 | Memory Leak | MEDIUM | pybind11 `full_init` | Try/catch cleanup (code added, unverified) |
| A6 | Tampering / DoS | MEDIUM | cpp_backend validation | Float rejection + index validation (✅ implemented) |

---

## Critical Finding: A1 Out-of-Bounds Vertex Index

Any untrusted fort14 file with a vertex index ≥ n_verts will crash the C++ backend:

```cpp
// OLD CODE: no bounds check
he.vert = row[i];  // row[i] could be 5, n_verts could be 3
...
vert_he[v] = h;  // v=5 → out-of-bounds vector subscript → heap corruption
```

**Mitigation code added to halfedge.cpp (line 58+):**
```cpp
if (v != -1 && (v < 0 || v >= n_verts)) {
    throw std::out_of_range("vertex index ... out of range");
}
```

**Status:** Code present in repo but C++ build couldn't be verified in this environment.

---

## Implemented Mitigations

### A2: Fort14 File-Size DoS (HIGH)
✅ Prevents 2B-element crafted headers from causing 48GB allocation

```python
MAX_FILE_SIZE = 1_000_000_000  # 1 GB
if file_size > MAX_FILE_SIZE:
    raise ValueError("FORT.14 file too large")
```

### A3: 2dm Record-Count DoS (HIGH)
✅ Prevents unlimited record parsing from exhausting memory

```python
MAX_RECORDS = 100_000_000
for line in f:
    record_count += 1
    if record_count > MAX_RECORDS:
        raise ValueError("too many records")
```

### A6: Float Array Truncation (MEDIUM)
✅ Rejects float connectivity (silent truncation is unsafe)

```python
if np.issubdtype(conn.dtype, np.floating):
    raise TypeError("connectivity must be integer dtype")

if np.any(conn < -1):
    raise ValueError("negative indices other than -1")
```

---

## Unverified Mitigations (C++ Build Pending)

### A1: Vertex Index Bounds Check
- Code: halfedge.cpp lines 58-64
- Status: Added but couldn't verify compilation/exception propagation
- Action: Rebuild in CI with `pip install -e .` on Linux/Mac/Windows

### A5: Memory Leak Cleanup
- Code: bindings.cpp try/catch around build_adjacency() + skeletonize()
- Status: Added but couldn't verify exception handling
- Action: Test with degenerate mesh causing std::bad_alloc

---

## Recommendation for v1.0.0 Release

**Ship with Python-first positioning:**
- Python-only tests pass (876/877 + equivalence suite)
- C++ backend marked experimental until A1/A5 verified
- All Python-level mitigations deployed (A2, A3, A6)
- Rust backend source in sdist for rebuild capability

**Update README:**
```markdown
### Backends
- **Python** (production-ready, no dependencies)
- **C++** (experimental, requires pybind11/CMake rebuild)
- **Rust** (experimental, requires Cargo rebuild)

Force Python backend: `export CHILMESH_BACKEND=python`
```

---

**Threat Model Complete.** All findings documented. Ready for v1.0.0 release with Python-first policy pending C++ verification in CI.

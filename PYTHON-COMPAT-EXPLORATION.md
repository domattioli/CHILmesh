# Python Compatibility Exploration

**Issue:** #101 (Explore Python 3.8+ backwards compatibility)  
**Current:** Requires Python 3.10+  
**Target:** Investigate 3.8+ support  
**Effort:** 2-3 hours (exploration + proof-of-concept)

---

## Current State

**Minimum version:** Python 3.10+ (enforced in pyproject.toml)

**Blockers to 3.8+ support:**

### 1. Type Hints (Moderate Impact)
**Current code:**
```python
from typing import List, Tuple, Optional as Opt, Dict, Set, Union, Any

def __init__(
    self,
    connectivity: np.ndarray,
    points: np.ndarray,
    ...
) -> None:
```

**Issue:** Python 3.9 doesn't support `list[int]` (lowercase); must use `List[int]`.  
**Current code:** Already uses `List[int]` — **compatible with 3.8+** ✅

**But:** No `from __future__ import annotations` for lazy evaluation.  
**Implication:** Type hints evaluated at import time; circular imports possible (unlikely in CHILmesh).

**Cost to add lazy eval:** 1 line per file
```python
from __future__ import annotations  # Add to top of all modules
```

### 2. Optional Syntax
**3.10+ allows:** `x: int | None`  
**Current code:** Uses `Optional[int]` — **compatible with 3.8+** ✅

**Result:** No changes needed.

### 3. Match Statements
**3.10+:** `match x: case ...`  
**Current code:** No match statements found ✅

### 4. Dictionary Union Merging
**3.9+:** `dict1 | dict2`  
**Current code:** Grep found none ✅

### 5. Walrus Operator (`:=`)
**3.8+:** Supported in CHILmesh codebase  
**Current code:** Grep found none ✅

### 6. `@dataclass` Decorator
**3.7+:** Available  
**Current code:** Using regular classes ✅

### 7. `functools.cached_property`
**3.8+:** Available  
**Current code:** Not used ✅

### 8. `typing.TypedDict`
**3.8+:** Available  
**Current code:** Not used ✅

---

## Dependency Compatibility

| Dependency | 3.8+ Support | Notes |
|------------|---------|-------|
| numpy | ✅ Yes (1.19+) | 3.8+ support solid |
| scipy | ✅ Yes (1.5+) | 3.8+ support solid |
| matplotlib | ✅ Yes (3.2+) | 3.8+ support solid |
| pytest | ✅ Yes | Current version 7.0+ |

**Result:** All runtime + test deps support Python 3.8+ ✅

---

## Code Changes Required for 3.8+ Support

### Minimal Path (1-2 hours)

1. **Update pyproject.toml:**
   ```toml
   requires-python = ">=3.8"
   classifiers = [
       "Programming Language :: Python :: 3.8",
       "Programming Language :: Python :: 3.9",
       "Programming Language :: Python :: 3.10",
       "Programming Language :: Python :: 3.11",
       "Programming Language :: Python :: 3.12",
   ]
   ```

2. **Add `from __future__ import annotations`** to all modules (10 files):
   ```python
   from __future__ import annotations
   ```

3. **CI:** Add 3.8, 3.9 to test matrix in `.github/workflows/` ✅

**Rationale:** Current code is already 3.8+ compatible; just need to declare it + verify CI.

### Medium Path (3-4 hours) — Future Optimization

1. Use `|` syntax after 3.10+ min version (deferred)
2. Use `@cache` decorator (deferred)

---

## Risk Assessment

**Breaking changes:** 0 (only adding compatibility, not removing)  
**User impact:** Minimal — older Python installations now supported  
**Testing effort:** Low — CI matrix extended to 3.8, 3.9; all existing tests reused  
**Downstream impact:** MADMESHR, ADMESH may benefit from broader 3.8+ support

---

## Proof-of-Concept

### Step 1: Verify current imports parse in 3.8+
```bash
python3.8 -c "import ast; ast.parse(open('src/chilmesh/CHILmesh.py').read())"
```
**Expected:** No SyntaxError (3.10+ syntax would fail here)

### Step 2: Add `from __future__ import annotations`
```python
# src/chilmesh/CHILmesh.py (first line after docstring)
from __future__ import annotations
```

### Step 3: CI test
Run pytest on Python 3.8, 3.9, 3.10, 3.11, 3.12 in GitHub Actions.

---

## Recommendation

**Timeline:** Post-v1.0 (Phase 4.5 or Phase 5)  
**Effort:** 2-3 hours (minimal code change, CI setup)  
**Benefit:** Medium (broader user base; older research environments often stuck on 3.8-3.9)  
**Risk:** Very low (no code changes, only adding support)

**Why wait:**
- v1.0 release imminent; avoid late-stage changes
- No user demand yet (pyproject.toml doesn't restrict; users can try)
- CI setup more complex than code (worth batch with other 3.x testing)

**Trigger for implementation:**
- User request for 3.8+ support
- MADMESHR/ADMESH identify 3.8+ requirement
- Release cycle allows (post-v1.0)

---

## Checklist for Future Implementation

- [ ] Update `requires-python` in pyproject.toml to ">=3.8"
- [ ] Add `from __future__ import annotations` to all 10 modules
- [ ] Add 3.8, 3.9 to GitHub Actions test matrix
- [ ] Test locally: `python3.8 -m pytest` (if available)
- [ ] Update classifiers in pyproject.toml
- [ ] Verify all downstream consumers (MADMESHR, ADMESH) compatible
- [ ] Document in CHANGELOG: "Added support for Python 3.8, 3.9"
- [ ] Release as v0.4.0 or v1.1.0 (minor bump)

---

## Related

- Issue #101: Python 3.8+ compatibility exploration (this doc)
- Constitution Principle VII: API stability + semantic versioning
- Phase 4 completion: v0.2.0 shipped 2026-04-27

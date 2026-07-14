# Plan: Prebuilt C++ Binary Wheels on PyPI

**Date:** 2026-07-14
**Status:** Plan — not yet implemented
**Tracking:** [#229](https://github.com/domattioli/CHILmesh/issues/229) (PyPI ships pure-Python), [#225](https://github.com/domattioli/CHILmesh/issues/225) (macOS/Windows runner billing), [#234](https://github.com/domattioli/CHILmesh/issues/234)
**Related:** [`docs/RUST_EVALUATION.md`](../RUST_EVALUATION.md) (why C++, not Rust, is the wheel to ship)

## Goal — the consumer outcome

Make the C++ speedup the **default** for a plain install, with **no toolchain on the
user's side**:

```bash
pip install chilmesh[cpp]     # downloads a prebuilt binary; no compiler needed
```

Today `pip install chilmesh` is pure-Python (the compiled extension must be built
from source — a C++ compiler + CMake barrier that most users, especially on
managed/HPC/Windows machines, will not clear). Prebuilding the wheel once per
platform in CI and publishing it to PyPI removes that barrier entirely: `pip`
downloads the already-compiled binary matching the user's OS + Python.

**Non-goal:** Rust wheels. The Rust backend is **frozen** (see
[`RUST_EVALUATION.md`](../RUST_EVALUATION.md)); it earns no perf niche over C++ and
is not part of the distribution story.

## Current state (accurate as of this plan)

| Piece | State |
|---|---|
| `chilmesh` (main pkg) | pure-Python, setuptools, `requires-python >=3.10`, PyPI = no compiled ext (#229) |
| `chilmesh-cpp` (`src/chilmesh_cpp`) | separate package, scikit-build-core + pybind11, version `0.6.0.dev0`, **not on PyPI** |
| `build-cpp-wheels.yml` | **build-only**: cibuildwheel manylinux **x86_64 only**, `cp310–312`, musllinux skipped, asserts `import chilmesh_cpp` — uploads **artifacts only, no publish**; manual dispatch |
| `publish-pypi.yml` | publishes the **main** pkg on `release`; `twine` + `PYPI_API_TOKEN` secret |
| Backend selection | `chilmesh.backend_info()` auto-detects `import chilmesh_cpp` and picks it — **no code change needed** once the wheel is installable |

So the runtime detection already works; the entire gap is **packaging + distribution**.

## Design decision — distribution model

**Model A (recommended): publish `chilmesh-cpp` as its own binary-wheel package; add a
`chilmesh[cpp]` extra that depends on it.**

- `chilmesh` stays pure-Python (universal default + fallback; sdist/uncovered
  platforms still work).
- `chilmesh[cpp]` pulls the prebuilt `chilmesh-cpp` binary → auto-selected at import.
- Minimal churn: matches the existing two-package layout and the existing
  `import chilmesh_cpp` detection. No change to `chilmesh`'s build backend.

**Model B (rejected for now): bundle the C++ extension into the main `chilmesh` wheel.**
Would turn `chilmesh` itself into a per-platform binary package (change its build
backend to scikit-build-core, ship a pure-Python sdist fallback). More invasive, and
it couples the pure-Python reference release cadence to the compiled build matrix.
Revisit only if Model A's two-package version-sync proves painful.

## Phased plan

### Phase 0 — decide + pin (design, no CI yet)
- Ratify Model A.
- Define the **version-compatibility contract** between `chilmesh` and `chilmesh-cpp`
  (they version independently today: `2.0.0` vs `0.6.0.dev0`). Pin a compatible range
  in the extra, e.g. `chilmesh[cpp]` → `chilmesh-cpp>=0.6,<0.7`, and bump them together
  on any extension-API change. Record the contract in `CONTRIBUTING.md`.
- Graduate `chilmesh-cpp` off `.dev0` to a real release version (e.g. `0.6.0`).

### Phase 1 — expand the build matrix (`build-cpp-wheels.yml`)
- Add **macOS** (`x86_64` + `arm64` / universal2) and **Windows** (`AMD64`) alongside
  manylinux `x86_64`. Resolve the runner-billing question in **#225** first (macOS/
  Windows are paid GitHub runners).
- Keep `cp310 cp311 cp312`; add `cp313` once pybind11/deps support it. (pybind11 needs a
  wheel **per Python version** — this matrix is inherently `platforms × versions`. This
  is the one place Rust's maturin `abi3` single-wheel would be simpler — noted, but out
  of scope since C++ wins on perf.)
- Decide musllinux (currently skipped) and linux `aarch64` (emulated builds are slow) —
  defer both unless a consumer needs them.
- Strengthen `CIBW_TEST_COMMAND`: import **and** run one real `full_init` on a tiny mesh,
  not just `hasattr(full_init)`.

### Phase 2 — add a publish job
- New job (or new workflow `publish-cpp-wheels.yml`) that, on a `chilmesh-cpp` **release
  tag**, builds the full matrix + sdist and `twine upload`s to PyPI.
- Auth: reuse the `PYPI_API_TOKEN` pattern from `publish-pypi.yml`, **or** migrate to
  PyPI **Trusted Publishing** (OIDC, `id-token: write`) — preferred, no long-lived token.
- Guard strictly on the tag/release event so nothing publishes from `development`.

### Phase 3 — wire the consumer-facing extra
- In `chilmesh`'s `pyproject.toml`:
  ```toml
  [project.optional-dependencies]
  cpp = ["chilmesh-cpp>=0.6,<0.7"]
  ```
- No backend code change — `backend_info()` already imports `chilmesh_cpp`.
- README + docs: document `pip install chilmesh[cpp]`; update the "PyPI installs are
  pure-Python" note (#229) to "plain install = pure-Python; `chilmesh[cpp]` = prebuilt
  C++ where a wheel exists, pure-Python fallback elsewhere."

### Phase 4 — validate + close out
- Post-publish smoke: in a **clean** environment on each platform,
  `pip install chilmesh[cpp]` then assert `backend_info()['selected'] == 'cpp'` and run
  `tests/test_backend_equivalence.py` against the **installed wheel** (not `src/` — avoid
  the namespace-stub shadowing from #163: run from outside the repo tree).
- Update `docs/BENCHMARK.md` / README wheel-availability notes; close **#229**.

## Risks & considerations
- **Runner cost (#225):** macOS/Windows runners are billed; Phase 1 is gated on that
  decision. manylinux (Linux) is free on `ubuntu-latest` and can ship first.
- **Matrix size:** pybind11 → one wheel per (platform × Python version). Manageable at
  3–4 Python versions × 3 platforms, but it grows.
- **macOS arm64:** cross-build/`universal2` needs care in cibuildwheel; test on Apple
  silicon if possible.
- **glibc target:** pick a manylinux baseline (e.g. `manylinux_2_28`) broad enough for
  HPC/older distros.
- **Version sync:** the `chilmesh` ↔ `chilmesh-cpp` compatibility pin must be maintained
  on every extension-API change (Phase 0 contract).
- **Test isolation:** always validate against the installed wheel, never the `src/`
  source dir on `sys.path` (#163 false-positive `CPP_AVAILABLE`).

## Sequencing
1. Phase 0 (design/pin) — cheap, unblocks everything.
2. Phase 1 Linux-only publish first (free runners) → real prebuilt wheel on PyPI for the
   largest user base, fastest win.
3. Resolve #225, then add macOS/Windows (Phase 1 remainder) + Phase 2 publish job.
4. Phase 3 extra + docs, Phase 4 validation, close #229.

---

## Appendix — ready-to-apply snippets

> **PROPOSED — not yet applied.** These are copy-paste-ready but intentionally
> *inert* (they are documentation, not live workflow files). Apply them only after
> Model A is ratified (Phase 0) and the macOS/Windows runner-billing question
> (#225) is resolved. The `macos-latest` / `windows-latest` matrix legs and the
> publish job are the parts that incur runner cost or perform the **irreversible,
> outward** act of uploading to PyPI — those need operator sign-off before they go
> live. The Linux-only build leg is free and safe to enable first.

### A1 — Phase 1: expand `build-cpp-wheels.yml` (matrix + real smoke test)

Replace the single `build-manylinux` job with a matrixed build. macOS/Windows legs
are gated on #225 — drop them from `matrix.os` to ship Linux-only first.

```yaml
jobs:
  build-wheels:
    name: cibuildwheel chilmesh_cpp (${{ matrix.os }})
    runs-on: ${{ matrix.os }}
    timeout-minutes: 45
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]  # macos/windows gated on #225
    steps:
      - uses: actions/checkout@v5
      - uses: actions/setup-python@v6
        with:
          python-version: "3.12"
      - name: Install cibuildwheel
        run: python -m pip install cibuildwheel==2.21.3
      - name: Build wheels
        run: python -m cibuildwheel --output-dir wheelhouse src/chilmesh_cpp
        env:
          CIBW_BUILD: "cp310-* cp311-* cp312-*"
          CIBW_ARCHS_LINUX: "x86_64"
          CIBW_ARCHS_MACOS: "x86_64 arm64"
          CIBW_ARCHS_WINDOWS: "AMD64"
          CIBW_SKIP: "*-musllinux*"
          # Real smoke: build a mesh, not just check the symbol exists.
          CIBW_TEST_REQUIRES: "numpy"
          CIBW_TEST_COMMAND: >
            python -c "import numpy as np, chilmesh_cpp;
            pts=np.array([[0.,0.],[1.,0.],[0.,1.],[1.,1.]]);
            conn=np.array([[0,1,2],[1,3,2]],dtype=np.int32);
            m=chilmesh_cpp.full_init(pts, conn);
            assert m.n_elems==2, m.n_elems; print('cpp full_init OK', m.n_verts, m.n_elems)"
      - uses: actions/upload-artifact@v4
        with:
          name: chilmesh-cpp-wheels-${{ matrix.os }}
          path: wheelhouse/*.whl
          if-no-files-found: error
```

### A2 — Phase 2: `publish-cpp-wheels.yml` (tag/release-gated publish)

New workflow. Prefer **Trusted Publishing** (OIDC, no long-lived token); the
token variant (mirroring `publish-pypi.yml`) is shown as a fallback comment. The
`build` job reuses A1; only `publish` is new.

```yaml
name: publish-cpp-wheels
on:
  release:
    types: [published]     # publishes ONLY on a real GitHub release; never on push
  workflow_dispatch:
permissions:
  contents: read
jobs:
  build:
    # ... reuse the A1 matrix build; uploads per-OS wheel artifacts ...
  sdist:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v5
      - run: pipx run build --sdist --outdir dist src/chilmesh_cpp
      - uses: actions/upload-artifact@v4
        with: { name: chilmesh-cpp-sdist, path: dist/*.tar.gz }
  publish:
    needs: [build, sdist]
    runs-on: ubuntu-latest
    environment: pypi              # protect with required reviewers in repo settings
    permissions:
      id-token: write              # PyPI Trusted Publishing (OIDC) — no secret
    steps:
      - uses: actions/download-artifact@v4
        with: { path: dist, merge-multiple: true }
      - uses: pypa/gh-action-pypi-publish@release/v1
        # Token fallback (if not using Trusted Publishing):
        #   with: { password: ${{ secrets.PYPI_API_TOKEN }} }
```

### A3 — Phase 3: consumer-facing extra + docs

`chilmesh/pyproject.toml` (add alongside the existing `dev` extra) — **land this in
the same change that publishes `chilmesh-cpp`, never before**, or
`pip install chilmesh[cpp]` resolves to a package that isn't on PyPI yet:

```toml
[project.optional-dependencies]
cpp = ["chilmesh-cpp>=0.6,<0.7"]   # pin bumps with any extension-API change (Phase 0 contract)
```

README / docs one-liner:

```bash
pip install chilmesh          # pure-Python, runs everywhere
pip install chilmesh[cpp]     # + prebuilt C++ acceleration where a wheel exists
```

No backend code changes — `chilmesh.backend_info()` already imports `chilmesh_cpp`
and auto-selects it.

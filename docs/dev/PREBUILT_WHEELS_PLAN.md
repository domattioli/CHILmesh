# Plan: Prebuilt C++ Binary Wheels on PyPI

**Date:** 2026-07-14
**Status:** Plan, not yet implemented
**Tracking:** [#229](https://github.com/domattioli/CHILmesh/issues/229) (PyPI ships pure-Python), [#225](https://github.com/domattioli/CHILmesh/issues/225) (macOS/Windows runner billing), [#234](https://github.com/domattioli/CHILmesh/issues/234)
**Related:** [`docs/RUST_EVALUATION.md`](../RUST_EVALUATION.md) (why C++, not Rust, is the wheel to ship)

## Summary

Publish `chilmesh-cpp` on PyPI as its own binary-wheel package (Model A). Add a `chilmesh[cpp]` extra that depends on it.
A plain install then gets the C++ speedup by default, with no compiler on the user's machine.
Runtime detection already works, so the remaining gap is packaging and distribution.

```bash
pip install chilmesh[cpp]     # downloads a prebuilt binary; no compiler needed
```

## Goal

Today `pip install chilmesh` is pure-Python.
The compiled extension must be built from source.
That needs a C++ compiler and CMake, which many users cannot get, especially on managed, HPC, or Windows machines.

CI prebuilds the wheel once per platform and publishes it to PyPI.
Then `pip` downloads the compiled binary that matches the user's OS and Python.

**Non-goal:** Rust wheels.
The Rust backend is frozen (see [`RUST_EVALUATION.md`](../RUST_EVALUATION.md)).
It has no performance niche over C++, so it is outside the distribution plan.

## Current state

| Piece | State |
|---|---|
| `chilmesh` (main pkg) | pure-Python, setuptools, `requires-python >=3.10`, PyPI = no compiled ext (#229) |
| `chilmesh-cpp` (`src/chilmesh_cpp`) | separate package, scikit-build-core + pybind11, version `0.6.0.dev0`, **not on PyPI** |
| `build-cpp-wheels.yml` | **build-only**: cibuildwheel manylinux **x86_64 only**, `cp310–312`, musllinux skipped, asserts `import chilmesh_cpp`; uploads **artifacts only, no publish**; manual dispatch |
| `publish-pypi.yml` | publishes the **main** pkg on `release`; `twine` + `PYPI_API_TOKEN` secret |
| Backend selection | `chilmesh.backend_info()` auto-detects `import chilmesh_cpp` and picks it; **no code change needed** once the wheel is installable |

## Design decision: distribution model

**Model A (recommended):** publish `chilmesh-cpp` as its own binary-wheel package, and add a `chilmesh[cpp]` extra that depends on it.

- `chilmesh` stays pure-Python. It is the universal default and fallback, so the sdist and uncovered platforms still work.
- `chilmesh[cpp]` pulls the prebuilt `chilmesh-cpp` binary. The backend is auto-selected at import.
- Churn is minimal. The layout matches the existing two-package structure and the existing `import chilmesh_cpp` detection.
- The build backend of `chilmesh` does not change.

**Model B (rejected for now):** bundle the C++ extension into the main `chilmesh` wheel.

- `chilmesh` itself would become a per-platform binary package.
- Its build backend would change to scikit-build-core, with a pure-Python sdist fallback.
- It is more invasive.
- It couples the release cadence of the pure-Python reference to the compiled build matrix.
- Revisit it only if the two-package version sync of Model A proves painful.

## Phased plan

### Phase 0: decide and pin (design, no CI yet)

1. Ratify Model A.
2. Define the version-compatibility contract between `chilmesh` and `chilmesh-cpp`. They version independently today: `chilmesh` 1.4.x and `chilmesh-cpp` `0.6.0.dev0`.
3. Pin a compatible range in the extra, for example `chilmesh[cpp]` → `chilmesh-cpp>=0.6,<0.7`.
4. Bump both versions together on any extension-API change.
5. Record the contract in `CONTRIBUTING.md`.
6. Graduate `chilmesh-cpp` off `.dev0` to a real release version (for example `0.6.0`).

**Phase 0 decisions (2026-10-02, operator-ratified for #256):**

- Model A ratified.
- `chilmesh-cpp` graduated `0.6.0.dev0` to `0.6.0`; its `requires-python` raised to `>=3.10` to match `chilmesh`.
- **Compat pin:** `chilmesh` 1.5.x pairs with `chilmesh-cpp>=0.6,<0.7`, enforced by the `cpp` extra in the root `pyproject.toml`. Bump the pin and both versions together on any extension-API change. Release target: `chilmesh` 1.5.0 and `chilmesh-cpp` 0.6.0 on one GitHub release.
- Platforms: Linux x86_64, macOS x86_64 + arm64, Windows AMD64; CPython 3.10 to 3.14 (cibuildwheel 3.3.0; 3.1.0 made cp314 a default build target).
- The Intel macOS leg uses `macos-15-intel` because `macos-13` was retired by GitHub.

### Phase 1: expand the build matrix (`build-cpp-wheels.yml`)

1. Resolve the runner-billing question in #225 first. macOS and Windows are paid GitHub runners.
2. Add macOS (`x86_64` + `arm64` / universal2) and Windows (`AMD64`) next to manylinux `x86_64`.
3. Keep `cp310 cp311 cp312`. Add `cp313` once pybind11 and the dependencies support it.
4. Defer musllinux (currently skipped) and linux `aarch64` (emulated builds are slow) unless a consumer needs them.
5. Strengthen `CIBW_TEST_COMMAND`: import the module and run one real `full_init` on a tiny mesh, not just `hasattr(full_init)`.

pybind11 needs one wheel per Python version, so the matrix is `platforms × versions`.
Rust's maturin `abi3` single wheel would be simpler here.
That is out of scope, since C++ wins on performance.

### Phase 2: add a publish job

1. Add a new job, or a new workflow `publish-cpp-wheels.yml`.
2. On a `chilmesh-cpp` release tag, build the full matrix and the sdist, then `twine upload` to PyPI.
3. Authenticate with PyPI Trusted Publishing (OIDC, `id-token: write`). This is preferred because it needs no long-lived token.
4. As an alternative, reuse the `PYPI_API_TOKEN` pattern from `publish-pypi.yml`.
5. Guard strictly on the tag or release event, so nothing publishes from `development`.

### Phase 3: wire the consumer-facing extra

1. In the `pyproject.toml` of `chilmesh`, add the extra:

  ```toml
  [project.optional-dependencies]
  cpp = ["chilmesh-cpp>=0.6,<0.7"]
  ```

2. Make no backend code change. `backend_info()` already imports `chilmesh_cpp`.
3. Document `pip install chilmesh[cpp]` in the README and docs.
4. Replace the "PyPI installs are pure-Python" note (#229) with this: plain install = pure-Python; `chilmesh[cpp]` = prebuilt C++ where a wheel exists, pure-Python fallback elsewhere.

### Phase 4: validate and close out

1. In a clean environment on each platform, run `pip install chilmesh[cpp]`.
2. Assert `backend_info()['selected'] == 'cpp'`.
3. Run `tests/test_backend_equivalence.py` against the installed wheel, not `src/`.
4. Run from outside the repo tree to avoid the namespace-stub shadowing from #163.
5. Update the wheel-availability notes in `docs/BENCHMARK.md` and the README.
6. Close #229.

## Risks and considerations

- **Runner cost (#225):** macOS and Windows runners are billed, so Phase 1 is gated on that decision. manylinux (Linux) is free on `ubuntu-latest` and can ship first.
- **Matrix size:** pybind11 gives one wheel per (platform × Python version). That is manageable at 3–4 Python versions × 3 platforms, but it grows.
- **macOS arm64:** cross-build or `universal2` needs care in cibuildwheel. Test on Apple silicon if possible.
- **glibc target:** pick a manylinux baseline (for example `manylinux_2_28`) that is broad enough for HPC and older distros.
- **Version sync:** maintain the `chilmesh` ↔ `chilmesh-cpp` compatibility pin on every extension-API change (Phase 0 contract).
- **Test isolation:** always validate against the installed wheel, never the `src/` source dir on `sys.path` (#163 false-positive `CPP_AVAILABLE`).

## Sequencing

1. Phase 0 (design and pin). It is cheap and unblocks everything.
2. Phase 1, Linux-only publish first (free runners). This puts a real prebuilt wheel on PyPI for the largest user base, the fastest win.
3. Resolve #225, then add macOS and Windows (Phase 1 remainder) and the Phase 2 publish job.
4. Phase 3 extra and docs, then Phase 4 validation, then close #229.

---

## Appendix: ready-to-apply snippets

> **PROPOSED: not yet applied.** These are copy-paste-ready but intentionally
> *inert* (they are documentation, not live workflow files). Apply them only after
> Model A is ratified (Phase 0) and the macOS/Windows runner-billing question
> (#225) is resolved. The `macos-latest` / `windows-latest` matrix legs and the
> publish job are the parts that incur runner cost or perform the **irreversible,
> outward** act of uploading to PyPI; those need operator sign-off before they go
> live. The Linux-only build leg is free and safe to enable first.

### A1: Phase 1: expand `build-cpp-wheels.yml` (matrix + real smoke test)

Replace the single `build-manylinux` job with a matrixed build. macOS/Windows legs
are gated on #225: drop them from `matrix.os` to ship Linux-only first.

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

### A2: Phase 2: `publish-cpp-wheels.yml` (tag/release-gated publish)

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

### A3: Phase 3: consumer-facing extra + docs

`chilmesh/pyproject.toml` (add alongside the existing `dev` extra), **land this in
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

No backend code changes: `chilmesh.backend_info()` already imports `chilmesh_cpp`
and auto-selects it.

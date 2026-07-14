# chilmesh_core (Rust backend) — STATUS: FROZEN

**Status:** ❄️ **FROZEN** (2026-07-14, operator-directed)
**Evaluation:** [`docs/RUST_EVALUATION.md`](../../docs/RUST_EVALUATION.md)

## What "frozen" means here

This Rust quad-edge backend is **kept but not developed further**:

- **Kept** — the crate, its `rust-equivalence` CI job, the cross-backend equivalence
  tests, and the measured benchmark data all stay. The backend remains buildable
  (`maturin build --release`) and output-equivalent to Python (76/76 equivalence tests).
- **Not evolved** — no new features, and it is **not** kept in lockstep as the Python
  API grows. Treat it as a reference/experimental artifact, not a maintained backend.

## Why (summary; full detail + numbers in the evaluation)

A measured, like-for-like comparison found that Rust:

- full-inits ~3–5× faster than pure Python but **~2–5× slower than C++** (≈5× behind on
  Block_O) — it earns **no performance niche over C++**, which is the recommended,
  bit-identical accelerator;
- adds a second compiled-backend maintenance surface (a 2nd toolchain, a 2nd wheel
  matrix, a 2nd equivalence lane, a 2nd binding set) for no speed C++ doesn't already
  provide.

The one real query-path defect (`get_vertex_edges` was O(_n_)/call) was fixed before the
freeze — queries are now O(1) — so the frozen state is *correct and fast on queries*, just
not the acceleration path. See the evaluation for the before/after numbers and the
"could Rust replace Python as the default?" analysis (answer: no — ship prebuilt **C++**
wheels instead; see [`docs/dev/PREBUILT_WHEELS_PLAN.md`](../../docs/dev/PREBUILT_WHEELS_PLAN.md)).

## If you are tempted to revive it

Only two greenfield niches could justify a compiled backend that C++ doesn't already
own — Phase-5 spatial indexing and topology-mutation `splice`. Even there, extend the
existing **C++ half-edge** backend first; unfreezing Rust needs a fresh, measured case.

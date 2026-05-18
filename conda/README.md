# Conda-forge recipe

Submission scaffold for [conda-forge](https://conda-forge.org/). Canonical
recipe lives at
[`conda-forge/chilmesh-feedstock`](https://github.com/conda-forge/staged-recipes)
once accepted; this directory tracks the in-repo source-of-truth.

## First-time submission

1. Publish the matching release to PyPI (`/pypi-publish`).
2. Fetch the sdist sha256:
   ```bash
   curl -sL https://pypi.org/pypi/chilmesh/json \
     | python -c "import sys,json; print(json.load(sys.stdin)['releases']['0.4.0'][0]['digests']['sha256'])"
   ```
3. Paste the hash into `meta.yaml` under `source.sha256`.
4. Fork [`conda-forge/staged-recipes`](https://github.com/conda-forge/staged-recipes).
5. Copy `meta.yaml` to `recipes/chilmesh/meta.yaml` in the fork.
6. Open a PR. Conda-forge bots run CI; maintainers approve.

## Updates after acceptance

Once accepted, a `chilmesh-feedstock` repo is auto-created. Future updates
land there via either:

- **Auto-bot.** `regro-cf-autotick-bot` opens an update PR within hours of a
  new PyPI release. Approve + merge.
- **Manual.** Edit `recipe/meta.yaml` in the feedstock; bump `version` and
  refresh `sha256`.

## Dependency notes

- Uses `matplotlib-base` (not `matplotlib`) — matches conda-forge convention
  to avoid pulling in GUI backend deps unless requested.
- `python >=3.8` mirrors `pyproject.toml`.

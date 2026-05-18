# Zenodo DOI setup

One-time auth + per-release auto-archive workflow. After setup, every
GitHub release auto-deposits to Zenodo and mints a fresh DOI.

## One-time setup (5 minutes)

1. **Sign in to Zenodo with GitHub.** Visit
   <https://zenodo.org/account/settings/github/> and authorise the OAuth
   app.
2. **Enable archiving for `domattioli/CHILmesh`.** Toggle the switch next
   to the repo. Zenodo subscribes to GitHub's release webhook.
3. **Cut a fresh GitHub release** (`/github-release`) — this is the first
   one Zenodo will see. The archive happens automatically; check
   <https://zenodo.org/account/settings/github/> after a few minutes.
4. **Copy the two DOIs** from the Zenodo upload page:
   - **Concept DOI** — always resolves to the *latest* release; cite this
     when you want a stable evergreen link.
   - **Version DOI** — pins this specific version (v0.4.0); cite this
     in papers that want reproducibility.

## Wire DOIs back into the repo

5. Edit `CITATION.cff` — uncomment the `identifiers:` block and paste
   both DOIs.
6. Add a Zenodo badge to `README.md` (top, near the PyPI badge):

   ```markdown
   [![DOI](https://zenodo.org/badge/<REPO_ID>.svg)](https://zenodo.org/badge/latestdoi/<REPO_ID>)
   ```

   Zenodo prints the exact markdown on the upload page; just paste.
7. Commit + push the DOI metadata to `daily-issue-fixing`. Open a PR.

## Per-release maintenance

- `/github-release` is the only manual step. Zenodo runs on the webhook.
- Bump `version` and `date-released` in `CITATION.cff` *before* cutting
  the release so the archived metadata is correct.
- The Concept DOI never changes; only Version DOI rotates per release.

## Notes

- Zenodo archives the GitHub release's source tarball, not the PyPI
  wheel. Anyone with the DOI can re-download the exact code, even if
  the GitHub repo disappears.
- JOSS submissions reference the Concept DOI in the review form.

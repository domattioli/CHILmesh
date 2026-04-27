# Decision Tree Templates for instructions.sh

Reference templates for the `{WORKFLOW_DECISION_TREE}` placeholder in `instructions.sh`.
Select the template matching the repo's workflow type. Adapt specifics (branch names,
test commands, target paths) to the actual repo — don't use placeholders verbatim.

---

## autonomous-routine

Use for: long-running training loops, data pipelines, mesh generation, RL experiments,
any session where a human should never need to intervene.

```markdown
## Decision tree

On session start:
1. Read CLAUDE.md and .specify/memory/constitution.md (if present).
2. Check git status — if dirty with unstaged changes from a prior run, commit or stash
   before proceeding. Never silently lose work.
3. Run the primary task (defined below).
4. After each major step, emit a one-line status: what ran, what was produced/changed.
5. At session end, commit any changes, emit session summary (cost, steps, outcomes).

Primary task loop:
- IF there is an explicit task spec (spec.md / plan.md in specs/), read it first.
- ELSE fall back to the repo-level mission in CLAUDE.md.
- Work through the task list top-to-bottom. Check off each item as complete.
- Do NOT ask for human input. If a decision is ambiguous and not covered by the
  constitution, pick the more conservative option and document the choice.
- Checkpoint state (commit + push) after every N=5 completed items or 30 min of
  wall time, whichever comes first.

On unexpected error:
- Log the error with full context.
- If recoverable (network timeout, transient file lock): retry once, then continue.
- If unrecoverable (missing required file, hard stop triggered): stop, commit partial
  work with a PARTIAL: prefix, emit summary, exit.
```

---

## issue-triage

Use for: reviewing, labeling, commenting on, or closing GitHub issues autonomously.

```markdown
## Decision tree

On session start:
1. Read CLAUDE.md and .specify/memory/constitution.md (if present).
2. Authenticate gh CLI: verify `gh auth status` returns OK. If not, stop and report.
3. Fetch open issues: `gh issue list --state open --limit 20`.
4. For each issue (up to 20 per run):
   a. Read the full issue body and all comments.
   b. Classify: bug / feature-request / question / duplicate / stale.
   c. Apply label matching the classification.
   d. If duplicate: comment with the canonical issue number and close.
   e. If stale (no activity > 60 days): comment "Closing as stale — reopen if
      still relevant." and close.
   f. If question and answer is known: post the answer as a comment.
   g. Otherwise: leave as open, add triage label, add brief comment summarising
      what investigation is needed.
5. Emit session summary: N issues triaged, N closed, N labeled, N commented.

Never:
- Close an issue that has the `human-review` label.
- Add or remove the `human-review` label.
- Close more than 10 issues per run without human review step.
```

---

## data-curation

Use for: ETL pipelines, data quality checks, sync routines, file transformation.

```markdown
## Decision tree

On session start:
1. Read CLAUDE.md and .specify/memory/constitution.md (if present).
2. Confirm source data is accessible (check paths / API endpoints). If not, stop.
3. Run the ETL / curation task:
   a. Read source data.
   b. Apply transforms / quality checks.
   c. Write output to the STAGING location (never overwrite source directly).
   d. Validate output: row count sanity check, schema check, null-rate check.
4. If validation passes: promote staging to production target.
5. If validation fails: leave staging in place, log failures, do NOT promote, stop.
6. Emit summary: rows in, rows out, rows failed validation, rows promoted.

Never:
- Delete or overwrite source data.
- Write directly to production without staging first.
- Silently drop rows — log every dropped row with reason.
```

---

## deployment

Use for: publishing to PyPI/npm, creating GitHub releases, deploying to cloud.

```markdown
## Decision tree

On session start:
1. Read CLAUDE.md and .specify/memory/constitution.md (if present).
2. Confirm CI is green: `gh run list --limit 5` — check the most recent run status.
   If any required checks are failing, stop. Do NOT deploy on red.
3. Confirm version bump is present:
   - For Python: check pyproject.toml / setup.cfg version > last published version.
   - For npm: check package.json version > last published version.
   If no version bump, stop and log.
4. Run tests locally: `{TEST_CMD}`. If any tests fail, stop.
5. Build the distribution artifact.
6. Publish: `{PUBLISH_CMD}`.
7. Create GitHub release with changelog notes.
8. Emit summary: version published, release URL, time taken.

Never:
- Push directly to main without passing CI.
- Publish without a version bump.
- Overwrite an existing published version.
- Deploy on a branch other than `main` (unless explicitly configured).
```

---

## custom

Use for: user-defined workflows not covered by the above templates.

```markdown
## Decision tree

{USER_DEFINED_TASK_DESCRIPTION}

General principles:
1. Read CLAUDE.md and .specify/memory/constitution.md (if present) before starting.
2. Work top-to-bottom through the task. Do NOT ask for human input.
3. Pick the more conservative option when a decision is ambiguous.
4. Checkpoint (commit + push or write to disk) after every major step.
5. Emit a one-line status after each step and a full summary at session end.

On error:
- Recoverable: retry once, then continue.
- Unrecoverable: stop, save partial work, emit summary.
```

#!/bin/bash
# instructions_on_start.sh — On-start health check for CHILmesh (downstream consumer of DomI)
#
# USAGE:
#   bash "$(git rev-parse --show-toplevel)"/scripts/instructions_on_start.sh
#
# CONTRACT:
#   1. DomI drift check (HARD STOP if behind) — see "DomI Sync Contract" in CLAUDE.md
#   2. CLAUDE.md presence
#   3. Git hygiene (branch, dirty, stash)
#   4. Repo-local extension via scripts/onstart_local.sh (optional)
#
# ENV:
#   DOMI_OWNER  (default: domattioli)
#   DOMI_REPO   (default: DomI)
#   DOMI_BRANCH (default: main)
#   DOMI_BLOCK_ON_DRIFT (default: 1; set 0 for read-only sessions)

set -euo pipefail

# ============================================================================
# DomI drift check — pasted from the sync-from-domi skill template
# (~/.claude/plugins/cache/DomI/sync-from-domi/<ver>/skills/sync-from-domi/templates/downstream_startup_hook.sh)
# ============================================================================
# --- locate sync-from-domi ---
# Search order:
#   1. Plugin cache (after `claude plugin install sync-from-domi@DomI`)
#      Path pattern: ~/.claude/plugins/cache/DomI/sync-from-domi/<version>/skills/sync-from-domi
#   2. Plugin marketplace clone
#   3. Vendored copy in this repo
#   4. Global ~/.claude/skills
DOMI_SKILL_PATH="${DOMI_SKILL_PATH:-}"
if [ -z "$DOMI_SKILL_PATH" ]; then
  # Glob-match the versioned cache path first (preferred)
  for cached in "${HOME}"/.claude/plugins/cache/DomI/sync-from-domi/*/skills/sync-from-domi; do
    if [ -f "${cached}/scripts/check_pin.sh" ]; then
      DOMI_SKILL_PATH="$cached"
      break
    fi
  done
fi
if [ -z "$DOMI_SKILL_PATH" ]; then
  for candidate in \
    "${HOME}/.claude/plugins/marketplaces/domattioli/DomI/plugins/sync-from-domi/skills/sync-from-domi" \
    "./plugins/sync-from-domi/skills/sync-from-domi" \
    "./skills/sync-from-domi" \
    "${HOME}/.claude/skills/sync-from-domi"; do
    if [ -f "${candidate}/scripts/check_pin.sh" ]; then
      DOMI_SKILL_PATH="$candidate"
      break
    fi
  done
fi

if [ -z "$DOMI_SKILL_PATH" ] || [ ! -f "${DOMI_SKILL_PATH}/scripts/check_pin.sh" ]; then
  echo "⚠ sync-from-domi skill not found locally; skipping DomI drift check"
  echo "  → install via: claude plugin marketplace add domattioli/DomI && claude plugin install sync-from-domi@DomI"
else
  set +e
  bash "${DOMI_SKILL_PATH}/scripts/check_pin.sh"
  DOMI_DRIFT_RC=$?
  set -e

  case $DOMI_DRIFT_RC in
    0)
      : # synced; continue
      ;;
    1)
      echo ""
      echo "============================================================"
      echo "🛑 HARD STOP: downstream is BEHIND DomI"
      echo "============================================================"
      echo "Invoke the sync-from-domi skill before any write work:"
      echo "  > sync from DomI"
      echo "Or run manually:"
      echo "  bash ${DOMI_SKILL_PATH}/scripts/update_pin.sh"
      echo "  (then commit .domi-pin and any updated skills)"
      echo "============================================================"
      if [ "${DOMI_BLOCK_ON_DRIFT:-1}" = "1" ]; then
        exit 1
      fi
      ;;
    2)
      echo "ⓘ First-time DomI pin needed; will create on next sync"
      ;;
    3)
      echo ""
      echo "============================================================"
      echo "🛑 HARD STOP: DomI pin FORKED (manifest hash mismatch)"
      echo "============================================================"
      echo "Local edits to vendored DomI artifacts suspected."
      echo "Operator must resolve manually before continuing."
      echo "============================================================"
      exit 1
      ;;
    4)
      echo "⚠ DomI drift check skipped (gh unavailable); continuing"
      ;;
  esac
fi
# --- end DomI drift check ---

# ============================================================================
# Audit checks
# ============================================================================
REPO_ROOT="${REPO_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null || echo ".")}"
cd "$REPO_ROOT" || exit 1

START_TIME=$(date +%s)
ISSUES=0
BLOCKERS=()

echo "=== On-Start Health Check ==="
echo "Repo: $REPO_ROOT"
echo "Started: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ""

# 1. CLAUDE.md presence
echo "Checking CLAUDE.md..."
if [ ! -f "$REPO_ROOT/CLAUDE.md" ] && [ ! -f "$REPO_ROOT/.claude/CLAUDE.md" ]; then
  echo "  ❌ CLAUDE.md missing — every Claude-driven repo must have one"
  echo "     Bootstrap: /maintain-claude-md init"
  BLOCKERS+=("CLAUDE.md missing")
  ISSUES=$((ISSUES + 1))
else
  echo "  ✓ CLAUDE.md present"
fi
echo ""

# 2. Git hygiene
echo "Checking git hygiene..."
BRANCH="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")"
DIRTY="$(git status --porcelain 2>/dev/null | wc -l | tr -d ' ')"
STASH_COUNT="$(git stash list 2>/dev/null | wc -l | tr -d ' ')"
echo "  Branch: $BRANCH | Uncommitted: $DIRTY | Stashed: $STASH_COUNT"
if [ "$BRANCH" = "main" ] || [ "$BRANCH" = "master" ]; then
  echo "  ⚠ On default branch — work on a feature branch per CLAUDE.md"
fi
echo ""

# 3. Repo-local extension (optional)
if [ -f "$REPO_ROOT/scripts/onstart_local.sh" ]; then
  echo "Running repo-local extension (scripts/onstart_local.sh)..."
  if bash "$REPO_ROOT/scripts/onstart_local.sh"; then
    echo "  ✓ Local extension passed"
  else
    echo "  ❌ Local extension failed"
    ISSUES=$((ISSUES + 1))
  fi
  echo ""
fi

# ============================================================================
# Summary
# ============================================================================
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo "=== Summary ==="
if [ ${#BLOCKERS[@]} -gt 0 ]; then
  echo "Status: 🛑 BLOCKED"
  echo "Blockers:"
  printf "  - %s\n" "${BLOCKERS[@]}"
  echo "Runtime: ${DURATION}s"
  exit 1
elif [ $ISSUES -gt 0 ]; then
  echo "Status: ⚠ ISSUES ($ISSUES found)"
  echo "Runtime: ${DURATION}s"
  exit 0
else
  echo "Status: ✓ HEALTHY"
  echo "Runtime: ${DURATION}s"
  exit 0
fi

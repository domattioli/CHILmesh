#!/bin/bash
# scripts/instructions_on_start.sh — session startup health check
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || echo ".")" 
cd "$REPO_ROOT" || exit 1

echo "=== Session Start: CHILmesh ==="
echo "Branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null) | Dirty: $(git status --porcelain 2>/dev/null | wc -l | tr -d ' ') files"
echo ""

# DomI drift check (plugin cache → skills marketplace → vendored)
_find_check_pin() {
  for d in "$HOME/.claude/plugins/cache/DomI/sync-from-domi" \
            "$HOME/.claude/skills/sync-from-domi" \
            "$REPO_ROOT/plugins/sync-from-domi"; do
    local f; f=$(find "$d" -name "check_pin.sh" -maxdepth 5 2>/dev/null | head -1)
    [ -n "$f" ] && echo "$f" && return 0
  done; return 1
}

CHECK_PIN=$(_find_check_pin 2>/dev/null || true)
if [ -n "$CHECK_PIN" ]; then
  set +e; bash "$CHECK_PIN"; rc=$?; set -e
  case $rc in
    0) echo "✓ DomI pin current" ;;
    1|3) echo "HARD STOP: DomI drift (exit $rc). Run '/sync-from-domi' before write work." >&2; exit 1 ;;
    2) echo "⚠ .domi-pin absent — run update_pin.sh to initialize" ;;
    4) echo "⚠ gh unavailable — DomI drift check skipped" ;;
  esac
else
  echo "⚠ sync-from-domi not installed. Run: claude plugin install sync-from-domi@DomI"
fi
echo ""

echo "=== ✓ Health check passed ==="

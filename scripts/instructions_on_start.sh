#!/bin/bash
# scripts/instructions_on_start.sh — session startup health check
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || echo ".")" 
cd "$REPO_ROOT" || exit 1

echo "=== Session Start: CHILmesh ==="
echo "Branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null) | Dirty: $(git status --porcelain 2>/dev/null | wc -l | tr -d ' ') files"
echo ""

# Remote health check — auto-heal dead local git proxy by falling back to github.com
# Pattern: cloud sessions sometimes route `origin` through a local proxy (127.0.0.1:<port>)
# that dies mid-session. Detect dead proxy + rewrite remote when GITHUB_TOKEN available.
# See DomI #48 (push-via-mcp) for the workaround that motivated this.
_remote_url=$(git config --get remote.origin.url 2>/dev/null || echo "")
if [[ "$_remote_url" =~ ^http://(.+@)?127\.0\.0\.1:([0-9]+)/ ]]; then
  _port="${BASH_REMATCH[2]}"
  if ! (echo > /dev/tcp/127.0.0.1/"$_port") 2>/dev/null; then
    if [ -n "${GITHUB_TOKEN:-}" ]; then
      git remote set-url origin "https://x-access-token:${GITHUB_TOKEN}@github.com/domattioli/CHILmesh.git"
      echo "⚠ Local git proxy dead on :$_port — switched origin to github.com via GITHUB_TOKEN"
      echo "⚠ proxy dead — manual recovery: git fetch origin development && git status"
    else
      echo "⚠ Local git proxy dead on :$_port and no GITHUB_TOKEN — push will fail" >&2
    fi
  fi
fi
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

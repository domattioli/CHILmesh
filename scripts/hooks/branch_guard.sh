#!/bin/bash
# PreToolUse:Bash hook — enforces branch policy from CLAUDE.md.
# Blocks: creating claude/* branches, commit/push on claude/*, --no-verify,
# force-push to main, refspec colon-rename (issue #31 antipattern).
# Override: CLAUDE_BRANCH_OVERRIDE=1 (logged to ~/.claude/hook-bypass.log).
# Exit 2 with stderr -> blocks tool call, surfaces reason to model.
# Synced from DomI upstream. See domattioli/DomI scripts/hooks/branch_guard.sh.

set -uo pipefail

input="$(cat)"
cmd="$(echo "$input" | jq -r '.tool_input.command // ""' 2>/dev/null)"
[ -z "$cmd" ] && exit 0

echo "$cmd" | grep -qE '(^|[^[:alnum:]_])git([^[:alnum:]_]|$)' || exit 0

branch="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "")"
override="${CLAUDE_BRANCH_OVERRIDE:-}"

log_bypass() {
  local reason="$1"
  mkdir -p "$HOME/.claude" 2>/dev/null
  printf '%s  branch_guard  bypass=%s  branch=%s  reason=%s  cmd=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$override" "$branch" "$reason" "$cmd" \
    >> "$HOME/.claude/hook-bypass.log" 2>/dev/null
}

if echo "$cmd" | grep -qE 'git[[:space:]]+(checkout[[:space:]]+-b|switch[[:space:]]+-c|branch)[[:space:]]+claude/'; then
  if [ "$override" = "1" ]; then
    log_bypass "create-claude-branch"
  else
    echo "BLOCKED (branch_guard): creating claude/* branch violates branch policy. Override: CLAUDE_BRANCH_OVERRIDE=1" >&2
    exit 2
  fi
fi

if echo "$cmd" | grep -qE 'git[[:space:]]+(commit|push)([[:space:]]|$).*--no-verify'; then
  echo "BLOCKED (branch_guard): --no-verify forbidden by policy" >&2
  exit 2
fi

if echo "$cmd" | grep -qE 'git[[:space:]]+push([[:space:]]|$).*(--force\b|--force-with-lease\b|[[:space:]]-f\b)'; then
  if [ "$branch" = "main" ] || [ "$branch" = "master" ]; then
    echo "BLOCKED (branch_guard): force-push to $branch refused" >&2
    exit 2
  fi
fi

if echo "$cmd" | grep -qE 'git[[:space:]]+push[[:space:]]+[^[:space:]]+[[:space:]]+[^[:space:]]+:[^[:space:]]+'; then
  echo "BLOCKED (branch_guard): refspec colon-rename (local:remote) creates new remote branch. See issue #31." >&2
  exit 2
fi

if echo "$cmd" | grep -qE 'git[[:space:]]+(commit|push)([[:space:]]|$)' && \
   echo "$branch" | grep -qE '^claude/'; then
  if [ "$override" = "1" ]; then
    log_bypass "commit-push-on-claude-branch"
  else
    echo "BLOCKED (branch_guard): on $branch (claude/* branch). Switch to working branch. Override: CLAUDE_BRANCH_OVERRIDE=1" >&2
    exit 2
  fi
fi

exit 0

#!/bin/bash
# PreToolUse:Bash hook — validates git commit message format.
# Requires '^(fix|feat|docs|chore|refactor|test): ' prefix.
# Synced from DomI upstream. See domattioli/DomI scripts/hooks/commit_msg_guard.sh.

set -uo pipefail

input="$(cat)"
cmd="$(echo "$input" | jq -r '.tool_input.command // ""' 2>/dev/null)"
[ -z "$cmd" ] && exit 0

echo "$cmd" | grep -qE 'git[[:space:]]+commit\b' || exit 0
echo "$cmd" | grep -qE '(^|[[:space:]])--amend\b' && exit 0
echo "$cmd" | grep -qE '<<.*EOF' && exit 0

msg=""
if [[ "$cmd" =~ -m[[:space:]]+\"([^\"]+)\" ]]; then
  msg="${BASH_REMATCH[1]}"
elif [[ "$cmd" =~ -m[[:space:]]+'([^']+)' ]]; then
  msg="${BASH_REMATCH[1]}"
fi

[ -z "$msg" ] && exit 0

echo "$msg" | grep -qE '^(Merge |Revert ")' && exit 0

if echo "$msg" | grep -qiE '^(wip|fixup!|squash!|tmp|test commit)([: ]|$)'; then
  echo "BLOCKED (commit_msg_guard): forbidden prefix: '$msg'" >&2
  exit 2
fi

if ! echo "$msg" | grep -qE '^(fix|feat|docs|chore|refactor|test): '; then
  echo "BLOCKED (commit_msg_guard): must match '^(fix|feat|docs|chore|refactor|test): '. Got: '$msg'" >&2
  exit 2
fi

exit 0

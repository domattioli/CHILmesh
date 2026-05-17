#!/bin/bash
# Stop hook — emits introspect reminder when session made changes.
# Advisory only: exit 0. Cron-safe: gated on CLAUDE_INTERACTIVE=1.
# Synced from DomI upstream. See domattioli/DomI scripts/hooks/stop_introspect.sh.

set -uo pipefail

input="$(cat)"
stop_active="$(echo "$input" | jq -r '.stop_hook_active // false' 2>/dev/null)"
[ "$stop_active" = "true" ] && exit 0

[ "${CLAUDE_INTERACTIVE:-}" != "1" ] && exit 0

git rev-parse --git-dir >/dev/null 2>&1 || exit 0

made_changes=0
[ -n "$(git status --porcelain 2>/dev/null)" ] && made_changes=1
ahead="$(git rev-list --count '@{upstream}..HEAD' 2>/dev/null || echo 0)"
[ "$ahead" -gt 0 ] 2>/dev/null && made_changes=1

[ "$made_changes" = "0" ] && exit 0

echo "Session made changes. Run /introspect to capture lessons before exit." >&2
exit 0

#!/bin/bash
# SessionStart hook — runs instructions_on_start.sh if present.
# Synced from DomI upstream. See domattioli/DomI scripts/hooks/session_start.sh.

set -uo pipefail

[ "${CLAUDE_HOOK_SESSION_STARTED:-}" = "1" ] && exit 0
export CLAUDE_HOOK_SESSION_STARTED=1

repo_root="$(git rev-parse --show-toplevel 2>/dev/null || echo "")"
[ -z "$repo_root" ] && exit 0
[ -f "$repo_root/scripts/instructions_on_start.sh" ] || exit 0

bash "$repo_root/scripts/instructions_on_start.sh" 1>&2
exit $?

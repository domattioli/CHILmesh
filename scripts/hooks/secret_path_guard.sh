#!/bin/bash
# PreToolUse:Write|Edit hook — blocks writes to secret-bearing paths.
# Synced from DomI upstream. See domattioli/DomI scripts/hooks/secret_path_guard.sh.

set -uo pipefail

input="$(cat)"
path="$(echo "$input" | jq -r '.tool_input.file_path // ""' 2>/dev/null)"
[ -z "$path" ] && exit 0

case "$path" in
  *.env.example|*.env.template|*.env.test|*.env.sample) exit 0 ;;
  *token-rotation*|*REFRESH_TOKEN_DOCS*|*token_rotation*) exit 0 ;;
  *secret-example*|*credentials.example*) exit 0 ;;
esac

if echo "$path" | grep -qE '(^|/)\..env$|\.pem$|/credentials\.|/secrets?(\.|$)|/token(\.|$)|\.gpg$|id_rsa$|id_ed25519$'; then
  echo "BLOCKED (secret_path_guard): $path matches secret-bearing pattern. Rename to *.example or *.template variant." >&2
  exit 2
fi

exit 0

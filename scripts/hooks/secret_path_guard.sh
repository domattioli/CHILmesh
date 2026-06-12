#!/bin/bash
# PreToolUse:Write|Edit hook — blocks writes to secret-bearing paths.
# Pattern synced from DomI canonical copy 2026-06-12 (fixed .env + root-relative bypass).

set -uo pipefail

input="$(cat)"
path="$(echo "$input" | jq -r '.tool_input.file_path // ""' 2>/dev/null)"
[ -z "$path" ] && exit 0

# Exemptions — return 0 fast
case "$path" in
  *.env.example|*.env.template|*.env.test|*.env.sample) exit 0 ;;
  *token-rotation*|*REFRESH_TOKEN_DOCS*|*token_rotation*) exit 0 ;;
  *secret-example*|*credentials.example*) exit 0 ;;
esac

# Block paths matching secret patterns.
# `.env` matches both top-level `.env` and nested `*/foo.env` per CLAUDE.md
# hard stop ("Never commit *.env"). Exemptions above (.env.example etc.) run
# first.
if echo "$path" | grep -qE '\.env$|\.pem$|(^|/)credentials\.|(^|/)secrets?(\.|$)|(^|/)token(\.|$)|\.gpg$|id_rsa$|id_ed25519$'; then
  echo "BLOCKED (secret_path_guard): $path matches secret-bearing pattern. If intentional, rename to *.example or *.template variant." >&2
  exit 2
fi

exit 0

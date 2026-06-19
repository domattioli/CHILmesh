#!/bin/bash
# PreToolUse:Read|Bash hook — blocks raw reads of large mesh files.
# Reroutes to lazy 'chilmesh summary <file>' CLI to avoid pulling multi-MB
# node/element arrays into context (issue #201 context-management pattern).
# Override: CHILMESH_MESH_READ_GUARD_BYPASS=1 (logged to ~/.claude/hook-bypass.log).
# Exit 2 with stderr -> blocks tool call, surfaces reason to model.

set -uo pipefail

input="$(cat)"
tool="$(echo "$input" | jq -r '.tool_name // ""' 2>/dev/null)"
[ -z "$tool" ] && exit 0

# Config
max_kb="${CHILMESH_MESH_READ_MAX_KB:-64}"
bypass="${CHILMESH_MESH_READ_GUARD_BYPASS:-}"

log_bypass() {
  local reason="$1"
  mkdir -p "$HOME/.claude" 2>/dev/null
  printf '%s  mesh_read_guard  bypass=%s  reason=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$bypass" "$reason" \
    >> "$HOME/.claude/hook-bypass.log" 2>/dev/null
}

# Helper: is path a mesh file?
is_mesh_file() {
  local p="$1"
  [ -z "$p" ] && return 1
  local base
  base="$(basename "$p")"
  # Match .14 .grd .2dm .13 .msh .npy .npz or fort.NNN
  echo "$base" | grep -qiE '\.(14|grd|2dm|13|msh|npy|npz)$|^fort\.[0-9]+$'
}

# Helper: is file size over threshold?
over_threshold() {
  local p="$1"
  [ ! -f "$p" ] && return 1
  local size_kb
  size_kb=$(($(stat -c%s "$p" 2>/dev/null || echo 0) / 1024))
  [ "$size_kb" -gt "$max_kb" ]
}

# Branch on tool type
case "$tool" in
  Read)
    path="$(echo "$input" | jq -r '.tool_input.file_path // ""' 2>/dev/null)"
    if [ -z "$path" ]; then
      exit 0
    fi
    if is_mesh_file "$path" && over_threshold "$path"; then
      if [ "$bypass" = "1" ]; then
        log_bypass "bypass=1"
        exit 0
      fi
      size_kb=$(($(stat -c%s "$path" 2>/dev/null || echo 0) / 1024))
      echo "BLOCKED (mesh_read_guard): $path (${size_kb}KB) is mesh data — reading raw mesh arrays wastes context. Use 'chilmesh summary <path>' (header-only, lazy) instead, or set CHILMESH_MESH_READ_GUARD_BYPASS=1 to override." >&2
      exit 2
    fi
    exit 0
    ;;
  Bash)
    cmd="$(echo "$input" | jq -r '.tool_input.command // ""' 2>/dev/null)"
    [ -z "$cmd" ] && exit 0

    # Allow chilmesh summary (reroute target)
    if echo "$cmd" | grep -q 'chilmesh summary'; then
      exit 0
    fi

    # Extract first token (leading command)
    first_token="$(echo "$cmd" | sed -E 's/^[[:space:]]*//; s/[[:space:]].*//')"

    # Check if it's a read command
    case "$first_token" in
      cat|head|tail|less|more|od|xxd|strings|nl|bat|hexdump)
        # Read command; tokenize and check each token for mesh file
        tokens="$(echo "$cmd" | tr ' ' '\n' | grep -v '^$')"
        while IFS= read -r token; do
          if is_mesh_file "$token" && over_threshold "$token"; then
            if [ "$bypass" = "1" ]; then
              log_bypass "bypass=1"
              exit 0
            fi
            size_kb=$(($(stat -c%s "$token" 2>/dev/null || echo 0) / 1024))
            echo "BLOCKED (mesh_read_guard): $token (${size_kb}KB) is mesh data — reading raw mesh arrays wastes context. Use 'chilmesh summary <path>' (header-only, lazy) instead, or set CHILMESH_MESH_READ_GUARD_BYPASS=1 to override." >&2
            exit 2
          fi
        done <<< "$tokens"
        exit 0
        ;;
      *)
        # Not a read command (e.g., echo, cat for write, etc.)
        exit 0
        ;;
    esac
    ;;
  *)
    exit 0
    ;;
esac

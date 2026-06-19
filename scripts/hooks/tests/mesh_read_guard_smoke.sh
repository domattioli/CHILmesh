#!/bin/bash
# Smoke test for mesh_read_guard.sh PreToolUse hook.
# Creates temp mesh files and crafted JSON, pipes through hook, asserts exit codes.

set -uo pipefail

fail=0
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# Create test files
head -c 71680 /dev/zero > "$tmp/big.14"
head -c 1024 /dev/zero > "$tmp/small.14"
head -c 71680 /dev/zero > "$tmp/code.py"
head -c 71680 /dev/zero > "$tmp/big.npy"

run() {
  local json="$1"
  echo "$json" | bash /home/user/CHILmesh/scripts/hooks/mesh_read_guard.sh 2>/dev/null
  echo "$?"
}

expect() {
  local desc="$1"
  local got="$2"
  local want="$3"
  if [ "$got" = "$want" ]; then
    echo "PASS: $desc"
  else
    echo "FAIL: $desc (got $got want $want)"
    fail=1
  fi
}

# Scenario 1: Read big.14 => exit 2
result=$(run "{\"tool_name\": \"Read\", \"tool_input\": {\"file_path\": \"$tmp/big.14\"}}")
expect "Read big.14" "$result" "2"

# Scenario 2: Read small.14 => exit 0
result=$(run "{\"tool_name\": \"Read\", \"tool_input\": {\"file_path\": \"$tmp/small.14\"}}")
expect "Read small.14" "$result" "0"

# Scenario 3: Read code.py => exit 0 (wrong extension)
result=$(run "{\"tool_name\": \"Read\", \"tool_input\": {\"file_path\": \"$tmp/code.py\"}}")
expect "Read code.py" "$result" "0"

# Scenario 4: Read big.npy => exit 2
result=$(run "{\"tool_name\": \"Read\", \"tool_input\": {\"file_path\": \"$tmp/big.npy\"}}")
expect "Read big.npy" "$result" "2"

# Scenario 5: Bash cat big.14 => exit 2
result=$(run "{\"tool_name\": \"Bash\", \"tool_input\": {\"command\": \"cat $tmp/big.14\"}}")
expect "Bash cat big.14" "$result" "2"

# Scenario 6: Bash head big.14 => exit 2
result=$(run "{\"tool_name\": \"Bash\", \"tool_input\": {\"command\": \"head $tmp/big.14\"}}")
expect "Bash head big.14" "$result" "2"

# Scenario 7: Bash chilmesh summary big.14 => exit 0 (reroute target)
result=$(run "{\"tool_name\": \"Bash\", \"tool_input\": {\"command\": \"chilmesh summary $tmp/big.14\"}}")
expect "Bash chilmesh summary big.14" "$result" "0"

# Scenario 8: Bash cat code.py => exit 0 (not a mesh file)
result=$(run "{\"tool_name\": \"Bash\", \"tool_input\": {\"command\": \"cat $tmp/code.py\"}}")
expect "Bash cat code.py" "$result" "0"

# Scenario 9: Bash echo (write) big.14 => exit 0 (write, not read)
result=$(run "{\"tool_name\": \"Bash\", \"tool_input\": {\"command\": \"echo hi > $tmp/big.14\"}}")
expect "Bash echo write big.14" "$result" "0"

# Scenario 10: Read big.14 with CHILMESH_MESH_READ_GUARD_BYPASS=1 => exit 0
result=$(CHILMESH_MESH_READ_GUARD_BYPASS=1 run "{\"tool_name\": \"Read\", \"tool_input\": {\"file_path\": \"$tmp/big.14\"}}")
expect "Read big.14 with bypass=1" "$result" "0"

# Scenario 11: Bash cat small.14 => exit 0 (under threshold)
result=$(run "{\"tool_name\": \"Bash\", \"tool_input\": {\"command\": \"cat $tmp/small.14\"}}")
expect "Bash cat small.14" "$result" "0"

if [ "${fail:-0}" = 1 ]; then
  echo ""
  echo "SMOKE FAILED"
  exit 1
else
  echo ""
  echo "ALL SMOKE PASS"
  exit 0
fi

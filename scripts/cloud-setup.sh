#!/bin/bash
# Cloud-environment setup script for CHILmesh.
#
# Pinned to the planning-optimize_modernize branch — does not create new
# branches and does not run spec-kit candidate scanning. Robust against
# missing specs/ directories under `set -euo pipefail`.
#
# Paste this into your cloud environment's "setup script" field.

set -euo pipefail

REPO_URL="${1:-https://github.com/domattioli/CHILmesh.git}"
REPO_DIR="${2:-/workspace/$(basename "$REPO_URL" .git)}"
BASE_BRANCH="main"
TARGET_BRANCH="planning-optimize_modernize"

echo "===> Routine setup starting"
echo "===> Repository: $REPO_URL"
echo "===> Target directory: $REPO_DIR"
echo "===> Pinned branch: $TARGET_BRANCH"

# --------------------------------------------------------------------
# 1. Prerequisites
# --------------------------------------------------------------------
if ! command -v node >/dev/null 2>&1; then
  echo "===> Installing Node.js 20.x"
  curl -fsSL https://deb.nodesource.com/setup_20.x | bash -
  apt-get update
  apt-get install -y nodejs
fi

# --------------------------------------------------------------------
# 2. Fresh clone (full history so we can check out feature branches)
# --------------------------------------------------------------------
echo "===> Cloning repository"
rm -rf "$REPO_DIR"
mkdir -p "$(dirname "$REPO_DIR")"
git clone "$REPO_URL" "$REPO_DIR"
cd "$REPO_DIR"

# Append the stream-timeout block to CLAUDE.md (creates the file if absent).
if ! grep -q "Stream Timeout Prevention" CLAUDE.md 2>/dev/null; then
  echo "===> Appending stream-timeout block to CLAUDE.md"
  cat >> CLAUDE.md <<'CLAUDEMD'

## Stream Timeout Prevention

1. Do each numbered task ONE AT A TIME. Complete one task fully,
   confirm it worked, then move to the next.
2. Never write a file longer than ~150 lines in a single tool call.
   If a file will be longer, write it in multiple append/edit passes.
3. Start a fresh session if the conversation gets long (20+ tool calls).
   The error gets worse as the session grows.
4. Keep individual grep/search outputs short. Use flags like
   --include and -l (list files only) to limit output size.
5. If you do hit the timeout, retry the same step in a shorter form.
   Don't repeat the entire task from scratch.
CLAUDEMD
fi

# --------------------------------------------------------------------
# 3. Checkout the pinned planning branch — NEVER create new branches.
# --------------------------------------------------------------------
echo "===> Checking out pinned branch: $TARGET_BRANCH"
git fetch --quiet origin "$TARGET_BRANCH"
if git show-ref --verify --quiet "refs/heads/$TARGET_BRANCH"; then
  git checkout "$TARGET_BRANCH"
  git pull --ff-only origin "$TARGET_BRANCH" || \
    echo "===> WARN: fast-forward pull failed; staying on local copy"
else
  git checkout -B "$TARGET_BRANCH" "origin/$TARGET_BRANCH"
fi

# --------------------------------------------------------------------
# 4. Python deps
# --------------------------------------------------------------------
PIP_ARGS="--retries 5 --timeout 120 --prefer-binary"
pip_install() {
  local attempt=1
  local max_attempts=3
  while [ $attempt -le $max_attempts ]; do
    if python3 -m pip install $PIP_ARGS "$@"; then
      return 0
    fi
    echo "===> pip install attempt $attempt failed; sleeping $((attempt * 10))s"
    sleep $((attempt * 10))
    attempt=$((attempt + 1))
  done
  echo "===> pip install gave up after $max_attempts attempts: $*" >&2
  return 1
}

{
  pip_install --upgrade pip || true
  if [ -f pyproject.toml ]; then
    echo "===> Installing CHILmesh[dev,publish]"
    pip_install -e ".[dev,publish]" \
      || pip_install -e ".[dev]" \
      || pip_install -e . \
      || echo "===> WARN: pyproject install failed; continuing"
  fi
} || echo "===> WARN: python deps step had failures; session can install lazily"

# --------------------------------------------------------------------
# 5. GSD
# --------------------------------------------------------------------
gsd_install() {
  local attempt=1
  while [ $attempt -le 3 ]; do
    if npx --yes get-shit-done-cc@latest --claude --global; then
      return 0
    fi
    echo "===> GSD install attempt $attempt failed; sleeping $((attempt * 10))s"
    sleep $((attempt * 10))
    attempt=$((attempt + 1))
  done
  echo "===> GSD install failed after 3 attempts" >&2
  return 1
}
gsd_install || echo "===> WARN: GSD not installed; /gsd-* commands will be missing"

# --------------------------------------------------------------------
# 6. Superpowers plugin
# --------------------------------------------------------------------
mkdir -p "$HOME/.claude"
if command -v claude >/dev/null 2>&1 && claude plugin --help >/dev/null 2>&1; then
  echo "===> Installing superpowers plugin"
  claude plugin marketplace add obra/superpowers-marketplace || true
  claude plugin install superpowers@obra/superpowers-marketplace || true
else
  echo "===> Pre-seeding plugins.json"
  cat > "$HOME/.claude/plugins.json" <<'JSON'
{
  "marketplaces": {
    "obra/superpowers-marketplace": {
      "source": "github:obra/superpowers-marketplace"
    }
  },
  "plugins": {
    "superpowers": {
      "marketplace": "obra/superpowers-marketplace",
      "enabled": true
    }
  }
}
JSON
fi

# --------------------------------------------------------------------
# 7. Sanity output
# --------------------------------------------------------------------
cd "$REPO_DIR"
BRANCH=$(git rev-parse --abbrev-ref HEAD)
COMMIT=$(git rev-parse --short HEAD)
GSD_COUNT=$(ls "$HOME/.claude/skills/" 2>/dev/null | grep -c '^gsd-' || true)
ACTIVE_SPEC="$( { ls -d specs/[0-9][0-9][0-9]-* 2>/dev/null || true; } | tail -1)"
[ -z "$ACTIVE_SPEC" ] && ACTIVE_SPEC="none"

echo ""
echo "===> Setup complete!"
echo "===> Repository: $(basename "$REPO_URL" .git)"
echo "===> Branch: $BRANCH @ $COMMIT"
echo "===> Directory: $REPO_DIR"
echo "===> GSD skills available: $GSD_COUNT"

if command -v CHILmesh >/dev/null 2>&1; then
  echo "===> Validating bundled manifest..."
  CHILmesh validate || echo "===> WARN: bundled manifest validation failed"
else
  echo "===> WARN: CHILmesh CLI not on PATH; in-session install needed"
fi

echo "===> Active spec folder: $ACTIVE_SPEC"
echo "===> Guidelines: $REPO_DIR/CLAUDE.md and $REPO_DIR/.specify/memory/constitution.md"

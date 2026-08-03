#!/usr/bin/env bash
# Measure the pixi lockfile forward/backward compatibility matrix.
#
# Motivating question: if we adopt `requires-pixi` as a LOWER BOUND, a contributor
# bisecting through old commits keeps whatever (newer) pixi they have installed. That
# is only workable if a newer pixi can consume an older lockfile without rewriting it.
# Conversely, if someone on an old pixi pulls a branch whose lock was written by a
# newer pixi, we need to know exactly what they see.
#
# pixi 0.68.0 bumped the lock format v6 -> v7, so 0.67.2 (this repo's pin) and 0.75.0
# (current latest) straddle the boundary and make a clean natural experiment.
#
# Run:  bash experiments/claude/pixi_lockfile_compat_matrix.sh 2>&1 \
#         | tee experiments/claude/logs/pixi_lockfile_compat_matrix.log
set -uo pipefail

OLD_PIXI="$(command -v pixi)"
WORK="${TMPDIR:-/tmp}/pixi_compat_$$"
NEW_VER="0.75.0"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

echo "### old pixi: $OLD_PIXI ($($OLD_PIXI --version))"

# --- fetch the new pixi into an isolated dir (never touches ~/.pixi) -----------------
echo
echo "### downloading pixi $NEW_VER"
curl -fsSL -o "$WORK/pixi-new.tar.gz" \
  "https://github.com/prefix-dev/pixi/releases/download/v${NEW_VER}/pixi-x86_64-unknown-linux-musl.tar.gz"
mkdir -p "$WORK/newbin"
tar -xzf "$WORK/pixi-new.tar.gz" -C "$WORK/newbin"
NEW_PIXI="$WORK/newbin/pixi"
chmod +x "$NEW_PIXI"
echo "### new pixi: $($NEW_PIXI --version)"

# A deliberately tiny workspace: we are testing lockfile FORMAT handling, not solving.
make_ws() {
  local dir="$1"
  mkdir -p "$dir"
  cat >"$dir/pixi.toml" <<'EOF'
[workspace]
name = "compat"
channels = ["conda-forge"]
platforms = ["linux-64"]

[dependencies]
xz = "*"
EOF
}

lock_version() { grep -m1 '^version:' "$1/pixi.lock" 2>/dev/null || echo "version: <none>"; }

# ------------------------------------------------------------------ case 1: v6 lock
echo
echo "======================================================================"
echo "CASE 1: lock written by OLD pixi (expect v6), then read by NEW pixi"
echo "======================================================================"
make_ws "$WORK/c1"
(cd "$WORK/c1" && "$OLD_PIXI" lock >/dev/null 2>&1)
echo "--- lock as written by old pixi: $(lock_version "$WORK/c1")"
cp "$WORK/c1/pixi.lock" "$WORK/c1/pixi.lock.orig"

echo "--- NEW pixi: does 'lock --check' accept the v6 lock as up to date?"
(cd "$WORK/c1" && "$NEW_PIXI" lock --check 2>&1 | tail -5)
echo "    exit=$?"

echo "--- NEW pixi: did merely reading it rewrite the file?"
if cmp -s "$WORK/c1/pixi.lock" "$WORK/c1/pixi.lock.orig"; then
  echo "    UNCHANGED (good: bisecting would not dirty the working tree)"
else
  echo "    REWRITTEN -> now $(lock_version "$WORK/c1")"
  echo "    (this would show up as a modified pixi.lock during a bisect)"
fi

echo "--- NEW pixi: does an install --locked against the v6 lock succeed?"
(cd "$WORK/c1" && "$NEW_PIXI" install --locked 2>&1 | tail -8)
echo "    exit=$?"
echo "--- lock after install --locked: $(lock_version "$WORK/c1")"
if cmp -s "$WORK/c1/pixi.lock" "$WORK/c1/pixi.lock.orig"; then
  echo "    still UNCHANGED"
else
  echo "    CHANGED by install --locked"
fi

# ------------------------------------------------------------------ case 2: v7 lock
echo
echo "======================================================================"
echo "CASE 2: lock written by NEW pixi (expect v7), then read by OLD pixi"
echo "======================================================================"
make_ws "$WORK/c2"
(cd "$WORK/c2" && "$NEW_PIXI" lock >/dev/null 2>&1)
echo "--- lock as written by new pixi: $(lock_version "$WORK/c2")"

echo "--- OLD pixi reading a v7 lock (this is the error a stale contributor sees):"
(cd "$WORK/c2" && "$OLD_PIXI" install --locked 2>&1 | tail -15)
echo "    exit=$?"

# -------------------------------------------------- case 3: requires-pixi guard text
echo
echo "======================================================================"
echo "CASE 3: requires-pixi lower bound -- exact wording of the guard"
echo "======================================================================"
make_ws "$WORK/c3"
sed -i 's/^platforms = .*/&\nrequires-pixi = ">=0.75.0"/' "$WORK/c3/pixi.toml"
echo "--- OLD pixi (0.67.2) against requires-pixi = \">=0.75.0\":"
(cd "$WORK/c3" && "$OLD_PIXI" run echo hi 2>&1 | tail -15)
echo "    exit=$?"

echo
echo "--- does the guard fire on a plain 'pixi install' too, or only 'run'?"
(cd "$WORK/c3" && "$OLD_PIXI" install 2>&1 | tail -6)
echo "    exit=$?"

echo
echo "### done"

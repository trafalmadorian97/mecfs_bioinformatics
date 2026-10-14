# Renovate pixi workflow-version automation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Have Renovate keep the pixi `pixi-version` pins in `.github/workflows/*.yml` current, atomically with the `constraints.pixi` pin in `renovate.json`, and let pixi bumps automerge.

**Architecture:** Add one regex customManager to `renovate.json` that matches `pixi-version: vX.Y.Z` in every workflow file. Because it emits the same `depName`/`datasource`/`currentValue` as the existing `constraints.pixi` customManager, Renovate groups both edits onto one branch (see spec section 2). Flip the existing `prefix-dev/pixi` packageRule to `automerge: true` and rewrite its now-obsolete description.

**Tech Stack:** Renovate config (JSON), regex customManager, GitHub releases datasource.

**Spec:** `experiments/claude/design_specs/2026-08-02-renovate-pixi-workflow-version.md`

## Global Constraints

- `constraints.pixi` in `renovate.json` MUST stay an EXACT version (`0.67.2` today), never a range. Renovate resolves the lock-generation pixi as `config.constraints?.pixi ?? requires-pixi`; a range would resolve to newest and could emit a lock CI cannot read.
- Do NOT add `requires-pixi` to `pyproject.toml`. Do NOT add any CI reconciliation check. (Explicitly out of scope per spec.)
- The `v` prefix in `pixi-version: vX.Y.Z` MUST stay OUTSIDE the regex capture group, so it is preserved on rewrite and the captured value is byte-identical to the bare `constraints.pixi` value.
- The new customManager must use `datasourceTemplate: github-releases`, `depNameTemplate: prefix-dev/pixi`, `extractVersionTemplate: ^v(?<version>.+)$` — identical to the existing `constraints.pixi` customManager — so both updates share one branch.

---

### Task 1: Add workflow-pixi customManager, enable pixi automerge

**Files:**
- Modify: `renovate.json` (add one `customManagers` entry; edit the `prefix-dev/pixi` `packageRules` entry)

**Interfaces:**
- Consumes: the existing `constraints.pixi` customManager (`renovate.json` lines ~20-30) as the pattern to mirror for datasource/depName/extractVersion.
- Produces: nothing other tasks depend on (terminal task).

- [ ] **Step 1: Verify the regex captures all three pins and preserves `v` on rewrite**

The repo already has a probe script from the design phase. Run it to confirm the capture group and rewrite behavior before editing config:

```bash
python3 - <<'EOF'
import re, pathlib, glob
# Mirror Renovate's matchString; Renovate replaces ONLY the capture group.
pat = re.compile(r"pixi-version:\s*v(?P<currentValue>\d+\.\d+\.\d+)")
files = sorted(glob.glob(".github/workflows/*.yml"))
hits = 0
for f in files:
    txt = pathlib.Path(f).read_text()
    for m in pat.finditer(txt):
        hits += 1
        new = txt[:m.start("currentValue")] + "0.99.0" + txt[m.end("currentValue"):]
        line = new[new.rfind("\n",0,m.start())+1 : new.find("\n", m.start())]
        assert "v0.99.0" in line, f"v prefix lost in {f}: {line!r}"
        print(f"{f}: captured={m.group('currentValue')!r} -> {line.strip()}")
assert hits == 3, f"expected 3 pins, found {hits}"
print("OK: 3 pins, v preserved on rewrite")
EOF
```

Expected: three lines, each rewritten to `pixi-version: v0.99.0`, then `OK: 3 pins, v preserved on rewrite`. If it does not print OK, STOP — the regex or the workflow files changed and the plan needs revisiting.

- [ ] **Step 2: Add the customManager to `renovate.json`**

In the `customManagers` array, add this entry after the existing `renovate.json`-constraints manager (order does not matter to Renovate; keep it adjacent for readability):

```json
{
  "customType": "regex",
  "managerFilePatterns": ["/^\\.github/workflows/.+\\.yml$/"],
  "matchStrings": ["pixi-version:\\s*v(?<currentValue>\\d+\\.\\d+\\.\\d+)"],
  "datasourceTemplate": "github-releases",
  "depNameTemplate": "prefix-dev/pixi",
  "extractVersionTemplate": "^v(?<version>.+)$"
}
```

- [ ] **Step 3: Flip `prefix-dev/pixi` to automerge and rewrite its description**

Find the `packageRules` entry for `prefix-dev/pixi` (currently `"automerge": false` with a description about coordinating workflow pins by hand). Replace the whole entry with:

```json
{
  "description": "pixi bumps move the three workflow pixi-version pins and constraints.pixi together on one Renovate branch (same depName+datasource+currentValue). Automerge is on: every step is lock-format safe (newer pixi reads older locks), and a contributor on stale pixi gets pixi's native self-update message. See experiments/claude/design_specs/2026-08-02-renovate-pixi-workflow-version.md.",
  "matchPackageNames": ["prefix-dev/pixi"],
  "automerge": true
}
```

- [ ] **Step 4: Validate the JSON**

```bash
python3 -c "import json; json.load(open('renovate.json')); print('renovate.json valid JSON')"
```

Expected: `renovate.json valid JSON`. (Renovate config is JSON5-tolerant, but keeping it strict-JSON-valid avoids surprises; this file is plain JSON today.)

- [ ] **Step 5: Confirm exactly one `automerge: true` change and the exact-pin constraint survived**

```bash
grep -n '"constraints"' -A2 renovate.json          # must still show "pixi": "0.67.2" (exact, no range)
grep -n 'prefix-dev/pixi' renovate.json            # depName present in customManager + packageRule
grep -c '"github-actions/.\+\.yml"' renovate.json || true
```

Expected: `constraints.pixi` is still the exact string `0.67.2`; `prefix-dev/pixi` appears in both a customManager (`depNameTemplate`) and the packageRule (`matchPackageNames`). If `constraints.pixi` is anything but an exact version, STOP — the Global Constraints are violated.

- [ ] **Step 6: Commit**

```bash
git add renovate.json experiments/claude/design_specs/2026-08-02-renovate-pixi-workflow-version.md experiments/claude/design_specs/2026-08-02-renovate-pixi-workflow-version-plan.md experiments/claude/pixi_lockfile_compat_matrix.sh experiments/claude/logs/pixi_lockfile_compat_matrix.log
git commit -m "Renovate: keep workflow pixi-version pins current, automerge pixi

Add a regex customManager that updates pixi-version in all workflow files;
it shares a branch with the constraints.pixi update (same depName/datasource/
currentValue), so CI pins and the lock-generation constraint move atomically.
Flip prefix-dev/pixi to automerge:true.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01CxyJX62EdLsBRjEuF5r4pR"
```

---

## Rollout verification (post-merge, not a code step)

These confirm the design's two live claims. They cannot run locally — they need the self-hosted Renovate workflow — so they happen after the PR merges to `main`.

1. **One-branch grouping (the ~96% claim from spec §2).** Dispatch the Renovate workflow with `dryRun: true`, `logLevel: debug`. In the log, confirm Renovate plans a SINGLE branch touching both `renovate.json` and the workflow files for `prefix-dev/pixi` — not two `prefix-dev/pixi` branches. Two branches means the grouping assumption broke; investigate before enabling further.
2. **Automerge path.** The first real pixi bump PR should update all four pins on one branch, pass CI (new pixi vs. existing older-format lock), and automerge. Watch the first one end-to-end.

## Self-Review notes

- **Spec coverage:** §Design 1 (customManager) → Task 1 Step 2. §Design 2 (one-branch, why safe) → verified Step 1 + Rollout 1. §Design 3 (exact constraints.pixi) → Global Constraints + Step 5 guard. §Design 4 (automerge:true + rewritten description) → Step 3. §"What we are NOT doing" → Global Constraints. All covered.
- **No placeholders:** every step has runnable commands / literal JSON.
- **Consistency:** the customManager JSON matches the spec verbatim; `depName`/`datasource`/`extractVersion` identical to the existing constraints manager, which is the property that guarantees one branch.

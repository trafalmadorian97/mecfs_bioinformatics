# Renovate: automate the pixi version used by CI workflows

Date: 2026-08-02
Status: approved, ready for implementation

## Problem

The pixi version is pinned in four places, all edited by hand today:

- `.github/workflows/on_pr.yml` -> `pixi-version: v0.67.2`
- `.github/workflows/deploy_docs.yml` -> `pixi-version: v0.67.2`
- `.github/workflows/check_links.yml` -> `pixi-version: v0.67.2`
- `renovate.json` -> `constraints.pixi: "0.67.2"` (the version Renovate uses to run `pixi lock`)

Renovate already updates `constraints.pixi` via a customManager, but never the three
workflow pins. So the CI pins drift and are bumped manually. This spec adds Renovate
coverage for the workflow pins so all four move together, and settles the question of
whether to also enforce a contributor-local pixi floor.

## Decision: do NOT adopt `requires-pixi`

We considered adding `requires-pixi` to `pyproject.toml` as a contributor-facing version
floor. We are deliberately NOT doing this. The reasoning rests on measured facts, not
guesses (experiment: `experiments/claude/pixi_lockfile_compat_matrix.sh`, log in
`experiments/claude/logs/pixi_lockfile_compat_matrix.log`):

1. **Lockfile compatibility is asymmetric.** A newer pixi reads an older lockfile; an
   older pixi cannot read a newer one. Measured:

   | pixi | reads v6 lock | reads v7 lock |
   |------|---------------|---------------|
   | 0.67.2 (current pin) | yes | no -- hard fail |
   | 0.75.0 (latest)      | yes (silently rewrites to v7) | yes |

2. **The only thing that has ever broken a contributor is a lock-format bump, and it is
   rare.** Across the last 100 pixi releases (`v0.23.0` May 2024 -> `v0.75.0` Jul 2026),
   spanning 54 distinct minor series, the lock format bumped exactly once (v6 -> v7 at
   0.68.0). pixi has no major versions to hook automation onto -- it has been 0.x for its
   entire life -- and "move the floor every minor" would force ~54 contributor upgrades to
   cover one real incompatibility.

3. **pixi's native too-new-lock error is already better than what `requires-pixi` would
   produce.** The native error names the fix command; the `requires-pixi` guard does not:

   ```
   x Lock-file version 7 is newer than supported
   help: Maximum supported version: 6 (pixi v0.67.2)
         Try running `pixi self-update` to update to the latest version.
   ```

   vs. the `requires-pixi` guard:

   ```
   x this project requires pixi '>=0.75.0', but you have pixi 0.67.2
   help: update pixi to a version that satisfies '>=0.75.0'
   ```

So adding `requires-pixi` would add machinery, a second place to keep in sync, and a
worse error message, to guard a failure mode that occurs about once every two years and
that pixi already reports well. We omit it.

## Design

### 1. New customManager for the workflow pins

Add to `renovate.json` `customManagers`:

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

- A glob (`.+\.yml`) rather than the three literal filenames, so a future workflow that
  uses `setup-pixi` is covered without anyone remembering to update this.
- The `v` prefix is left OUTSIDE the capture group, so it is preserved on rewrite
  (`v0.67.2` -> `v0.75.0`). Verified by simulating the capture-group replacement against
  all three workflow files.

### 2. Atomic branch property (why this is safe)

The captured `currentValue` from the workflows (`0.67.2`) is byte-identical to the value
the existing `constraints.pixi` customManager captures from `renovate.json` (`0.67.2`).
Both customManagers therefore emit the same `depName` (`prefix-dev/pixi`) and the same
`newValue`, so Renovate groups them onto ONE branch. The CI pins and the
lock-generation constraint move atomically -- they can never disagree.

**Why one branch, from the Renovate source (confidence ~96%):**

- Renovate creates one branch per unique `branchName`, and the default template is
  `{{branchPrefix}}{{additionalBranchPrefix}}{{branchTopic}}`, where `branchTopic` is
  `{{depNameSanitized}}-{{newMajor}}...{{newMinor}}.x...` (from
  `lib/config/options/index.ts`). The MANAGER NAME appears nowhere in it, and
  `additionalBranchPrefix` defaults to `''`. So the branch key is a pure function of
  `depName` plus the resolved version fields.
- `lib/workers/repository/updates/branchify.ts` buckets updates into a dict keyed by
  `branchName` (`branchUpgrades[update.branchName] = ...`), and the dedup filters by
  `packageFile`, `depName`, and `currentValue` -- NOT by manager. Updates from different
  managers that share a `branchName` land in the same branch.
- Both managers emit identical `depName` + `datasource` (github-releases) +
  `currentValue`, so they resolve to identical `newMajor`/`newMinor`/`isPatch`, hence an
  identical `branchName`. Every input to the key is the same.
- Strongest evidence: the single workflow customManager already matches `pixi-version` in
  all THREE workflow files, and those three edits land on one branch. The cross-manager
  case runs the same code path (the manager-free `branchName` key), so it groups for the
  same reason.

**Residual risk (~4%), all observable or under our control:** a future `packageRules`
entry setting `additionalBranchPrefix`/`separateMinorPatch`/`groupName` differently (we
have none; both managers match the same `prefix-dev/pixi` rule), or the two
`depNameTemplate`s drifting apart (byte-identical today). The definitive, zero-cost check
is built into rollout: the first dispatched dry-run lists the branches Renovate WOULD
create. Two `prefix-dev/pixi` branches instead of one would reveal a split before any real
bump. Confirm this on the first dry-run.

### 3. Keep `constraints.pixi` exact

`constraints.pixi` in `renovate.json` stays an exact pin, and its customManager stays.
This is load-bearing: Renovate resolves the pixi tool version as
`config.constraints?.pixi ?? requires-pixi` (see the pixi manager's `artifacts.ts`). An
exact pin forces Renovate to run `pixi lock` with exactly the intended version. A range
here would resolve to the newest satisfying pixi and could silently emit a v7 lock that
the pinned CI cannot read.

### 4. automerge: true for pixi

Change the existing `prefix-dev/pixi` packageRule from `automerge: false` to
`automerge: true`, and rewrite its now-obsolete description (which claimed the workflow
pins must be coordinated by hand -- this spec automates that).

Tradeoff, recorded deliberately: with automerge on, the eventual lock-format transition
happens with no human watching the merge. This is acceptable because every step of the
flow below is safe -- the only party ever pushed to upgrade is a contributor on a stale
pixi, who receives the well-worded native `self-update` message.

## End-to-end flow of a pixi bump

1. Renovate opens ONE PR bumping `constraints.pixi` + all three workflow `pixi-version`
   pins.
2. That PR changes no dependencies, so it triggers no lock regeneration. CI runs the NEW
   pixi against the EXISTING (older-format) lock -- fine, since newer pixi reads older
   locks. PR automerges.
3. The next ordinary dependency PR regenerates the lock under the new constraint. If the
   new pixi bumped the lock format, that PR is the first to carry the new-format lock. CI
   is already new enough to read it. Green.
4. A contributor still on old pixi who pulls a new-format lock hits pixi's native error,
   which names `pixi self-update`. `pixi self-update` is NOT blocked by an unsatisfiable
   workspace state (verified), so there is no catch-22.

## What we are NOT doing

- Not adding `requires-pixi` to `pyproject.toml`.
- Not adding a CI check that reconciles lock-format version against a version floor.
- Not changing how contributors install or upgrade pixi (`pixi self-update` remains the
  manual, well-messaged path).

## Out of scope / follow-up

- CI currently pins pixi (`v0.67.2`) BELOW what a contributor might upgrade to. Because
  newer pixi reads older locks but not vice versa, an exact CI pin means any contributor
  who upgrades and regenerates the lock turns CI red until the pin catches up. This spec
  keeps CI pinned (so `constraints.pixi` can stay exact and control lock generation), and
  relies on Renovate to keep that pin current. Letting CI float to newest is a larger
  change with its own tradeoffs and is not part of this work.

## Verification

- `experiments/claude/pixi_lockfile_compat_matrix.sh` -- reproduces the compatibility
  matrix, the requires-pixi guard wording, and the self-update-not-blocked result.
- Regex rewrite simulated against all three workflow files: capture is `0.67.2`, rewrite
  preserves the `v` prefix.
- Renovate pixi manager source confirmed: `requires-pixi` is parsed into the schema but
  never emitted as an updatable dependency (so part B genuinely needs the customManager),
  and `artifacts.ts` resolves the tool version as
  `config.constraints?.pixi ?? requires-pixi`.
- One-branch grouping confirmed against Renovate source: default `branchName`/`branchTopic`
  in `lib/config/options/index.ts` contain no manager component; `branchify.ts` buckets by
  `branchName` and dedups by `packageFile`/`depName`/`currentValue`, not by manager. See
  "Why one branch" above. Verify empirically on the first dispatched dry-run.

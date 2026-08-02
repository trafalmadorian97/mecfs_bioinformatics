# Self-Hosted Renovate — Design

**Date:** 2026-08-01
**Status:** Approved, not yet implemented
**Repo:** `trafalmadorian97/mecfs_bioinformatics` (public, single repo)

## Problem

Renovate 43.288.0 gated `pixi lock` behind the `allowedUnsafeExecutions` global option
(renovatebot/renovate#44939), on the grounds that conda package hooks are arbitrary code
execution. `allowedUnsafeExecutions` is settable only by a self-hosted administrator, and we
run the Mend-hosted app.

The result: Renovate updates `pyproject.toml` but leaves `pixi.lock` untouched, so CI's
`pixi install --locked` fails with `lock-file not up-to-date with the workspace`. It fails
silently — a `logger.once.warn` in Renovate's job log, no PR comment, no "Artifact update
problem" section. First observed on PR #983 (ty 0.0.65); fixed by hand in commit `e6ca20c5`.

Every future Renovate PR touching `[tool.pixi.*]` fails the same way, and the Saturday
`lockFileMaintenance` run is now a no-op.

Self-hosting Renovate via GitHub Actions lets us set `allowedUnsafeExecutions: ['pixi']` and
restore the previous behaviour.

## Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Auth | GitHub App | Replaces a GitHub App (Mend's), so PR attribution and `gitAuthor` semantics are unchanged. No expiry cliff. Renovate's documented recommendation for bot accounts. |
| Cadence | Hourly | Fast feedback while the setup is new; matches the hosted app's responsiveness. Actions minutes are free on a public repo. May be reduced later — see the `lockFileMaintenance` note. |
| Cutover | Dry run, then switch | Validates against a live repo without opening PRs. |
| Approvals | Approve step using `GITHUB_TOKEN` | `renovate-approve` does not work self-hosted (see below). No extra App or hosted service. |
| Lock updates | Renovate does them in-process | Avoids a second workflow, `gitIgnoredAuthors`, and the branch-modified problem entirely. |

Rejected: a separate lockfile-refresh workflow with a lower-privilege token. Its apparent
privilege-separation benefit is notional — the refresh job still needs push rights to
branches whose contents CI executes — and it costs two extra moving parts plus a
`gitIgnoredAuthors` entry that widens what automerges.

Rejected: raw `docker run` of the Renovate image. Hand-rolls what
`renovatebot/github-action` already does correctly.

## Architecture

### Components

| Piece | New? | Purpose |
|---|---|---|
| Renovate GitHub App | new | Bot identity. Installed on this repo only. Supplies a 1-hour token per run. |
| `.github/workflows/renovate.yml` | new | Hourly cron + `workflow_dispatch`. Mints token, runs Renovate, approves PRs. |
| `.github/renovate-global.js` | new | Self-hosted-only config (~10 lines). |
| `renovate.json` | edited | Dependency policy. Two small edits. |
| Repo secrets | new | `RENOVATE_APP_ID`, `RENOVATE_APP_PRIVATE_KEY`. |

### Data flow

```
hourly cron  (or manual workflow_dispatch)
  -> actions/create-github-app-token  -> ghs_... token, valid 1h
  -> renovatebot/github-action (pinned)
       reads renovate.json from the repo
       resolves deps; for pixi changes runs `pixi lock` in-process
         (enabled by allowedUnsafeExecutions: ['pixi'])
       pushes branch + opens PR as <app>[bot], gitAuthor matching
  -> approve step: GITHUB_TOKEN approves open PRs authored by <app>[bot]
  -> on_pr checks run (App-token push, so workflow triggers fire normally)
  -> GitHub native auto-merge merges once checks + approval are green
```

Renovate commits `pixi.lock` in the same commit as the `pyproject.toml` change, carrying its
own `gitAuthor`. `isBranchModified()` therefore stays false and automerge is never blocked.

### Config split

Governing rule: `renovate-global.js` holds **only** what cannot live in `renovate.json`.
Everything about which dependencies update and how stays in `renovate.json`.

```js
module.exports = {
  platform: 'github',
  repositories: ['trafalmadorian97/mecfs_bioinformatics'],
  allowedUnsafeExecutions: ['pixi'],
  gitAuthor:
    'mecfs-bio-renovate[bot] <311934930+mecfs-bio-renovate[bot]@users.noreply.github.com>',
};
```

`renovate.json` edits:

1. Widen `lockFileMaintenance.schedule` from `"before 4am on saturday"` to all-day Saturday
   (`["* * * * 6"]`). Renovate schedules are filters, not triggers: they apply only if a run
   occurs inside the window.

   At hourly cadence this is **not load-bearing** — runs at 00:00–03:00 UTC on Saturday fall
   inside the existing narrow window regardless of scheduling delays. It is kept as insurance
   against a later cadence reduction: at daily cadence a delayed run would miss a pre-04:00
   window and silently skip a week, and that interaction is easy to forget when changing the
   cron. Cost of keeping it is zero.

   Renovate accepts cron syntax here but **requires `*` in the minutes field** — verified in
   `lib/workers/repository/update/branch/schedule.ts`, which rejects anything else with
   `Invalid schedule: ... has cron syntax, but doesn't have * as minutes`. `"* * * * 6"`
   satisfies this.
2. Add a `customManager` to keep the pinned `renovate-version` updated (see below).

3. Add `"custom.regex"` to `enabledManagers`. Discovered during planning: `enabledManagers`
   "disables all other managers", and custom managers are named `custom.regex` (the config
   migration maps legacy `regex` to it). The current value `["pixi", "github-actions"]` omits
   it, so the **existing pixi-version `customManager` has never run** — `constraints.pixi` is
   pinned at `0.67.2` while pixi is at `0.75.0`, and the only commit ever touching that value
   is `da0a1a5b` (a manual bump, not a Renovate PR). Without this edit the new
   `renovate-version` customManager would be equally dead.

   Consequence: once enabled, Renovate will propose `constraints.pixi` 0.67.2 -> 0.75.0, and
   the existing `minor`/`patch` rule would automerge it. Since that changes which pixi
   resolves `pixi.lock` in CI, decide deliberately whether to let it automerge or add a
   `"automerge": false` rule for `prefix-dev/pixi` first.

`constraints.pixi`, `packageRules`, `automerge`, `enabledManagers`, `minimumReleaseAge`,
`prHourlyLimit` and the existing pixi-version `customManager` carry over unchanged.

Rationale for the strict split: dependency policy stays reviewable where it already lives,
the diff stays small, and rollback never requires rewriting policy — `renovate.json` remains
valid for both the hosted app and the self-hosted runner.

## Workflow

### Triggers

Land with `workflow_dispatch` only; the `schedule:` block is added in a **second** commit
after dry-run validation. Merging the workflow therefore cannot itself trigger a live run.

```yaml
on:
  workflow_dispatch:
    inputs:
      dryRun: { type: boolean, default: false }
      logLevel: { type: string, default: info }
  # schedule:            # second commit
  #   - cron: "0 * * * *"

concurrency: { group: renovate, cancel-in-progress: false }
permissions:
  contents: read          # checkout the config file
  pull-requests: write    # the approve step
timeout-minutes: 60
```

`permissions` constrains only `GITHUB_TOKEN`. Renovate's own rights come from the App token
and are configured on the App. The two are deliberately independent, keeping the approve
step's authority narrow.

### Steps

```yaml
- uses: actions/create-github-app-token@v3
  id: app-token
  with:
    app-id: ${{ secrets.RENOVATE_APP_ID }}
    private-key: ${{ secrets.RENOVATE_APP_PRIVATE_KEY }}

- uses: actions/checkout@v7

- uses: actions/cache@v4
  with:
    path: /tmp/renovate/cache
    key: renovate-cache-${{ github.run_id }}
    restore-keys: renovate-cache-

- uses: renovatebot/github-action@v46.2.0
  with:
    configurationFile: .github/renovate-global.js
    token: ${{ steps.app-token.outputs.token }}
    renovate-version: 44.5.3
  env:
    RENOVATE_REPOSITORY_CACHE: enabled
    RENOVATE_DRY_RUN: ${{ inputs.dryRun && 'full' || '' }}
    LOG_LEVEL: ${{ inputs.logLevel }}

- name: Approve Renovate PRs
  env: { GH_TOKEN: ${{ secrets.GITHUB_TOKEN }} }
  run: |  # filtered + idempotent
```

Caching works because the action mounts `/tmp:/tmp` by default, so the container's
`/tmp/renovate/cache` is the host's. The `run_id` key plus `restore-keys` prefix is required
rotation: `actions/cache` will not overwrite an existing key, so a static key freezes on day
one.

### Version pinning

`renovate-version` must be pinned to an exact version. The action's own default is `'44'` — a
**floating major**, which is exactly how 43.288.0's silent behaviour change would have
reached us regardless of pinning the action tag.

The `github-actions` manager understands only `uses:` refs; it has no support for
`# renovate:` annotation comments on input values (it handles `ratchet` comments only). So
the action tag auto-updates via the existing manager, but `renovate-version` needs a small
`customManager`, mirroring the one already in `renovate.json` for the pixi version.

### Approve step

Two required properties:

- **Filtered to the App's PRs.** `gh pr list --json number,author`, selecting
  `author.login == "app/mecfs-bio-renovate"`. Unfiltered, it would approve human PRs.

  Note the login formats differ by field and this is easy to get wrong: `gh` renders PR
  authors that are Apps as `app/<slug>`, but *review* authors as a plain login with neither
  the `app/` prefix nor the `[bot]` suffix (e.g. `renovate-approve`). A `<slug>[bot]` filter
  matches nothing, silently.
- **Idempotent.** It runs hourly and will see the same open PR repeatedly; re-approving
  errors. It queries existing reviews and skips PRs already approved by
  `github-actions[bot]`, rather than suppressing errors with `|| true`, which would also
  swallow real failures.

It runs inside the scheduled workflow, not on a `pull_request`/`pull_request_target`
trigger. On a public repo the latter is a known privilege-escalation footgun; running on
cron means no attacker-controlled input is in scope.

## Approvals

`renovate-approve` will **not** carry over. Its README states: "For self-hosted Renovate,
you'll need to run one or more of your own Approve bots with appropriate permissions, as
GitHub users and bots are not able to self-approve." It approves PRs from Mend's hosted
`renovate[bot]`; our PRs will be authored by `app/mecfs-bio-renovate`, a different actor. Left
unhandled, every PR would stall unapproved against the repo's one-approval requirement.

`GITHUB_TOKEN` approvals do count toward required approvals — "Allow GitHub Actions reviews
to count towards required approval" is enabled by default — and `github-actions[bot]` is a
distinct actor from the App, so this is not self-approval.

This means the one-approval gate is satisfied by a bot. That is already true today via
`renovate-approve`; this changes which bot, not whether the gate is real.

Renovate's own automerge only runs when Renovate runs. GitHub native auto-merge
(`platformAutomerge`, default on) sidesteps this, but `lockFileMaintenance` sets
`platformAutomerge: false`, so those weekly PRs wait for the next Renovate run — about an
hour at hourly cadence. Worth revisiting if cadence is ever reduced.

## Cutover

### Phase 0 — GitHub UI (manual)

**Status: App registered as `mecfs-bio-renovate` (bot user id `311934930`).** Remaining:
confirm it is installed on `mecfs_bioinformatics`, and add the two repo secrets.

Register the App: webhook **unchecked** (cron-driven, nothing listens), installable on this
account only. Permissions:

| Permission | Level | Why |
|---|---|---|
| Contents | RW | push branches, commit `pixi.lock` |
| Pull requests | RW | open/update/merge PRs |
| Issues | RW | dependency dashboard |
| Commit statuses | RW | read CI results for automerge |
| Checks | RW | same |
| Workflows | RW | **required** — the `github-actions` manager edits `.github/workflows/` |
| Administration | read | read branch protection to decide automerge |
| Metadata | read | mandatory baseline |

Generate a private key, install on `mecfs_bioinformatics`, add `RENOVATE_APP_ID` and
`RENOVATE_APP_PRIVATE_KEY` as repo secrets.

`gitAuthor` was derived rather than hand-typed:
`gh api /users/mecfs-bio-renovate%5Bbot%5D --jq .id` returned `311934930`, giving
`311934930+mecfs-bio-renovate[bot]@users.noreply.github.com`. This is the design's most
error-prone value; a mismatch silently disables automerge.

The private key must be pasted directly into the repo secret by the repo owner. It should
not be transmitted through any other channel.

### Phases 1–4

1. **Land dispatch-only.** One PR: `renovate.yml` (no `schedule:`), `renovate-global.js`,
   `renovate.json` edits. Merging cannot trigger anything.
2. **Dry run.** Dispatch with `dryRun: true, logLevel: debug`. The key check is a negative:
   the **absence** of `` `pixi lock` was requested to run, but `pixi` is not permitted in the
   allowedUnsafeExecutions `` in the log. That string is the exact failure being fixed, so
   its absence proves the fix landed.
3. **Cut over.** Uninstall Mend Renovate *and* `renovate-approve`, then a second PR adds the
   cron. In that order — two live bots means duplicate PRs on identical branch names.
4. **Verify.** Watch the first scheduled run end to end: PR opens, is approved, checks pass,
   merges.

**Prerequisite: zero open Renovate PRs at cutover.** Renovate finds PRs by branch name, so
the new App would adopt existing `renovate/*` branches — but those PRs are authored by Mend's
`renovate[bot]`, and the approve step filters on our App's login, so they would never be
approved and would sit indefinitely.

- **#983** (ty 0.0.65) — merged 2026-08-01. Cleared.
- **#955** (Lock file maintenance, branch `renovate/lock-file-maintenance`) — **close it.**
  Open since 2026-07-25, `mergeStateStatus: BEHIND`, last checks green on 2026-07-25. It
  cannot self-heal: `updatePixiLockfile` returns `null` at the `allowedUnsafeExecutions`
  gate *before* doing any work, and a lock-maintenance PR's only content is that lock
  update — so a rebase under the current hosted app produces nothing. Once self-hosting is
  live, the next Saturday window generates a fresh, up-to-date replacement, which strictly
  dominates a week-stale one.

  Merging it instead would be safe but pointless: its diff is 589 lines of `pixi.lock` only
  and does not touch `ty`, so there is no risk of reverting the recent bump — it is simply
  redundant with the refresh that self-hosting will produce.

## Failure modes

| Symptom | Cause | Detection |
|---|---|---|
| PRs open, never merge | `gitAuthor` does not match App identity | `branch.isModified() = true` in debug log |
| `pixi.lock` still stale | `allowedUnsafeExecutions` not applied | the `logger.once.warn` string above |
| Action-bump PRs fail | missing Workflows permission | 403 in log |
| Human PRs approved | approve-step filter bug | test the `gh pr list` filter standalone before wiring it up |
| Cache step fails | container UID vs `runner` | drop caching; do **not** escalate to `docker-user: root` |

The existing `prHourlyLimit: 10` caps the blast radius if a cold cache on the first run is
mistaken for a pile of new work.

## Testing

Workflows cannot be unit-tested, but `on_pr` already runs **actionlint**, which validates
workflow YAML and shellchecks inline `run:` blocks — so the approve script is linted on the
PR that introduces it. Beyond that, verification is the dry run plus the first live run.

## Rollback

Fast path: disable the workflow in the Actions UI, reinstall the Mend app and
`renovate-approve`. Back to current behaviour (including the manual `pixi.lock` chore) in
about two minutes.

Full path: revert both PRs, delete the secrets, uninstall the App.

Reversibility comes from the config split: `renovate.json` stays valid for both the hosted
app and the self-hosted runner, so rolling back never means rewriting dependency policy.

## Accepted risk

`allowedUnsafeExecutions: ['pixi']` means conda package hooks — code from dependencies —
execute in a job holding a repo-write token. Bounded by the 1-hour token lifetime, the
single-repo installation, and the permission set above, and comparable to `on_pr` already
running `pixi install`. It is nonetheless deliberate un-gating of something upstream gated,
and is accepted knowingly.

## Watch for

Upstream may relax this gate: `mise` was gated identically and then got safe-mode locking in
43.282.0 (renovatebot/renovate#44749). No pixi equivalent is filed. If one lands, this whole
setup could be retired in favour of the hosted app.

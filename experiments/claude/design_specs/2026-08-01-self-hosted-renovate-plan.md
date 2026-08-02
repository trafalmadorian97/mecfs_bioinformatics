# Self-Hosted Renovate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run Renovate ourselves via GitHub Actions with the `mecfs-bio-renovate` GitHub App, so `allowedUnsafeExecutions: ['pixi']` can be set and `pixi.lock` is updated again.

**Architecture:** An hourly workflow mints a 1-hour GitHub App installation token, runs the Renovate container against a checked-in global config, then approves the resulting PRs with `GITHUB_TOKEN`. `renovate.json` remains the source of truth for dependency policy.

**Tech Stack:** GitHub Actions, `renovatebot/github-action@v46.2.0`, `actions/create-github-app-token@v3`, Renovate 44.5.3, `gh` CLI, actionlint.

Spec: `experiments/claude/design_specs/2026-08-01-self-hosted-renovate-design.md`

## Global Constraints

- App: `mecfs-bio-renovate`, bot user id `311934930`.
- `gitAuthor` must be exactly `mecfs-bio-renovate[bot] <311934930+mecfs-bio-renovate[bot]@users.noreply.github.com>`. A mismatch silently disables automerge via `isBranchModified()`.
- `renovate-version` pinned to an exact version. The action's default `'44'` is a floating major and must not be relied on.
- Renovate `schedule` cron requires `*` in the minutes field.
- Never escalate the Renovate container to `docker-user: root`. If caching fails, delete the cache step.
- The App private key is pasted directly into the repo secret by the repo owner and transmitted nowhere else.
- All repo commands run via `pixi r` per CLAUDE.md.
- Cutover requires zero open Renovate PRs (#983 merged, #955 closed — already satisfied).

---

## File Structure

| File | Responsibility |
|---|---|
| `.github/renovate-global.js` (create) | Self-hosted-only config: platform, repo list, `allowedUnsafeExecutions`, `gitAuthor`. Nothing about dependency policy. |
| `.github/workflows/renovate.yml` (create) | Trigger, token minting, Renovate invocation, PR approval. |
| `renovate.json` (modify) | Dependency policy. Three edits: enable `custom.regex`, widen the Saturday window, add the `renovate-version` customManager. |

---

### Task 1: Repo secrets and App installation

Manual GitHub UI work. Nothing later works without it.

**Files:** none.

**Interfaces:**
- Produces: repo secrets `RENOVATE_APP_ID`, `RENOVATE_APP_PRIVATE_KEY`; App installed on `trafalmadorian97/mecfs_bioinformatics`.

- [ ] **Step 1: Confirm the App is installed on the repo**

Registering and installing are separate. Go to
`https://github.com/settings/apps/mecfs-bio-renovate/installations`.
It must show `mecfs_bioinformatics` under "Install App". If it does not, click Install and
select **Only select repositories** → `mecfs_bioinformatics`.

- [ ] **Step 2: Copy the App ID**

On `https://github.com/settings/apps/mecfs-bio-renovate`, under "About", copy the numeric
**App ID** (a 6–8 digit number). This is *not* the `311934930` bot user id — that is a
different number and is already baked into `gitAuthor`.

- [ ] **Step 3: Add `RENOVATE_APP_ID`**

Repo → Settings → Secrets and variables → Actions → New repository secret.
Name: `RENOVATE_APP_ID`. Value: the number from Step 2, digits only, no quotes or spaces.

- [ ] **Step 4: Generate a private key**

On the App settings page, scroll to "Private keys" → **Generate a private key**. A `.pem`
file downloads. Treat it as a credential: it is the App's root secret and does not expire.

- [ ] **Step 5: Add `RENOVATE_APP_PRIVATE_KEY`**

New repository secret. Name: `RENOVATE_APP_PRIVATE_KEY`. Value: the **entire** contents of
the `.pem` file, including both delimiter lines:

```
-----BEGIN RSA PRIVATE KEY-----
MIIEow...
...many lines...
-----END RSA PRIVATE KEY-----
```

This is the step that most often goes wrong, and the failure is confusing rather than
obvious — a truncated or reflowed key produces `error:1E08010C:DECODER routines::unsupported`
at token-minting time, not a clear "bad key" message. Open the `.pem` in a plain text editor,
select all, copy. Do not retype it, do not strip the `BEGIN`/`END` lines, and do not worry
about the trailing newline (GitHub handles it).

- [ ] **Step 6: Verify both secrets exist**

```bash
gh secret list --repo trafalmadorian97/mecfs_bioinformatics
```

Expected: rows for `RENOVATE_APP_ID` and `RENOVATE_APP_PRIVATE_KEY`. Values are never
displayed — presence is all that can be checked here. Correctness is proven in Task 5.

- [ ] **Step 7: Delete the local `.pem`**

```bash
rm ~/Downloads/mecfs-bio-renovate.*.private-key.pem
```

It is now in the secret; a copy sitting in Downloads is pure liability. If it is ever
needed again, generate a fresh key and update the secret.

---

### Task 2: `renovate.json` edits and the global config

**Files:**
- Modify: `renovate.json`
- Create: `.github/renovate-global.js`

**Interfaces:**
- Produces: `.github/renovate-global.js` (consumed by Task 3's `configurationFile` input); `custom.regex` enabled so Task 3's `renovate-version` pin is auto-updatable.

- [ ] **Step 1: Enable `custom.regex` in `enabledManagers`**

`enabledManagers` disables every manager not listed, and custom managers are named
`custom.regex`. Its omission means the existing pixi customManager has never run.

In `renovate.json`, replace lines 8–11:

```json
  "enabledManagers": [
    "pixi",
    "github-actions"
  ],
```

with:

```json
  "enabledManagers": [
    "pixi",
    "github-actions",
    "custom.regex"
  ],
```

- [ ] **Step 2: Widen the lock-maintenance window**

Replace `"schedule": ["before 4am on saturday"],` with:

```json
    "schedule": ["* * * * 6"],
```

Not load-bearing at hourly cadence; insurance against a later cadence reduction. Renovate
rejects cron here unless minutes is `*`.

- [ ] **Step 3: Add the `renovate-version` customManager**

In the `customManagers` array, after the existing pixi entry (line 25's `}`), add a comma and:

```json
    {
      "customType": "regex",
      "managerFilePatterns": ["/^\\.github/workflows/renovate\\.yml$/"],
      "matchStrings": [
        "renovate-version:\\s*(?<currentValue>\\d+\\.\\d+\\.\\d+)"
      ],
      "datasourceTemplate": "docker",
      "depNameTemplate": "ghcr.io/renovatebot/renovate"
    }
```

- [ ] **Step 4: Create `.github/renovate-global.js`**

```js
// Self-hosted (global) Renovate config. Dependency policy lives in renovate.json.
// Only options that CANNOT be set in renovate.json belong here.
module.exports = {
  platform: 'github',
  repositories: ['trafalmadorian97/mecfs_bioinformatics'],

  // `pixi lock` executes conda package hooks, which Renovate classes as an unsafe
  // execution and gates behind this option. Without it, pyproject.toml is updated but
  // pixi.lock is left stale and CI's `pixi install --locked` fails.
  allowedUnsafeExecutions: ['pixi'],

  // Must match the App's identity exactly, or isBranchModified() flags Renovate's own
  // branches as externally modified and automerge silently stops.
  gitAuthor:
    'mecfs-bio-renovate[bot] <311934930+mecfs-bio-renovate[bot]@users.noreply.github.com>',
};
```

- [ ] **Step 5: Validate the config**

```bash
docker run --rm -v "$PWD:/repo" -w /repo \
  ghcr.io/renovatebot/renovate:44.5.3 renovate-config-validator
```

Expected: `INFO: Config validated successfully`. If it reports an invalid schedule, the cron
minutes field is wrong (Step 2).

- [ ] **Step 6: Commit**

```bash
git add renovate.json .github/renovate-global.js
git commit -m "Add self-hosted Renovate global config

Enable custom.regex in enabledManagers: it was omitted, which silently
disabled the existing pixi-version customManager (constraints.pixi is
stuck at 0.67.2 while pixi is at 0.75.0)."
```

---

### Task 3: The workflow

**Files:**
- Create: `.github/workflows/renovate.yml`

**Interfaces:**
- Consumes: `RENOVATE_APP_ID`, `RENOVATE_APP_PRIVATE_KEY` (Task 1); `.github/renovate-global.js` (Task 2).
- Produces: a `workflow_dispatch`-triggerable workflow with `dryRun` and `logLevel` inputs. The `schedule:` block is deliberately absent until Task 6.

- [ ] **Step 1: Create the workflow**

```yaml
name: renovate

on:
  workflow_dispatch:
    inputs:
      dryRun:
        description: "Run without creating branches, PRs or merges"
        type: boolean
        default: false
      logLevel:
        description: "Renovate log level (info or debug)"
        type: string
        default: info
  # schedule: added in Task 6, after dry-run validation.

concurrency:
  group: renovate
  cancel-in-progress: false

permissions:
  contents: read
  pull-requests: write

jobs:
  renovate:
    name: Renovate
    runs-on: ubuntu-latest
    timeout-minutes: 60
    steps:
      - name: Mint App token
        id: app-token
        uses: actions/create-github-app-token@v3
        with:
          app-id: ${{ secrets.RENOVATE_APP_ID }}
          private-key: ${{ secrets.RENOVATE_APP_PRIVATE_KEY }}

      - name: Check out repository
        uses: actions/checkout@v7

      - name: Restore Renovate cache
        uses: actions/cache@v4
        with:
          path: /tmp/renovate/cache
          key: renovate-cache-${{ github.run_id }}
          restore-keys: renovate-cache-

      - name: Run Renovate
        uses: renovatebot/github-action@v46.2.0
        with:
          configurationFile: .github/renovate-global.js
          token: ${{ steps.app-token.outputs.token }}
          renovate-version: 44.5.3
        env:
          RENOVATE_REPOSITORY_CACHE: enabled
          RENOVATE_DRY_RUN: ${{ inputs.dryRun && 'full' || '' }}
          LOG_LEVEL: ${{ inputs.logLevel }}

      - name: Approve Renovate PRs
        if: ${{ !inputs.dryRun }}
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
        run: |
          set -euo pipefail

          pr_numbers=$(gh pr list --state open --json number,author \
            --jq '.[] | select(.author.login == "mecfs-bio-renovate[bot]") | .number')

          if [ -z "$pr_numbers" ]; then
            echo "No open PRs from mecfs-bio-renovate[bot]."
            exit 0
          fi

          while IFS= read -r pr; do
            [ -n "$pr" ] || continue
            already=$(gh pr view "$pr" --json reviews \
              --jq '[.reviews[]
                     | select((.author.login | startswith("github-actions"))
                              and .state == "APPROVED")] | length')
            if [ "$already" -gt 0 ]; then
              echo "PR #${pr}: already approved, skipping."
              continue
            fi
            echo "PR #${pr}: approving."
            gh pr review "$pr" --approve
          done <<< "$pr_numbers"
```

The approve step is skipped on dry runs because there is nothing to approve, and it is
filtered to the App's login so it can never approve a human PR.

- [ ] **Step 2: Lint the workflow**

```bash
pixi r invoke lint-actions
```

actionlint validates the YAML and shellchecks the inline `run:` block. Expected: no output,
exit 0. Fix any `SC####` findings before continuing — this is the only automated check this
workflow will ever get.

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/renovate.yml
git commit -m "Add self-hosted Renovate workflow (dispatch-only)"
```

---

### Task 4: Verify the approve filter before it can ever run

The filter is the only thing preventing the workflow from approving human PRs. Test it
standalone, against the live repo, before it has any authority.

**Files:** none (verification only).

- [ ] **Step 1: Confirm the filter selects nothing right now**

```bash
gh pr list --state open --json number,author \
  --jq '.[] | select(.author.login == "mecfs-bio-renovate[bot]") | .number'
```

Expected: empty. No PRs exist from the App yet, so anything else means the filter is wrong.

- [ ] **Step 2: Confirm the filter is actually discriminating**

```bash
gh pr list --state all --limit 20 --json number,author \
  --jq '.[] | [.number, .author.login] | @tsv'
```

Read the logins. Confirm the historical Renovate PRs show `renovate[bot]` (Mend's, which the
filter must **not** match) and your own PRs show your login. This proves Step 1's empty
result is a real filter and not a broken query returning nothing regardless of input.

- [ ] **Step 3: Confirm the already-approved query parses**

```bash
gh pr view 983 --json reviews \
  --jq '[.reviews[] | select((.author.login | startswith("github-actions")) and .state == "APPROVED")] | length'
```

Expected: `0` (#983 was approved by `renovate-approve`, not `github-actions`). A jq syntax
error here would, in the workflow, make `already` empty and crash the `-gt` comparison under
`set -u`.

---

### Task 5: Merge and dry-run validation

**Files:** none (validation only).

**Interfaces:**
- Consumes: everything from Tasks 1–3.

- [ ] **Step 1: Open and merge the PR**

```bash
git push -u origin self-hosted-renovate
gh pr create --title "Self-hosted Renovate (dispatch-only)" \
  --body "Spec: experiments/claude/design_specs/2026-08-01-self-hosted-renovate-design.md

Adds the workflow without a schedule trigger, so merging cannot start a live run.
The cron is added separately after dry-run validation."
```

Merge it once `on_pr` is green. No live run can result — there is no `schedule:` trigger.

- [ ] **Step 2: Dispatch a dry run**

```bash
gh workflow run renovate.yml -f dryRun=true -f logLevel=debug
```

- [ ] **Step 3: Watch it**

```bash
gh run watch "$(gh run list --workflow=renovate.yml --limit 1 --json databaseId --jq '.[0].databaseId')"
```

- [ ] **Step 4: Check the four things that matter**

```bash
gh run view --log "$(gh run list --workflow=renovate.yml --limit 1 --json databaseId --jq '.[0].databaseId')" > /tmp/renovate-dryrun.log
```

1. **Token minted and repo found** — `grep -c "mecfs_bioinformatics" /tmp/renovate-dryrun.log` is non-zero.
2. **The gate is open** — the decisive check, and it is a *negative*:
   ```bash
   grep -i "not permitted in the allowedUnsafeExecutions" /tmp/renovate-dryrun.log
   ```
   Expected: **no match.** That string is the exact failure being fixed; its absence proves
   `allowedUnsafeExecutions` took effect.
3. **`gitAuthor` accepted** — `grep -i "invalid gitAuthor\|gitAuthor" /tmp/renovate-dryrun.log` shows no error. Renovate validates the format at startup.
4. **Managers ran** — `grep -iE "custom.regex|pixi|github-actions" /tmp/renovate-dryrun.log` shows all three extracting.

- [ ] **Step 5: Expect a pixi 0.75 proposal**

With `custom.regex` newly enabled, Renovate will now see that `constraints.pixi` (`0.67.2`)
is far behind pixi `0.75.0`. In the dry-run log it will report this as a pending update.
That is correct and expected — it is a bug being fixed, not a malfunction. Note that under
the existing `minor`/`patch` automerge rule it **will automerge** once live. Decide before
Task 6 whether that is acceptable; if not, add a `packageRules` entry setting
`"automerge": false` for `prefix-dev/pixi` first. Changing the pixi version changes which
pixi resolves `pixi.lock` in CI, so this is worth a deliberate choice.

---

### Task 6: Cutover

**Files:**
- Modify: `.github/workflows/renovate.yml`

- [ ] **Step 1: Uninstall the Mend apps**

At `https://github.com/settings/installations`, uninstall **Renovate** and
**Renovate Approve** from `mecfs_bioinformatics`. Do this *before* enabling the cron: two
live bots target identical `renovate/*` branch names and will duplicate and fight.

- [ ] **Step 2: Confirm no Renovate PRs appeared in the meantime**

```bash
gh pr list --state open --json number,author --jq '.[] | [.number, .author.login] | @tsv'
```

Expected: no `renovate[bot]` rows. If any exist, close them — the new App's approve filter
will never approve them and they would sit forever.

- [ ] **Step 3: Enable the cron**

In `.github/workflows/renovate.yml`, replace:

```yaml
  # schedule: added in Task 6, after dry-run validation.
```

with:

```yaml
  schedule:
    - cron: "0 * * * *"
```

- [ ] **Step 4: Lint, commit, PR, merge**

```bash
pixi r invoke lint-actions
git add .github/workflows/renovate.yml
git commit -m "Enable hourly Renovate schedule"
```

Open and merge a PR for this as normal.

---

### Task 7: Verify the first live run

**Files:** none.

- [ ] **Step 1: Wait for or force a run**

```bash
gh workflow run renovate.yml
gh run watch "$(gh run list --workflow=renovate.yml --limit 1 --json databaseId --jq '.[0].databaseId')"
```

- [ ] **Step 2: Confirm the PR identity is right**

```bash
gh pr list --state open --json number,author,title --jq '.[] | [.number, .author.login, .title] | @tsv'
```

Expected: any new PRs show `mecfs-bio-renovate[bot]`. If they show your login, `gitAuthor` or
the token is wrong.

- [ ] **Step 3: Confirm `pixi.lock` is included**

For any PR touching `pyproject.toml`:

```bash
gh pr diff <N> --name-only
```

Expected: **both** `pyproject.toml` and `pixi.lock`. `pyproject.toml` alone means the
`allowedUnsafeExecutions` gate is still closed — the original bug, unfixed.

- [ ] **Step 4: Confirm approval and automerge**

```bash
gh pr view <N> --json reviews,autoMergeRequest \
  --jq '{approvedBy: [.reviews[] | select(.state == "APPROVED") | .author.login], autoMerge: .autoMergeRequest}'
```

Expected: `github-actions` in `approvedBy`, and `autoMerge` non-null. Then confirm it
actually merges once `on_pr` goes green.

- [ ] **Step 5: If a PR opens but never merges**

Symptom of a `gitAuthor` mismatch. Confirm with:

```bash
gh run view --log <run-id> | grep -i "isModified\|branch.isModified"
```

`branch.isModified() = true` means Renovate does not recognise its own commits. Compare the
`gitAuthor` in `.github/renovate-global.js` against
`gh api /users/mecfs-bio-renovate%5Bbot%5D --jq .id` (must be `311934930`).

---

## Rollback

At any point: disable the workflow (Actions tab → renovate → ⋯ → Disable workflow), then
reinstall the Mend Renovate and Renovate Approve apps. That restores current behaviour,
including the manual `pixi.lock` chore, in about two minutes. `renovate.json` stays valid for
both, so rollback never requires rewriting dependency policy.

Full removal: revert the Task 2/3/6 commits, delete both secrets, uninstall the App.

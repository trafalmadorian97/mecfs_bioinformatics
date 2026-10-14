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

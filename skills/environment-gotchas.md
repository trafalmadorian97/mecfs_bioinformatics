# Environment gotchas

Recurring dependency and tool traps a contributor is likely to hit, with the fix already in place
and how to reason about it if it resurfaces. Most fixes live in `pyproject.toml` under pixi's
dependency sections.

## Dependency pins (why they exist — don't drop them blindly)

- **`rpy2` must come from conda-forge, not PyPI.** It is declared in
  `[tool.pixi.feature.analysis.dependencies]` (conda), not `pypi-dependencies`. PyPI publishes
  `rpy2-rinterface` as an sdist only, so pixi compiles it against the env's R; that build's link
  step needs the unversioned `lib/liblzma.so` dev symlink, which ships only in `liblzma-devel` — a
  package that fell out of the solve on the Python 3.13 upgrade. conda-forge ships rpy2 prebuilt
  against a matching r-base, so nothing compiles. If a future R/Python bump makes conda rpy2
  unavailable for the target combo, that is the constraint to solve — not the linker error.
- **`pysam` is overridden to `>=0.24,<0.25`** (`[tool.pixi.pypi-options.dependency-overrides]`).
  gwaslab pins `pysam==0.22.1`, which has no cp313 wheel; on Python 3.13 uv builds the sdist, which
  re-cythonizes and breaks `import pysam` with `AttributeError: ... has no attribute 'CMATCH'`,
  surfacing as ~10 pytest *collection* errors. The PyPI pysam fully shadows the conda one
  (`pixi list pysam` shows two rows — the tell). gwaslab only touches the stable
  `VariantFile`/tabix surface, so relaxing the pin is safe. General lesson: when a PyPI dep pins an
  extension package exactly and you bump Python, check `pixi list <pkg>` for a duplicate conda/pypi
  pair before believing the package itself is broken.
- **`setuptools` is pinned `>=80.10.1,<82`** (`[tool.pixi.pypi-dependencies]`). setuptools 82.0.0
  removed `pkg_resources`, which `zepid/datasets/__init__.py` imports unconditionally at import
  time, failing collection on every test that transitively imports `plot_mr_effect_measure_task`.
  If a future lockfile bump reintroduces `ModuleNotFoundError: No module named 'pkg_resources'`,
  check whether the pin was dropped or a dep forced setuptools past it. The durable alternative is
  to stop importing zepid (used only for `EffectMeasurePlot`).
- **`numpy` is capped `>=2.1,<2.4`** (dependency-override + a matching `renovate.json` rule with
  `allowedVersions: "<2.4"`). gwaslab's vendored LDSC rg estimator calls `float()` on size-1
  `ndim>0` arrays; numpy 2.4.0 turned that long-deprecated conversion into a hard `TypeError`,
  crashing the CT-LDSC genetic-correlation path. Lift the cap only once gwaslab extracts the scalar
  with `.item()` upstream — bumping to numpy 2.4+ before that reintroduces the crash; downgrading
  below 2.0 isn't an option on py3.13 (no numpy 1.x cp313 wheels).

## gwaslab runtime traps

- **gwaslab's rg error handler masks the real exception.** When cross-trait LDSC fails inside
  gwaslab, the *final* exception (`TypeError: '>=' not supported between instances of 'LinAlgError'
  and 'int'`) is noise from a logging bug (`estimate_rg` passes an exception where
  `traceback.format_exc` expects a `limit`). **Read the FIRST traceback in the chain** — that holds
  the real failure (e.g. `LinAlgError: SVD did not converge`, usually from non-finite Z). Bare
  LAPACK lines like `On entry to DLASCL parameter number 4 had an illegal value` mean NaN/Inf
  reached `np.linalg.lstsq` — check Z, BETA/SE, and N for non-finite values.
- **gwaslab 4.x `harmonize()` OOMs on datasets with long alleles.** It stores EA/NEA as pandas
  `category` dtype, and a whole-column `Categorical.astype(str)` in `_flip_allele_stats` materializes
  a fixed-width unicode array sized to the *longest* category over every row (rows × maxlen × 4
  bytes) — one 662 bp indel across ~17M rows blows up memory. MAF pre-filtering is the wrong lever
  (blowup scales with allele length, not frequency). The fix has **two required parts**: (1)
  `filter_indels=True` before harmonization, and (2) `.cat.remove_unused_categories()` — because
  filtering drops the *rows* but a categorical keeps its full `.categories` index, so the long
  category (and the `<UNNN>` width) survives until you prune it.

## Filesystem / platform

- **An asset store on a WSL DrvFs mount (e.g. `/mnt/d`) breaks `shutil.copy2` cross-device moves.**
  Finalizing an asset from a `/tmp` (ext4) scratch dir to a DrvFs asset root is a cross-device move;
  `shutil.copy2`/`copytree` call `copystat` (utime/chmod), which DrvFs forbids
  (`PermissionError: Operation not permitted`) — *after* the expensive task already produced its
  output, so the result is discarded. Fix: copy data only (`shutil.copyfile` for files, manual
  `os.walk` + `copyfile` for dirs), keeping the atomic `os.replace` swap (which works on DrvFs).
  Asset integrity is tracked by the content trace, not filesystem metadata, so dropping
  timestamps/perms is safe. Beyond this bug, DrvFs mounts have high per-file latency, so directory
  assets with thousands of small files are painfully slow there — prefer local ext4 for those.

## SBayesRC / polypwas run via Docker, not local conda R

SBayesRC runs via the `zhiliz/sbayesrc` Docker image (Ubuntu 24.04), not a local conda-R install,
because the SBayesRC R package fails to compile against this env's gcc 15.2 (a Boost.Math overload
mismatch), and the image also sidesteps an Ubuntu 22.04 OpenBLAS bug that corrupts its
eigen-decomposition. `polypwas train` shells out to SBayesRC via a Docker `Rscript` shim; `polypwas
assoc` needs no R. When adding related tasks, mirror the MiXeR subsystem's Docker invocation shape
rather than adding conda R dependencies.

# GWFM known-truth toy recipe — proven end-to-end (Task 12 recon)

Date: 2026-08-09. Author: controller spike driving real gctb 2.5.5 (static Linux binary
from the pinned GCTB_BINARY_URL). Everything below was executed and verified; the working
generator is experiments/claude/gwfm_toy_spike/gen_toy_gwfm.py and reproduces it from
scratch. This addendum supplements 2026-08-08-gwfm-recon-findings.md.

## Result

A tiny synthetic dataset (3 LD blocks x 150 SNPs, 2000 individuals) drives the full
four-step gctb GWFM pipeline to completion in ~4 s and recovers all three planted causal
SNPs (rs11, rs171, rs331) each as its own 1-SNP credible set with PIP = 1.0
(gwfm.snpRes PIP column and gwfm.lcs local credible sets). All four gctb invocations exit 0.

## The four commands (the same three templates Task 10 assembles, plus the reference build)

```
# 0. build block LD reference from PLINK genotypes (this is the toy stand-in for the
#    staged ukbEUR_13M_FullLDM.zip -> ldm13M/ bundle)
gctb --bfile toy --make-block-ldm --block-info blocks.txt --out ldm --thread 1
# 1. GCTB_MAKE_LDM_EIGEN_TEMPLATE
gctb --ldm ldm --gwas-summary toy.ma --make-ldm-eigen --ldm-eigen-cutoff 0.995 --thread 1 --out matched
# 2. GCTB_GWFM_TEMPLATE
gctb --gwfm RC --ldm-eigen matched --gwas-summary toy.ma --annot annot.txt --gene-map genemap.txt --thread 1 --out gwfm --chain-length 3000 --burn-in 1000
# 3. GCTB_CS_TEMPLATE
gctb --cs --pwld-file ldm/rsq0.5.pwld --pip 0.9 --pep 0.7 --gene-map genemap.txt --mcmc-samples gwfm --out gwfm
```

With >=2 blocks, `--make-block-ldm` writes the folder-level `snp.info`, `ldm.info`,
`rsq0.5.pwld` plus `block<N>.ldm.bin` directly — exactly the `ldm13M/` structure Task 10
expects. No separate `--merge-block-ldm-info` is needed (it errors "only one info file"
on a single block; for >=2 blocks the merge is implicit).

## Input file formats (confirmed against gctb source github.com/jianzeng/GCTB)

- **.ma** (COJO): `SNP A1 A2 freq b se p N`. For known truth, compute the LD matrix
  `R = corrcoef(genotypes)` and set `b_marginal = R @ b_true` (+ tiny noise), `se = 1/sqrt(N)`.
- **.bim column 3 (genetic position, cM) MUST be nonzero and increasing.** All-zero cM
  yields `GenPos = 0` for every SNP in snp.info; downstream tempered-Gibbs genetic-distance
  math then misbehaves. We use `cM = bp * 1e-5`.
- **block-info** (`--block-info`): header `Block Chr StartBP EndBP`, one row per block.
- **annotation** (`--annot`): header `SNP Intercept Anno1 ...`; `Intercept = 1` for every
  SNP (gctb requires every SNP to carry the intercept annotation), plus 0/1 columns.
- **gene-map** (`--gene-map`): EXACT 9-column header
  `Ensgid GeneName GeneType Chrom_hg38 Start_hg38 End_hg38 Chrom_hg19 Start_hg19 End_hg19`.

## Failure modes discovered (all were real crashes; each has a definitive root cause)

1. **SIGFPE in `Data::inputPairwiseLD` (data.cpp:6370, `sum/cnt`).** gwfm-RC calls
   `inputPairwiseLD(eigenDir+"/rsq0.5.pwld", 0.95)` with a HARDCODED rsq threshold of 0.95
   (main.cpp:228, "for TGS sampling"). It keeps only pairs with `ldcor^2 > 0.95`
   (|r| > 0.975); `cnt` counts SNPs with >=1 such friend. If the toy has NO pair that
   strong, `cnt = 0` and `sum/cnt` is integer division by zero -> SIGFPE. **Fix: inject at
   least one near-duplicate SNP pair (|r| > 0.975).** The generator copies a neighbor with
   ~0.5% random flips at three non-causal positions, yielding 3 rsq>0.95 pairs.
2. **SIGSEGV in `Data::readGeneMapFile` (data.cpp:446).** A gene-map missing any of the 9
   expected column names makes `header.getIndex("Ensgid")` return -1, and `colData[-1]`
   reads out of bounds -> segfault. **Fix: the exact 9-column header above.**
3. **make-ldm-eigen threading nondeterminism.** With `--thread 2` one spike run computed
   only block 1's eigen and skipped the folder-level merge; `--thread 1` was deterministic.
   Toy generation uses `--thread 1` (data is tiny, speed is irrelevant).

## --cs `--pwld-file` path quirk (note for Task 10 / the real run)

In `--cs` mode gctb reads the pwld as `opt.eigenMatrixFile + "/" + opt.pairwiseLDfile`
(main.cpp:373). In a bare `--cs` call `eigenMatrixFile` is empty, so `--pwld-file
ldm/rsq0.5.pwld` is read as the ABSOLUTE path `/ldm/rsq0.5.pwld` and fails to open. In the
spike the credible sets (gwfm.lcs/gwfm.gcs, incl. multi-SNP sets) were still produced
correctly from gwfm.snpRes despite this — the pwld only refines correlated-SNP grouping.
Task 10's step-3 template passes `--pwld-file work/ref/ldm13M/rsq0.5.pwld`, which becomes
`/work/ref/ldm13M/rsq0.5.pwld` (absolute). Whether that resolves depends on the container
mount layout; the LocalDockerRemoteExecutor stages under `/work/work/...` (see below), so
the absolute `/work/ref/...` will NOT match there. **The known-truth assertion should key
off gwfm.snpRes PIP (robust to this quirk), not off pwld-dependent grouping.**

## LocalDockerRemoteExecutor directory-output bug (MUST fix in Task 12)

`local_docker_remote_executor.py` retrieves each `output_files` entry with
`shutil.copy(staged_output, final_dest)` (line ~76). Task 10's output (`work/out`) is a
DIRECTORY, and `shutil.copy` raises `IsADirectoryError` on a directory. This mirrors the
FakeRemoteExecutor dir-aware fix already made in Task 10. **Fix: if `staged_output.is_dir()`
use `shutil.copytree`, else `shutil.copy`.** Also note the executor mounts `host_work_dir`
at `/work` and stages `input_files`/`s3_inputs` under their remote_dest verbatim; since
Task 10's remote_dest paths already start with `work/`, container paths end up under
`/work/work/...` — consistent because the gctb commands are relative to cwd `/work`.

## Reproduce

`pixi r python experiments/claude/gwfm_toy_spike/gen_toy_gwfm.py <outdir>` then run the four
commands above inside `<outdir>`. Verified 2026-08-09: 4/4 exit 0; rs11/rs171/rs331 PIP=1.0.

"""PROVEN toy generator for the GWFM known-truth system test (spike artifact).

Produces a tiny PLINK bfile + block-info + planted-causal .ma + annot + gene-map that
drive the full four-step gctb GWFM pipeline (make-block-ldm, make-ldm-eigen, gwfm RC, cs)
to completion in seconds, recovering the three planted causal SNPs (rs11, rs171, rs331)
each as its own 1-SNP credible set with PIP = 1.0.

Non-obvious requirements discovered by driving real gctb 2.5.5 (see recon addendum
2026-08-09-gwfm-toy-recipe.md):
  - .bim column 3 (genetic position, cM) MUST be nonzero/increasing, else GenPos=0 in
    snp.info (tempered-Gibbs genetic-distance math).
  - The toy MUST contain at least one SNP pair with rsq > 0.95 (|r| > 0.975), else the
    gwfm-RC tempered-Gibbs pwld reader (data.cpp inputPairwiseLD, hardcoded rsq 0.95)
    divides sum/cnt with cnt=0 -> SIGFPE. We inject near-duplicate tag SNPs.
  - The gene-map MUST have the exact 9-column header
    Ensgid GeneName GeneType Chrom_hg38 Start_hg38 End_hg38 Chrom_hg19 Start_hg19 End_hg19,
    else readGeneMapFile getIndex returns -1 -> out-of-bounds -> SIGSEGV.
  - make-ldm-eigen must run with --thread 1 (threaded runs nondeterministically dropped a
    block's eigen output in the spike).
"""

import sys
from math import erf, sqrt
from pathlib import Path

import numpy as np

outdir = Path(sys.argv[1])
outdir.mkdir(parents=True, exist_ok=True)
pref = outdir / "toy"
rng = np.random.default_rng(7)

n_ind = 2000
n_block = 3
spb = 150  # SNPs per block
n_snp = n_block * spb
chrom = 1

freq = rng.uniform(0.25, 0.5, n_snp)
latent = rng.standard_normal((n_ind, n_snp))
for j in range(1, n_snp):
    if j % spb == 0:
        continue  # block boundary breaks LD
    latent[:, j] = 0.9 * latent[:, j - 1] + np.sqrt(1 - 0.81) * latent[:, j]

geno = np.zeros((n_ind, n_snp), dtype=np.int8)
for j in range(n_snp):
    lo = np.quantile(latent[:, j], (1 - freq[j]) ** 2)
    hi = np.quantile(latent[:, j], 1 - freq[j] ** 2)
    col = np.zeros(n_ind, dtype=np.int8)
    col[latent[:, j] > lo] = 1
    col[latent[:, j] > hi] = 2
    geno[:, j] = col

# inject near-duplicate (rsq > 0.95) NON-causal pairs so TGS pwld reader has cnt > 0
for j in [3, spb + 3, 2 * spb + 3]:
    col = geno[:, j - 1].copy()
    flip = rng.random(n_ind) < 0.005
    col[flip] = (col[flip] + 1) % 3
    geno[:, j] = col
    freq[j] = geno[:, j].mean() / 2

# .fam .bim .bed
pref.with_suffix(".fam").write_text("\n".join(f"F{i} I{i} 0 0 0 -9" for i in range(n_ind)) + "\n")
bp = np.arange(n_snp) * 1000 + 1
pref.with_suffix(".bim").write_text(
    "\n".join(f"{chrom}\trs{j + 1}\t{bp[j] * 1e-5:.6f}\t{bp[j]}\tA\tG" for j in range(n_snp)) + "\n"
)
code_map = {0: 0b11, 1: 0b10, 2: 0b00}
with pref.with_suffix(".bed").open("wb") as f:
    f.write(bytes([0x6C, 0x1B, 0x01]))
    for j in range(n_snp):
        col = geno[:, j]
        buf = bytearray()
        for i0 in range(0, n_ind, 4):
            b = 0
            for k in range(4):
                i = i0 + k
                b |= (code_map[int(col[i])] if i < n_ind else 0) << (2 * k)
            buf.append(b)
        f.write(bytes(buf))

# block-info: header Block Chr StartBP EndBP
blk = ["Block\tChr\tStartBP\tEndBP"]
for bl in range(n_block):
    blk.append(f"{bl + 1}\t1\t{bl * spb * 1000 + 1}\t{(bl + 1) * spb * 1000}")
(outdir / "blocks.txt").write_text("\n".join(blk) + "\n")

# planted causals + consistent marginal betas b_marg = R @ b_true
R = np.corrcoef(geno.T.astype(float))
b_true = np.zeros(n_snp)
causal_idx = [10, spb + 20, 2 * spb + 30]  # rs11, rs171, rs331
b_true[causal_idx] = [0.15, 0.13, 0.14]
b_marg = R @ b_true
N = 100000
se = np.full(n_snp, 1.0 / np.sqrt(N))
b_marg = b_marg + rng.normal(0, se)
z = b_marg / se
p = np.clip([2 * (1 - 0.5 * (1 + erf(abs(zz) / sqrt(2)))) for zz in z], 1e-300, 1.0)
ma = ["SNP A1 A2 freq b se p N"] + [
    f"rs{j + 1} A G {freq[j]:.4f} {b_marg[j]:.6f} {se[j]:.6f} {p[j]:.3e} {N}" for j in range(n_snp)
]
(outdir / "toy.ma").write_text("\n".join(ma) + "\n")

# annotation: SNP Intercept Anno1  (Intercept=1 for all; Anno1 marks block 1)
an = ["SNP Intercept Anno1"] + [f"rs{j + 1} 1 {1 if j // spb == 0 else 0}" for j in range(n_snp)]
(outdir / "annot.txt").write_text("\n".join(an) + "\n")

# gene-map: EXACT 9-column header required by readGeneMapFile
gm = ["Ensgid\tGeneName\tGeneType\tChrom_hg38\tStart_hg38\tEnd_hg38\tChrom_hg19\tStart_hg19\tEnd_hg19"]
for bl in range(n_block):
    s, e = bl * spb * 1000 + 1, (bl + 1) * spb * 1000
    gm.append(f"ENSG{bl}\tGENE{bl + 1}\tprotein_coding\t1\t{s}\t{e}\t1\t{s}\t{e}")
(outdir / "genemap.txt").write_text("\n".join(gm) + "\n")

n_hi = int(((R**2 > 0.95) & (~np.eye(n_snp, dtype=bool))).sum() // 2)
print("causals:", [f"rs{i + 1}" for i in causal_idx])
print("max |off-diag r|:", float(np.nanmax(np.abs(R - np.eye(n_snp)))), "rsq>0.95 pairs:", n_hi)

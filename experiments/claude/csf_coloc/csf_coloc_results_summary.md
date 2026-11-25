# CSF protein coloc results summary

Coloc.abf of each trait GWAS against the Western et al. 2024 CSF pQTL aptamers (GRCh38).
All runs completed with no failed pairs in any trait.

| Trait | Loci | Tested | Not tested (budget) | Colocalized (PP.H4 ≥ 0.8) |
|---|---|---|---|---|
| Alzheimer's (positive control) | 54 | 130 | 399 | 33 |
| DecodeME (p < 5e-8) | 6 | 1 | 0 | 0 |
| DecodeME exploratory (p < 1e-6, MHC in) | 25 | 32 | 2 | 0 |
| DecodeME exploratory (p < 1e-5, MHC in) | 72 | 72 | 4 | 0 |
| Rheumatoid arthritis | 93 | 105 | 0 | 3 |
| LDL | 74 | 160 | 430 | 11 |
| IBD | 125 | 200 (budget reached) | 45 | 8 |
| Schizophrenia | 168 | 144 | 0 | 6 |
| Parkinson's | 26 | 65 | 12 | 16 |
| ALS | 14 | 11 | 0 | 0 |

"Loci" includes loci where no CSF protein had a nearby signal.

Including the MHC did not help: the exploratory DecodeME runs found no colocalizations even
with it in.

## DecodeME in the MHC

The only DecodeME locus there is chr6:25.7–27.7 Mb (the butyrophilin/BTN region; lead
p = 4e-9).

- 20 pairs were tested and 4 were cut by the per-locus cap.
- Every tested pair came out as PP.H3 ≥ 0.95. That includes the aptamers for genes inside the
  window: BTN3A2, BTN2A1, BTN2A2, BTN3A3 and SCGN.
- So in CSF, DecodeME's MHC signal and the local pQTLs sit on different causal variants.
- Caveat: MHC linkage disequilibrium (LD) is where coloc.abf's one-causal-variant assumption is
  weakest. SuSiE-coloc could still separate a secondary signal here, as with PRSS8 in
  Alzheimer's.

## Best non-MHC DecodeME pairs (p < 1e-5)

All are weak:

- VWA1 at chr1:1.6 Mb, cis: PP.H4 0.27.
- KCNE5: 0.20.
- TDGF1 at chr3:49 Mb: 0.14.

## Sources

Per-trait summaries come from `coloc_run_summary.py`; its output was saved as
`coloc_summary_*.log` and `coloc_top_*.log`, which are not committed.

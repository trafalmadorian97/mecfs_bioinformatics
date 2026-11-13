# Trait x UKB-PPP colocalization, v1 run summary (2026-10-29)

coloc.abf, EUR discovery cohort, HapMap3 hit table (p < 5e-8, 250 kb merge gap), trait windows
+-500 kb (merged), pair budget 30, at most 5 pairs per locus, MHC excluded. Trait inputs
harmonized to GRCh38 against UCSC hg38 + the 1000 Genomes EUR 30x panel (gnomAD v4.1 is not
built locally). Logs: experiments/claude/logs/*_ppp_coloc_rerun.log (uncommitted).

| Trait | Harmonization trust | Loci | Candidate pairs | Tested | Failed | No candidates | MHC | PP.H4 > 0.8 |
|---|---|---|---|---|---|---|---|---|
| Willer LDL | untrusted (original alleles) | 74 | 3,924 | 30 | 0 | 0 | 2 | PCSK9 (cis, 1.00), CELSR2 (cis, 0.9997) |
| DecodeME GWAS-1 | trusted | 6 | 11 | 8 | 0 | 3 | 1 | none (max RABGAP1L 0.076) |
| deCODE seropositive RA (MAF >= 0.01) | untrusted | 91 | 1,655 | 30 | 0 | 14 | 1 | CCL19 (trans, 0.9987) at chr1:113.7 Mb (PTPN22 locus) |
| Liu 2023 IBD (EUR) | trusted | 125 | 3,709 | 30 | 0 | 2 | 2 | MST1 (cis, 0.986) at chr3:49.7 Mb |

Positive control: PCSK9 cis PP.H4 = 1.00 with the same lead variant in both panels
(chr1:55,039,974, R46L). At large merged loci with several trait signals (APOE, LDLR region)
coloc.abf's single-causal-variant assumption pushes pairs to H3; SuSiE-coloc is the follow-up.

## R coloc cross-check (2026-11-13)

`experiments/claude/protein_coloc/r_coloc_crosscheck.py` reran every tested pair of
`decode_me_gwas_1_exploratory_2_ppp_coloc` (517 pairs) through R coloc 5.2.3 `coloc.abf`
(throwaway `pixi exec` env) on the exact saved windows and priors. Result: PASS.
- R nsnps == N_ALIGNED for all 517 pairs.
- max |dPP| over all pairs and hypotheses: 2.8e-15; over the 8 hits (PP.H4 > 0.8): 8.6e-16.
- R's prior SDs confirmed from source: cc with beta/varbeta uses 0.2 (s and MAF unused);
  quant with sdY = 1 uses 0.15. Same as ColocAbf.
This verifies the ABF/posterior arithmetic only, not the harmonization/alignment upstream.

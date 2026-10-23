# NET (SLC6A2) A457P — what has been learned since Shannon et al. 2000

Variant: SLC6A2 c.1369G>C, p.Ala457Pro (NM_001172501.3, MANE; same numbering on NM_001043),
rs121918126, GRCh38 16-55698005-G-C, GRCh37 16-55731917-G-C. Codon GCC -> CCC, in TM9.

Reproduce: `pixi r python experiments/claude/net_a457p_gnomad/query_gnomad.py > experiments/claude/net_a457p_gnomad/query_gnomad.json`
then `pixi r python experiments/claude/net_a457p_gnomad/summarize.py` (queried 2026-10-23).

## gnomAD

- A457P is absent from gnomAD v4.1 (exomes AN ~1.46M, genomes AN ~152k at this site) and from v2.1.
  The site is well covered: other alleles at the same position are called with full AN.
  Upper 95% bound on allele frequency is roughly 3 / 1.61M ~ 2e-6.
- Same codon, other alleles (v4 joint): A457T (G>A) AC=6 (5 NFE, 1 AFR, no homozygotes),
  A457S (G>T) AC=2 (EAS). So Ala457 tolerates rare conservative substitutions; Pro is never seen.
- Gene constraint, v4 (MANE ENST00000568943): pLI 0.96, LOF o/e 0.36 (LOEUF 0.50),
  missense o/e 0.72 (90% CI 0.68-0.77), mis_z 3.35. v2: LOEUF 0.48, mis o/e 0.70.
  Moderately intolerant to LoF; mildly missense-constrained.

## In silico scores for A457P (dbNSFP via myvariant.info; gnomAD only shows predictors for observed variants)

| Predictor | Score | Call |
|---|---|---|
| CADD v1.x PHRED | 28.1 | deleterious |
| REVEL | 0.788 | above 0.773 PP3-moderate threshold (Pejaver 2022) |
| AlphaMissense | 0.945 (MANE) | likely pathogenic |
| MetaRNN | 0.978 | D |
| ESM1b | -13.4 | D |
| VARITY_R | 0.956 | D |
| ClinPred | 0.990 | D |
| BayesDel (no AF) | 0.250 | D |
| MVP | 0.913 | |
| PolyPhen2 HDIV/HVAR | 0.972 / 0.924 | probably damaging |
| SIFT | 0.004 | deleterious |
| PrimateAI | 0.771 | T (the one dissenter) |
| phyloP 470-way mammal / 100-way vert | 11.76 / 5.17 | highly conserved |

For comparison A457T: CADD 23.8, PolyPhen 0.21; A457S: CADD 17.8, SIFT 1.0. SpliceAI/Pangolin ~0 at this position.

## ClinVar

VCV000014006: Pathogenic, "no assertion criteria provided", single submission (OMIM, last evaluated 2002-03-01),
condition orthostatic intolerance (OMIM 604715). No clinical lab submissions since.

## Literature since 2000 (PubMed: A457P/Ala457Pro/rs121918126 -> 11 hits)

Functional:
- Paczkowski 2002 Pharmacogenetics (PMID 11875370): COS-7; surface Bmax ~40% WT, but 50-fold higher Km for NE;
  argued loss of function is mainly reduced substrate affinity, not only trafficking. Altered cocaine/nisoxetine affinities.
- Hahn, Robertson & Blakely 2003 J Neurosci (PMID 12805287): impaired glycosylation/maturation, surface ~30% WT;
  A457P co-immunoprecipitates with WT NET and is dominant negative on WT uptake and surface expression,
  which explains a strong phenotype in heterozygotes.
- Hahn 2005 Mol Pharmacol (PMID 15894713): characterized 10 other NET coding SNPs; A369P and N292T also
  trafficking-defective / dominant negative. A457P used as reference null.
- Shirey-Rice 2013 Dis Model Mech (PMID 23580201): A457P knock-in mice. Heterozygotes: reduced brain/sympathetic NE uptake,
  raised plasma/urine NE, low DHPG:NE, tachycardia without BP or baroreflex change, plus anxiety-/cognition-related
  behaviours; homozygotes near-zero uptake. In vivo confirmation of dominant-negative, POTS-like phenotype.
- Brown, Shetty & Foster 2025 J Neurochem (PMID 39395208): NET is palmitoylated (DHHC1 raises WT NET);
  co-expressing DHHC1 did not rescue A457P but increased N-terminal fragments; C44 loss mimics A457P low expression.
- Avsar 2025 Nucleosides Nucleotides Nucleic Acids (PMID 39761170): bioinformatic survey of 504 NET missense SNPs;
  A457P among 9 predicted high-risk.

Population / clinical:
- Ivancsits 2003 Am J Med Sci (PMID 12589229): 0/14 young men with full OI syndrome carry A457P (exon 9 sequencing).
- Nickander 2005 Auton Neurosci (PMID 15946904, Mayo): 0/29 POTS patients carry A457P.
- Schroeder & Jordan 2012 AJP-Heart review (PMID 23023867): no additional A457P families identified by them or others.
- Bayles 2012 ATVB and Khan 2017 JCI Insight: POTS patients *without* A457P have reduced NET protein,
  attributed to epigenetic silencing (chromatin modification; let-7i/MeCP2), not coding variants.
- Qu 2025 Clin Auton Res (PMID 39964606): pediatric POTS GWAS (207 cases) + WES (87 cases); no genome-wide
  significant loci; A457P cited only as prior work.

## Bottom line

- Effect on the protein is well established: strongly damaging by every modern predictor except PrimateAI,
  and by multiple independent functional assays plus a knock-in mouse.
- Effect on disease rests on one family. It is ultra-rare (absent from ~807k gnomAD v4 individuals) and has not
  been found again in POTS/OI screens, so it is a private familial variant and explains essentially none of POTS
  in the population. The "NET deficiency" mechanism generalized instead through acquired/epigenetic reduction of
  NET expression.
- ClinVar has not been reassessed under ACMG; a modern curation would likely use PM2 (absent), PS3 (functional/mouse),
  PP3 (REVEL 0.788), and PP1 (single-family segregation). That probably lands at Likely Pathogenic for an autosomal-dominant
  OI/NET-deficiency phenotype, limited mainly by having only one family.

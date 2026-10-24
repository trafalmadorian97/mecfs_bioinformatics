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

Full text read (user-supplied PDF): Bayles et al. 2012 ATVB. Genetics arm was 20 POTS vs 65 controls: none of 16 known
NET coding variants (incl. A457P) in the POTS patients; synonymous rs5564 (exon 5 splice consensus) nominally associated
(P=0.0014), unreplicated. Epigenetic arm: leukocyte ChIP in only 6 vs 6 women (more H3K9me3/H3K27me3/MeCP2, less
H3K9/14ac/H3K4me3 at the SLC6A2 promoter in POTS); leukocyte NET protein 0.66x controls by Western blot; NET mRNA
in leukocytes was too low to quantify, which undercuts leukocytes as a proxy for sympathetic neurons.

## PheWAS-style survey (phewas_slc6a2.py, phewas_lookups.py; summaries via summarize_phewas.py, summarize_lookups.py)

Queried 2026-10-24. Outputs: phewas_slc6a2.json (committed) and genebass_slc6a2.parquet, finngen_slc6a2_lookups.parquet,
ot_gwas_credible_sets_containing_lookups.parquet (gitignored; regenerate with phewas_lookups.py).

### a) Variants that change the protein

- Genebass (UKB 394k exomes, SAIGE-GENE SKAT-O, ~4,500 phenotypes): no Bonferroni-significant phenotype for
  pLoF (16 variants; top: M99 biomechanical lesions p=2.7e-5, plateletcrit 1.6e-4, essential hypertension 2.2e-4 with
  carriers at higher risk, "trouble concentrating" 2.4e-4, "irritability" 2.1e-4) or missense|LC (280 variants; top
  MCV 1.9e-4). Counts of p<1e-3 (9 and 12 vs ~4.5 expected) are similar to the synonymous negative control (20),
  so no excess. Open Targets gene_burden (Genebass / AZ PheWAS / Regeneron) has no SLC6A2 entries, consistent.
- Common coding variants (gnomAD v4 AF >= 0.1%): T99I rs1805065 (1.7%), V245I rs1805066 (0.9%), C27Y rs11568343 (0.8%),
  start-lost p.Met1? rs552629380 (0.13%), T283R rs45564432 (0.1%, absent from FinnGen). FinnGen R12 PheWAS
  (2,470 endpoints, Bonferroni 2e-5): nothing significant; all top hits p >= 1.5e-4 and look like noise.
- MVP PheWAS (Verma 2024, via Open Targets credible sets L2G -> SLC6A2): a cluster of very rare (EAF ~0.04-0.1%) imputed
  variants with huge effects (log-OR 2-3.7) for unrelated phecodes: hyperaldosteronism/adrenal hyperfunction
  (rs564420311, p=2.8e-12), spinal cord injury, anaphylaxis, amputation stump, chondromalacia. The pattern (rare,
  imputed, implausible ORs, unrelated traits, PICS sets) looks artefactual; hyperaldosteronism is biologically curious but
  unconfirmed.

### b) Variants that change expression

- SLC6A2 cis-eQTLs exist (GTEx v8/v10, TwinsUK): strongest in skin (rs17841327 / rs9937434, alt lowers expression,
  p~1e-71..1e-94), adrenal gland (rs9937434, NES -0.45, p=6e-18), testis (rs187714), esophagus (rs55830311), plus
  sQTL/tuQTL. Caveat: the tissues where NET matters (sympathetic ganglia, locus coeruleus, adrenal medulla) are absent
  or diluted in GTEx (adrenal samples are mostly cortex), so these eQTLs are uncertain proxies for neuronal NET.
- Colocalisation with GWAS (Open Targets, COLOC/eCAVIAR, h4 >= 0.8):
  - testis eQTL (rs187714) <-> birth weight (Warrington 2019, lead p=3e-10; UKB p~2e-6), h4 0.99, L2G top gene SLC6A2.
  - testis/placenta tuQTL (rs42460, rs36007) <-> placental weight (Beaumont, p=7.6e-10; L2G top gene CES1).
  - esophagus eQTL (rs55830311) <-> height (Yengo), appendicular lean mass, plasma MMP-2 levels (h4 0.83-0.93; MMP2 is
    a neighbouring gene, so likely a shared regulatory region rather than NET biology).
  - The adrenal and skin eQTLs do not colocalise with any GWAS trait.
  None of these is an autonomic, cardiovascular or psychiatric trait.
- FinnGen R12 PheWAS of the eQTL leads: nothing significant. Sub-threshold: expression-lowering alleles of the
  skin/adrenal eQTL associate with otosclerosis (beta +0.11-0.12, p 5e-5 / 1e-4) and, weakly, with less essential
  hypertension (p~0.01); rs42460 with "other menopausal disorders" (p=4.5e-6, passes per-variant Bonferroni but not
  across the 11 variants looked up).
- rs3785151 (MVP sensorineural hearing loss, p=3.8e-9, PIP 0.63, L2G SLC6A2 0.83; only a weak skin eQTL) does not
  replicate in FinnGen hearing-loss endpoints (p>0.05).
- FinnGen note: its variant API returns statistics only for nominally associated endpoints (~5%) and nulls otherwise;
  Bonferroni/expected counts above use all 2,470 endpoints. FinnGen has no POTS/orthostatic intolerance endpoint, and
  ICD-coded biobanks capture POTS poorly in general.
- Mouse (MGI/IMPC via Open Targets): Slc6a2 knockouts show raised plasma noradrenaline/adrenaline and dopamine, raised
  heart rate (A457P knock-in het), lower body temperature and weight, anxiety/despair-like behaviour, enhanced
  cocaine/psychostimulant responses, increased antinociception, and in IMPC homozygotes small heart, abnormal gait and
  incomplete preweaning lethality.

### Context: other ways NET function has been changed in humans

- Pharmacological NET blockade reproduces the POTS phenotype: reboxetine in healthy volunteers (Schroeder 2002
  Circulation) raised tilt heart rate to 119 vs 84 bpm on placebo, with fewer vasovagal reactions.
- Common NET promoter/synonymous variants (rs28386840 -3081A/T, rs2242446 -182T/C, rs5569 G1287A) have a long, inconsistent
  candidate-gene literature for ADHD risk; a 2021 meta-analysis found rs28386840 T carriers responded better to
  methylphenidate (OR ~2), still needing replication. Open Targets has no ADHD (or any psychiatric) GWAS credible set
assigned to SLC6A2.

### Interpretation

Human population genetics has little to say about NET dosage beyond the A457P family: protein-altering variation is
too rare for gene-level tests to be well powered (no signal in UKB exomes), and the eQTLs that exist are in tissues of
doubtful relevance and colocalise only with growth traits (birth weight, placental weight, height). Most of the
phenotype evidence for reduced NET comes from pharmacology (NET inhibitors) and mouse models. Possible next steps:
look up the SLC6A2 region in DecodeME (ME/CFS) and POTS-adjacent GWAS available in this repo, and test
quantitative autonomic traits (resting/standing heart rate, BP) in Pan-UKBB summary statistics, which ICD-based
FinnGen endpoints cannot capture.

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

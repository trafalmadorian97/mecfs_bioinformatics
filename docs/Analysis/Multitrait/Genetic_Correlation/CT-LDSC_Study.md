---
hide:
- navigation
- toc
---
# CT-LDSC Study

I applied [cross-trait linkage disequilibrium score regression](../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md)[@bulik2015atlas] to estimate [genetic correlation](../../../Bioinformatics_Concepts/Genetic_Correlation.md) between a diverse phenotypes, including the ME/CFS phenotype defined in DecodeME[@genetics2025initial].  The genetic correlations are plotted in the heatmap below. Cells with asterixes are Bonferroni-corrected statistically significant.

<iframe src="../../../../_figs/initial_genetic_correlation_by_ct_ldsc_plot.html" style="width:100%; height:750px; border:none;"></iframe>

The traits are:

- **Asthma** from the GWAS of Han et al.[@han2020genome]  Several theories relating ME/CFS to the activation of mast cells (_todo: add citation_), but the evidence for these theories has is often anecdotal. I added asthma as a prototypical example of a common disease involving the IgE/ mast cell immune axis (See _Janeway's Immunobiology,_[@murphy2022janeway] _Chapter 14: Allergic Diseases_).
- Diastolic blood pressure (**DBP**) from the Keaton et al.[@keaton2024genome] Several authors theorize that ME/CFS is related to dysregulation of blood flow.  I included the Keaton study as a well-powered, well-measured GWAS of a blood-flow related trait.
- **Educational attainment** from Lee et al[@lee2018gene], an example of an extremely complex trait with many determinants, some which lie in the central nervous system.
- ME/CFS, from **DecodeME**[@genetics2025initial].
- **Multi-site pain**, from Johnston et al.[@johnston2019genome] The DecodeME preprint reported colocalization between the ME/CFS signal and a multi-site chronic pain signal at a locus on chromosome 17.  I added multi-site pain to investigate whether there was evidence of genome-wide matching between pain and ME/CFS, in addition to localized matching.
- **Inflammatory bowel disease** from Liu et al.[@liu2023genetic] Like asthma, IBD is an immunological disease, but it operates via a different subsystem of the immune system.  I was interested in observing how this contrast is reflected in genetic correlations.
- **Schizophrenia** from the 2022 PGC study.[@trubetskoy2022mapping] Unlike certain other common complex diseases of the central nervous system, schizophrenia has a severe and relatively distinctive phenotype, which reduces the likelihood of diagnostic error.  In addition to involving genes active in neurons, genetic risk for schizophrenia appears to also be driven by the immunological genes ( see _"Identification of Genes for Schizophrenia Highlights the Interplay of Rare and Common Risk Variants_" in Chapter 1 of Kandel et al.[@kandel2021principles] and also the [S-LDSC analysis of schizophrenia](../../PGC_2022_(SCH)/S-LDSC_SCH_Analysis.md) in this repo).  It may thus be reasonable to describe schizophrenia as a neuroimmune condition. ME/CFS has also been called a neuroimmune condition, though the phenotypes of the two diseases are vastly different.

## Comment on results

- The trait most strongly genetically correlated with ME/CFS is multi-site pain. This is consistent with the [observation](https://www.s4me.info/threads/genome-wide-association-study-of-multisite-chronic-pain-in-uk-biobank-2019-johnston.44973/post-674317) that the Manhattan plots for DecodeME and multi-site pain appear to match closely at several loci.  This intriguing finding  deserves follow-up.
- There are two other significant genetic correlations: with schizophrenia and with asthma.  As is often the case in genetic correlation studies, this finding is difficult to interpret. There are numerous possible theories.  For instance, the correlation between asthma and ME/CFS could reflect an IgE-related immune etiology, while the correlation between ME/CFS and Schizophrenia could reflect neurological etiology; or both correlations could reflect some common immune pathway. 

## Next steps

Potential next steps:

- Techniques like GenomicSEM[@grotzinger2019genomic], LCV[@o2018distinguishing], and MiXeR may elucidate the causal structure of correlations discovered so far.
-  The inclusion of additional well-chosen phenotypes in the correlation study may shed light on existing traits through triangulation.
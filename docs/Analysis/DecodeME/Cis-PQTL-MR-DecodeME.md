# Cis pQTL MR Analysis


I applied [Mendelian Randomization](../../Bioinformatics_Concepts/Mendelian_Randomization.md) to DecodeME GWAS 1[@genetics2025initial] using cis-pQTLs from the [UK Biobank Pharma Proteomics Project](../../Data_Sources/UKBB_PPP.md)[@sun2023plasma] as instruments.  I aimed to identify candidate proteins that may be causal in the DecodeME disease process.


Since each protein in the UKBB PPP has at most one cis-pQTL, I used the Wald Ratio method to estimate the causal effect of the proteins on ME/CFS.

I used the R package `TwoSampleMR`.


NOTE: since the controls for the DecodeME project come from the UK biobank, and since the pQTLs from the UKBB PPP are also derived from analysis of the UK biobank, we do not actually have two independent samples, so we violate some of the assumptions of `TwoSampleMR`.  Thus it is possible that p values generated from `TwoSampleMR` may be artificially low.

## Results


There was one protein for which the Wald Ratio method as implemented by `TwoSampleMR` produced a significant result: RABGAP1L.  This


| Assay Target   |        b |        pval | Function [CC]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
|:---------------|---------:|------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| RABGAP1L       | -1.01937 | 8.16913e-07 | FUNCTION: GTP-hydrolysis activating protein (GAP) for small GTPase RAB22A, converting active RAB22A-GTP to the inactive form RAB22A-GDP (PubMed:16923123). Plays a role in endocytosis and intracellular protein transport. Recruited by ANK2 to phosphatidylinositol 3-phosphate (PI3P)-positive early endosomes, where it inactivates RAB22A, and promotes polarized trafficking to the leading edge of the migrating cells. Part of the ANK2/RABGAP1L complex which is required for the polarized recycling of fibronectin receptor ITGA5 ITGB1 to the plasma membrane that enables continuous directional cell migration (By similarity). {ECO:0000250|UniProtKB:A6H6A9, ECO:0000269|PubMed:16923123}. |⏎


The Candidate Genes supplement to the DecodeME preprint notes that "RABGAP1L inhibits the entry and/or promotes the expulsion of bacteria or viruses. It is a viral restriction factor. Reduction in RABGAP1L expression would
be expected to enhance susceptibility to bacterial and viral infection, which often precedes
initial ME/CFS symptoms."
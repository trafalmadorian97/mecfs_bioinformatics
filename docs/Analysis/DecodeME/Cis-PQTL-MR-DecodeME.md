---
hide:
- toc
---
# cis-pQTL MR Analysis


I applied [Mendelian Randomization](../../Bioinformatics_Concepts/Mendelian_Randomization.md) to DecodeME GWAS-1[@genetics2025initial] using cis-pQTL from the [UK Biobank Pharma Proteomics Project](../../Data_Sources/UKBB_PPP.md)[@sun2023plasma] as instruments.  I aimed to identify candidate proteins that may be causal in the ME/CFS disease process.


Since each protein in the UKBB PPP has at most one cis-pQTL, I used the Wald Ratio method to estimate the causal effect of the proteins on ME/CFS.

I used the R package `TwoSampleMR`.


**NOTE**: since the controls for the DecodeME project come from the UK Biobank, and since the pQTLs from the UKBB PPP are also derived from analysis of the UK Biobank, we do not actually have two fully independent samples, so we partially violate the assumptions of `TwoSampleMR`.  This could introduce some error. 

## Results

### RABGAP1L

At a significance level of 0.01 (Bonferroni-corrected), `TwoSampleMR` identified one potentially causal plasma protein: RABGAP1L.


| Assay Target   |        b |        pval | UniProt Function [CC]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
|:---------------|---------:|------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| RABGAP1L       | -1.01937 | 8.16913e-07 | FUNCTION: GTP-hydrolysis activating protein (GAP) for small GTPase RAB22A, converting active RAB22A-GTP to the inactive form RAB22A-GDP (PubMed:16923123). Plays a role in endocytosis and intracellular protein transport. Recruited by ANK2 to phosphatidylinositol 3-phosphate (PI3P)-positive early endosomes, where it inactivates RAB22A, and promotes polarized trafficking to the leading edge of the migrating cells. Part of the ANK2/RABGAP1L complex which is required for the polarized recycling of fibronectin receptor ITGA5 ITGB1 to the plasma membrane that enables continuous directional cell migration (By similarity). \{ECO:0000250 \|UniProtKB:A6H6A9, ECO:0000269 \|PubMed:16923123}. |⏎


The mechanism by which RABGAP1L could affect ME/CFS is not immediately clear.  The[ Candidate Genes supplement to the DecodeME preprint](https://www.pure.ed.ac.uk/ws/portalfiles/portal/533352484/Candidate_Genes.pdf) notes that "RABGAP1L inhibits the entry and/or promotes the expulsion of bacteria or viruses. It is a viral restriction factor. Reduction in RABGAP1L expression would
be expected to enhance susceptibility to bacterial and viral infection, which often precedes
initial ME/CFS symptoms."


To further investigate, I generate two region plots.  The first shows the DECODE ME GWAS signal in region of the RABGAP1L gene:

![mecfs_rabgap1l_region_plot](https://github.com/user-attachments/assets/e2a71bbf-e803-4b57-9a7d-a22b4f9fd2f2)


The second shows the same region from the UK Biobank Pharma Proteomics Project GWAS of plasma RABGAP1L levels:

![rabgap1l_region_plot](https://github.com/user-attachments/assets/161f1b0e-e309-4c0e-a88a-d312b8dc80ee)

From visual inspection, it seems clear that at least the primary signals do not colocalize: the GWAS peaks are clearly in different locations.  However, it could be that there are secondary signals that do colocalize.  A colocalization algorithm may clarify the question.

### BTN1A1

If we relax the significance threshold to 0.05, one additional potentially causal protein is identified: BTN1A1

| Assay Target   |         b |        pval | UniProt Function [CC]                                                                                                                                                                                                                                                                                                                                                          |
|:---------------|----------:|------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| BTN1A1         |  0.284424 | 2.38381e-05 | FUNCTION: May function in the secretion of milk-fat droplets. May act as a specific membrane-associated receptor for the association of cytoplasmic droplets with the apical plasma membrane (By similarity). Inhibits the proliferation of CD4 and CD8 T-cells activated by anti-CD3 antibodies, T-cell metabolism and IL2 and IFNG secretion (By similarity). {ECO:0000250}. |⏎  


This protein appears to be anti-inflammatory, and higher levels are associated with a greater risk of ME/CFS.


Again, to further investigate, I compared region plots in the vicinity of BTN1A1 from DecodeME and from the UK Biobank Pharma Proteomics Project GWAS of plasma BTN1A1 levels.

The first plot is from DecodeME:

![decode-me-btn1a1-region](https://github.com/user-attachments/assets/71e71dde-84d9-4e8e-8d2d-a5e28e075914)


The second is from UK Biobank Pharma Proteomics Project:


![ukbb_ppp_btn1a1](https://github.com/user-attachments/assets/400039dc-0e31-4cfd-b219-c3cea2350a3c)


In this comparison, it seems even more obvious that the primary signals do not colocalize: they are separated by a recombination boundary.


## Follow-up questions

- Does co-localization analysis suggest that the DecodeME causal variant and the cis-pQTL co-localize?

- Do the cis-pQTLs identified above co-localize with trans-pQTLs for other proteins? If so, this could suggest causal protein regulatory chains.
# PPP Coloc
For further insights into the molecular mechanisms of multisite pain,, I ran [ABF Colocalization](../../../Bioinformatics_Concepts/Colocalization.md)[@giambartolomei2014bayesian] to search for loci at which the Johnston et al. GWAS of multisite pain[@johnston2019genome] colocalized with proteomics GWAS from the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md)[@sun2023plasma]

Locus/protein pairs with greater than 80 percent probability of colocalization are shown below:


{{
coloc_results_table("docs/_figs/johnston_et_al_pain_ppp_coloc_full_ppp_coloc_display.parquet","pain-ppp-coloc")
}}



Multisite pain has cis colocalizations with SUOX, CDH6, ECM1, SORT1, PGLYRP4, and C8B. Interestingly, the cis colocalization SUOX is [shared with fibromyalgia](../../Fibromyalgia/Kerrebijin_et_al/PPP_COLOC_Kerrebijin_et_al.md).  [CDH6](https://www.genecards.org/card/CDH6) appears to have a role in cell adhesion.


The results above should be taken with a grain of salt. Both the multisite pain and the PPP proteomic GWAS source participants from the [UK Biobank](../../../Data_Sources/UKBB.md).  This violates an assumption used in the derivation of ABF colocalization.

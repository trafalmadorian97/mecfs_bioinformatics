# PPP CT-LDSC

I applied [Cross Trait Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md) (CT-LDSC)[@bulik2015atlas] to estimate [genetic correlation](../../../Bioinformatics_Concepts/Genetic_Correlation.md) between the Han et al.[@han2020genome] GWAS of Asthma and the Olink proteomic [@wik2021proximity] GWAS of the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma].


## Results

As is standard for LDSC analysis, I restricted the summary statistics to Hapmap 3 variants, and excluded the MHC region.  I used the standard thousand genomes linkage disequilibrium scores provided by the authors of LDSC. To focus on trans effects, I excluded the cis regions from the proteomic GWAS.



The results are below:

{{ ppp_rg_data_table("docs/_figs/han_asthma_ppp_rg_cis_excluded_display_frame.parquet", id="ukbb-ppp-ldsc-han-asthma-rg" )}}



## Interpretation



- The top most-significant protein is [CLC (Charcot-Leyden Crystal Galectin)](https://www.genecards.org/card/CLC).  Given that CLC is known to be expressed on eosinophils, and eosinophils are known to be important to allergy and asthma, one can make argument for the plausibility of the genetic correlation between asthma and plasma CLC.

# PPP CT-LDSC

I applied [Cross Trait Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md) (CT-LDSC)[@bulik2015atlas] to estimate [genetic correlation](../../../Bioinformatics_Concepts/Genetic_Correlation.md) between the Johnston et al. GWAS of multisite pain[@johnston2019genome] and the Olink proteomic [@wik2021proximity] GWAS from the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma]


## Results

As is standard for LDSC analysis, I restricted the summary statistics to Hapmap 3 variants, and excluded the MHC region.  I used the standard thousand genomes linkage disequilibrium scores provided by the authors of LDSC. To focus on trans effects, I excluded the cis regions from the proteomic GWAS.


The results are below:

{{ ppp_rg_data_table("docs/_figs/johnston_et_al_pain_ppp_rg_cis_excluded_display_frame.parquet", id="ukbb-ppp-ldsc-pgc-pain-rg" )}}


## Interpretation

Interestingly, there are a very large number of plasma proteins that are genetically correlated with multi-site pain.  Interesting examples:

- LEP (Leptin), strongly involved in the control of body weight and appetite[^leptin_note].
- CFH (Complement Factor H), which is important to the innate immune system.
- GDF15 (Growth Differentiating factor 15), which was [also genetically correlated with myocardial infarction](../../Myocardial_Infarction/Verma_et_al/1_MV_MI_ppp_rg.md).  


It is important that we interpret these results carefully. They tell use that there is a correlation between genetic factors that predict multisite pain, and genetic factors that predict the blood proteins above.  They do not tell us the causal relationship between multisite pain and these blood proteins.  It will be interesting to see if we can learn more by applying causal inference techniques.


[^leptin_note]:  See Yeo 2018[@yeo2018gene] for a readable popular science account of the genetics of eating, which includes coverage of the discovery of leptin.
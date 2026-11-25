# PPP CT-LDSC

I applied [Cross Trait Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md) (CT-LDSC)[@bulik2015atlas] to estimate [genetic correlation](../../../Bioinformatics_Concepts/Genetic_Correlation.md) between the Kerrebejin et al. GWAS of Fibromyalgia[@kerrebijn2026genetic] and the Olink proteomic [@wik2021proximity] GWAS from the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma].



## Results

As is standard for LDSC analysis, I restricted the summary statistics to Hapmap[@international2005haplotype] 3 variants, and excluded the MHC region.  I used the standard thousand genomes linkage disequilibrium scores provided by the authors of LDSC. To focus on trans effects, I excluded the cis regions from the proteomic GWAS.


The results are below:


{{ ppp_rg_data_table("docs/_figs/kerrebeijin_et_al_pain_ppp_rg_cis_excluded_display_frame.parquet", id="ukbb-ppp-ldsc-pgc-fibro-rg" )}}

There are large number of proteins whose plasma levels are genetically correlated with fibromyalgia  It is interesting to observe a number of inflammation-related proteins with significant genetic correlations, such as


- **RARRES2**: ["Acts as a chemotactic factor for leukocyte populations expressing CMKLR1"](https://www.uniprot.org/uniprotkb/Q29RS5/entry)
- **IL1RN**:  Interleukin-1 receptor antagonist.  Important to the regulation of the immune response.
- **CFH**: [Modulates compliment activation.](https://www.uniprot.org/uniprotkb/P08603/entry)

Some significant proteins have neurological functions:

- **IGSF9**: [Involved in regulation of synapse maturation](https://www.genecards.org/card/IGSF9).
- **NTRK3**: ["neurotrophic receptor tyrosine kinase 3, a membrane-bound receptor tyrosine kinase that transduces neurotrophin signals""](https://www.genecards.org/card/NTRK3)

While other significant proteins have very diverse/unclear functions, so that it is difficult to associate them with a particular pathway.

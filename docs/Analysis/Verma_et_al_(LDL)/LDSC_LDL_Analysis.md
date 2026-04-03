
# LDSC Analysis

I applied [Linkage Disequilibrium Score Regression](../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] to the [Million Veterans](../../Data_Sources/Million_Veterans_Program.md) LDL GWAS of Verma et al.[@bulik2015ld] to estimate [heritability](../../Bioinformatics_Concepts/Heritability.md) and look for signs of population stratification or confounding.

The results are shown in the table below:

--8<--  "docs/_figs/million_veterans_ldl_ldsc_heritability_markdown.mdx"


The heritability of 0.1145 is low-to-moderate.  The LD score regression intercept of 1.586 is large, and indicates possible confounding or population stratification. This finding may be grounds to take any results that build on the Verma et al. summary statistics with a grain of salt.  At minimum, we should compare with other LDL GWAS.

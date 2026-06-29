# H-MAGMA
I applied H-MAGMA[@sey2020computational] to the [GWAS-by-subtraction](../../../../../Bioinformatics_Concepts/GWAS_By_Subtraction.md)[@demange2021investigating; @huang2024gwas] residual of [DecodeME](../../../../../Data_Sources/DecodeME.md)[@genetics2025initial] minus Johnston et al.'s GWAS of multisite pain[@johnston2019genome].


## Results

I used the 6 standard variant-to-gene assignment maps provided by the authors of H-MAGMA. Like in the other MAGMA analyses, I used the European 1000-genomes linkage disequilibrium reference. The Manhattan plots below illustrate the gene-level MAGMA results.

### Adult Brain


{{ plotly_embed("docs/_figs/decode_me_minus_johnston_ols_h_magma_adult_brain_gene_manhattan_plot.html", id="h-magma-adult-brain-decode-me-minus-johnston") }}


### Fetal Brain


{{ plotly_embed("docs/_figs/decode_me_minus_johnston_ols_h_magma_fetal_brain_gene_manhattan_plot.html", id="h-magma-fetal-brain-decode-me-minus-johnston") }}


### Cortical Neuron



{{ plotly_embed("docs/_figs/decode_me_minus_johnston_ols_h_magma_cortical_neuron_gene_manhattan_plot.html", id="h-magma-cortical-neuron-decode-me-minus-johnston") }}



### Midbrain


{{ plotly_embed("docs/_figs/decode_me_minus_johnston_ols_h_magma_midbrain_da_gene_manhattan_plot.html", id="h-magma-midbrain-decode-me-minus-johnston") }}



### IPSC-derived Astrocyte


{{ plotly_embed("docs/_figs/decode_me_minus_johnston_ols_h_magma_ipsc_derived_astro_gene_manhattan_plot.html", id="h-magma-astro-decode-me-minus-johnston") }}



### IPSC-derived Neuron




{{ plotly_embed("docs/_figs/decode_me_minus_johnston_ols_h_magma_ipsc_derived_neuro_gene_manhattan_plot.html", id="h-magma-ipc-neuron-decode-me-minus-johnston") }}


## Discussion

It is interesting to observe that the BTN1A1/BTN2A2 region of chromosome 6 is significant in many of the above results. This suggests that this region may be characteristic of ME/CFS and distinct from general multisite pain.
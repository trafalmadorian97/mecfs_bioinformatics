# H-MAGMA


I applied H-MAGMA[@sey2020computational] to the [DecodeME](../../../Data_Sources/DecodeME.md)[@genetics2025initial] GWAS-1 summary statistics.  H-MAGMA operates identically to standard [MAGMA](../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma], except that variants are assigned to genes not using proximity, but according to the results of [Hi-C](../../../Bioinformatics_Concepts/Epigenetics.md#hi-c) chromatin interaction experiments performed on a variety of neural and glial cells.



## Results
I used the 6 standard variant-to-gene assignment maps provided by the authors of H-MAGMA. Like in the other MAGMA analyses, I used the European 1000-genomes linkage disequilibrium reference. The Manhattan plots below illustrate the gene-level MAGMA results.



### Adult Brain



{{ plotly_embed("docs/_figs/decode_me_h_magma_adult_brain_gene_manhattan_plot.html", id="h-magma-adult-brain-decode-me-gene-manhattan") }}


### Fetal Brain


{{ plotly_embed("docs/_figs/decode_me_h_magma_fetal_brain_gene_manhattan_plot.html", id="h-magma-fetal-brain-decode-me-gene-manhattan") }}


### Cortical Neuron 


{{ plotly_embed("docs/_figs/decode_me_h_magma_cortical_neuron_gene_manhattan_plot.html", id="h-magma-cortical-neuron-decode-me-gene-manhattan") }}




### Midbrain


{{ plotly_embed("docs/_figs/decode_me_h_magma_midbrain_da_gene_manhattan_plot.html", id="h-magma-midbrain-decode-me-gene-manhattan") }}



### IPSC-derived Astrocyte


{{ plotly_embed("docs/_figs/decode_me_h_magma_ipsc_derived_astro_gene_manhattan_plot.html", id="h-magma-ipsc-astrocyte-decode-me-gene-manhattan") }}



### IPSC-derived Neuron


{{ plotly_embed("docs/_figs/decode_me_h_magma_ipsc_derived_neuro_gene_manhattan_plot.html", id="h-magma-ipsc-neuron-decode-me-gene-manhattan") }}



## Discussion

- As would be expected, switching from MAGMA to H-MAGMA preserved the broad pattern of significance across the genome.
- On the other hand, some of the individual significant genes change.  For instance, BTN1A1 is significant in a number of the H-MAGMA analyses, but was not significant in the original MAGMA analysis  [^footnote1].



[^footnote1]: The authors of H-MAGMA included non-protein-coding RNA transcripts in their annotation file, in addition to protein coding genes. This results in around 54k gene-like entities, and a Bonferoni threshold of $0.05/54000 \approx 9 \times 10^{-7}$. The significance threshold differs slightly between analyses because of variations in the number of genes with a minimum number of annotations.
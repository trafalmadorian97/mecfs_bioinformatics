# Linkage Disequilibrium
_Linkage disequilibrium_ (LD) refers to statistical dependence between genetic variants. LD is central to statistical genomics.

## Measures

- When we are interested in patterns of LD across an entire genomic region of length $n$, it is common to report the LD matrix $R\in\mathbb{R}^{n\times n}$, whose $(i,j)$ component is $r_{i,j}$ the [Pearson correlation](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) between variant $i$ and variant $j$.
- Note, however, that this matrix $R$ reflects only pairwise dependence, and so is not a complete characterization of LD.  In particular, there are many higher-order dependence structures consistent with any given $R$ matrix.

## Drivers of LD

There are two main physical processes that cause LD patterns: mutation and recombination.

### Mutation

For simplicity, first consider LD in the  absence of recombination, as in the case in mitochondrial DNA and certain regions of the Y chromosome.  In such recombination-free regions, the genomic distance between two variants is irrelevant to their LD. Instead, LD is a function of historical mutation patterns.  

Figure 7 from the Hapmap paper[@international2005haplotype] illustrates the concept:


![hapmap-mutation-fig](https://github.com/user-attachments/assets/d56eb383-5edf-4601-b4c5-a991acd25931)



### Recombination


Besides mutation, the other major driver of LD patterns in the Eukaryotic genome is recombination.


Some facts about LD:

- LD decays as the distance between variants increases, because the odds of an intervening [recombination event](https://en.wikipedia.org/wiki/Genetic_recombination) correspondingly increase. However, due to the complex structure of Eukaryotic DNA, the odds of recombination events are non-uniform across a chromosome. Thus, the rate of LD decay with genomic distance is not constant.
- Genetic variants that are relatively recent tend to have low frequency in the population, and thus low LD with all other variants, regardless of distance.


As an illustrative example, here is a plot of the absolute value of the correlation between genetic variants in a region of chromosome 1.  This plot was generated from the [UK Biobank LD matrices stored on AWS OpenData](https://registry.opendata.aws/ukbb-ld/).  In the plot, the x and y axes correspond genomic position, while color indicates absolute correlation.


![ld_example_plot](https://github.com/user-attachments/assets/a05681d5-91f3-4b89-8023-d3d50a22b8bd)

Consistent with the facts above, we observe irregularly spaced LD blocks. 


## Genomic Distance


## Genotyping
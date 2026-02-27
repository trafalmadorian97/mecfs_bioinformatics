
# GTEx Bulk RNA-seq

The Genotype Tissue Expression (GTEx) project[@gtex2020gtex] generates useful datasets via the analysis of diverse tissues from post-mortem human donors.

One important class of GTEx dataset records levels of RNA expression across a wide range of tissues.  These GTEx RNAseq datasets can help us discover the tissues in which a particular gene plays an important role.  For instance, if we find that across many donors, the RNA transcripts of a gene are over-expressed in the colon relative to other tissues, it is reasonable to assume that the gene plays its most important role in the colon.

GTEx bulk tissue RNA-seq datasets are used by the gene set analysis step of [MAGMA](../Bioinformatics_Concepts/MAGMA_Overview.md) to construct hypotheses about which tissues are most involved in a trait or disease.  Continuing the above example, MAGMA works on the principle that if most of the genes strongly associated with a disease are over-expressed in the colon according to GTEx RNAseq data, it is likely that the colon is central to the disease process.


[//]: # (## The RNA-seq Vector Space)

[//]: # ()
[//]: # ()
[//]: # (Besides allowing us to reason about which tissues are central to which disease processes, the GTEx bulk tissue RNA-seq datasets also generate a notion of similarity between tissues.  Two tissues are similar if, across most donors, their patterns of RNA expression are similar.)

[//]: # (RNAseq transcript-per-million measurements are distributed roughly according a zero-inflated lognormal distribution.  To compare tissues via standard statistical techniques, it is useful to pre-transform these measurements so that they follow a more conventional distribution. One approach is:)

[//]: # ()
[//]: # (1.  Compute the median level of RNA expression for each gene/tissue pair, measured on transcripts-per-million &#40;TPM&#41;.  Winsorize the levels to a maximum of 50 TPM.)

[//]: # (2. Log transform these median TPM values with a pseudocount of 1.  i.e.: $y=\log&#40;x+1&#41;$.)

[//]: # (3. Associate with each tissue the vector of these log-transformed median TPM values.  Two tissue are similar if their vectors in are close in this [vector space.]&#40;https://www.amazon.ca/Finite-Dimensional-Vector-Spaces-Paul-Halmos/dp/178139573X/ref=sr_1_2?crid=1NP4YJ625N57Q&dib=eyJ2IjoiMSJ9.RvuvtK5wnnXaXwjSyhb2f2Rew81JmSRnjAm5_9lOLvHKA8ao96xae_g4QQ_KrW7ae8ooj39H8M3cS_I45y-PJ0352qCDsvDR3iNLRKdx1IOD_7hx63eAVrzPWERq6ClWyxAXNsq1YpV0YPayj7MihW2ASQilGq76qolAt7bC1EeMiPlalsEA8gCET3a1CBzA0tb76Xt8IgF5PGgs9R0mS-R9sVOVgjTKYIf1bsUEU1VDiAkbZXBW2HwehlgzGozxuls4FRmWPT1HrygTdM1Uw1j34aaEHGuVXAgC6Gx7L_Y.2Bb5clC88ItUNoLhsRy6M94O8SUYlI67ubbO-Fekq3k&dib_tag=se&keywords=finite+dimensional+vector+spaces&qid=1762831919&sprefix=finite+dimensional+vector+space%2Caps%2C118&sr=8-2&#41;.)

[//]: # ()
[//]: # (This approach to defining tissue similarity is based on the one used by [FUMA]&#40;https://fuma.ctglab.nl/tutorial#snp-1gene&#41;.  )

[//]: # ()
## Artifacts in GTEx RNAseq data

Finucane et al.[@finucane2018heritability] applied Gene Ontology Enrichment Analysis to the genes expressed in lung tissue in the GTEx bulk RNAseq dataset.  They found strong enrichment for immune-related genes.  They hypothesized that this enrichment of immune genes in the lung was not a real property of lung tissue, but was instead a consequence of the presence of blood in the lungs of the post-mortem donors from which GTEx samples were collected. Users of GTEx datasets need to be aware of the presence of this kind of data artifact.

[//]: # (![rnaseq_pca]&#40;https://github.com/user-attachments/assets/02023612-3f4a-4714-940e-347acd20b054&#41;)



[//]: # (## Key Issues)

[//]: # ()
[//]: # (There are several objections that could be reasonably raised to the above:)

[//]: # (- It does not inevitably follow that the tissue in which a gene's RNA is mostly highly transcribed is the tissue in which the gene serves its most important role.  One could imagine that a gene is expressed at an extremely low level in some tissue, but nevertheless plays a vital role.)

[//]: # ()
[//]: # ()

## References

[FUMA Tutorial](https://fuma.ctglab.nl/tutorial#snp-1gene) (Discusses use of GTEx data for post-GWAS analysis)

[GTEx Bulk Tissue Expression Data](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression)

[GTEx paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7737656/)


[Alberts, Bruce, et al. Molecular biology of the cell: seventh edition. WW Norton & Company, 2022. See Chapter 6.](https://wwnorton.com/books/9780393884821)
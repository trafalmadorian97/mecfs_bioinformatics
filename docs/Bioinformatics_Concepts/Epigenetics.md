## Introduction

Neglecting somatic mutations, all non-gamete cells in the body contain the same DNA.  Despite this, the cells of multicellular organisms are diverse: adipocytes differ from neurons, which differ from macrophages, which differ from B cells, which differ from cardiomyocytes.

What explains this diversity? Chiefly, the answer is gene expression.  While all genes are present in the DNA of all cells, different genes are switched on or off in different cell types.


Gene expression is partly controlled via structural and chemical modifications of DNA that do not change the nucleotides sequence, but do affect transcription.  These structural and chemical modifications are referred to collectively as "Epigenetics".

## Histone marks

One of the key epigenetic regulatory systems involves histone marks.

Eukaryotic DNA is tightly wound around spool-like complexes of proteins called histones.  The resulting DNA-histone units are called nucleosomes.

![histone_image_from_wikipedia](https://github.com/user-attachments/assets/077abc51-5dfa-4891-8fa1-e83db80976eb)

[Image Source: David O Morgan](https://en.wikipedia.org/wiki/Nucleosome#/media/File:Basic_units_of_chromatin_structure.svg).


Histones can covalently bond to a variety of molecules.  These bonds are called "histone marks".  Histone marks affect the three-dimensional structure of chromatin.  For example, certain histone marks cause DNA to be tightly wound, so that the genes it contains cannot be transcribed, while other histone marks cause DNA to open up, allowing genes to be readily transcribed (see Alberts et al. pg 220[@alberts2022molecular]).


Here are some examples of histone marks with known effects:


| Histone Mark | Effect on Chromatin                                    |
|:-------------|:-------------------------------------------------------|
| H3K4me1      | [Associated with enhanced transcription](https://en.wikipedia.org/wiki/H3K4me1)|
| H4K4me3      | Open, accessible chromatin[@alberts2022molecular]                              |
| H3K9ac       | Open, accessible chromatin[@alberts2022molecular]                              |
| H3K9me3      | Closed, inaccessible chromatin[@alberts2022molecular]                          |
| [H3K27ac]      | [Associated with higher activation of transcription](https://en.wikipedia.org/wiki/H3K27ac) |
| H3K27me3     | Closed, inaccessible chromatin[@alberts2022molecular]                          |  


[//]: # (Thus, histone marks constitute part of the answer to the question of how gene expression is controlled.)


Assays that report histone marks are a key method for understanding epigenetic regulation.


## ATAC-seq

Besides assays that report histone marks, assays that directly measure chromatin accessibility are also informative about epigenetic regulation. ATAC-seq is an important example of such an accessibility assay.

ATAC-seq  repurposes [transposase](https://en.wikipedia.org/wiki/Transposase), the enzyme responsible for the movement of transposons, or "jumping genes".  It works on the principle that transposase can interact with open chromatin, but not closed chromatin.  This explains the name: Assay for Transposase Accessible Chromatin (ATAC).
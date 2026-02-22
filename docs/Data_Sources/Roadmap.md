# The Roadmap Epigenetics Project

## Overview

The Roadmap Epigenetics Project[@bernstein2010nih] was an NIH-funded initiative whose aim was to gather [epigenetic data](../Bioinformatics_Concepts/Epigenetics.md) such as:

- DNA methylation data
- [Histone modification data](../Bioinformatics_Concepts/Epigenetics.md#histone-marks)
- Chromatin accessibility data (via the DNase 1 chromatin accessibility assay),

from a wide variety of human tissue.

A key feature of the Roadmap Project was its adherence to the principle of *In Vivo Veritas*.  Whereas previous epigenetics studies had studied immortalized cell lines that may not faithfully recapitulate normal biology, the Roadmap Project studied normal cells retrieved from donors.  


## Utility

Roadmap datasets tell us which regions of the genome are under what kind of epigenetic control in which tissues.  They can thus be combined with GWAS results via tools like [S-LDSC](../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md) and [MAGMA](../Bioinformatics_Concepts/MAGMA_Overview.md). Essentially, if we find a GWAS hit in a particular region of the genome, we can use Roadmap to try to ask: "In which tissues is that region of the genome epigenetically up-regulated?"


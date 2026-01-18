# UKBB PPP

The UK Biobank Pharma Proteomics Project (PPP)[@sun2023plasma] is an effort to characterize the plasma proteomes of a sample of UK Biobank participants.  The PPP was funded an organized by a consortium consisting of the research arms of major pharmaceutical companies.  The aim was to produce datasets suitable for drug discovery research.

## Motivation

Often, GWAS hits are mechanistically mysterious: it is unclear what biological pathway links a given genetic polymorphism with an associated phenotype.  At a high level the PPP is motivated by the observation that since proteins are the fundamental biological effector molecules, they can help explain GWAS hits.  If we observe that a GWAS hit for a phenotype is also associated with changed plasma levels of a protein, we may reason that the protein lies along the biological pathway leading to the GWAS hit.

## Public Data
 
While access to the raw individual-level PPP data requires UKBB accreditation, there are a number of powerful PPP datasets that have been publicly released.

### pQTLs

A plasma pQTL for a protein is a genetic variation that has a statistical association with plasma levels of that protein. One can generate plasma pQTLs by running GWAS with the plasma level of a target protein as the phenotype, then applying fine-mapping to each hit.   The supplementary information of the original Sun et al. paper[@sun2023plasma] includes tables of pQTLs for 2525 different proteins (see Supplementary Tables 9 and 10).  These pQTLs are divided into cis-pQTLs (within 1Mb of the gene encoding the protein they affect), and trans-pQTLs.  

A key use of pQTLs is in [Mendelian randomization](../../Bioinformatics_Concepts/Mendelian_Randomization.md).  One runs an MR study with pQTLs for a protein as instruments, the level of the protein as the exposure, and a phenotype of interest as the outcome.  This approach can generate evidence that the level of the protein is causally related to the genesis of the phenotype.  In this kind of MR study, one must decide whether to use exclusively cis-pQTLs, or both trans- and cis- pQTLs.  One the on hand, using only cis-pQTLs reduces the risk of violations of the MR assumptions due to pleiotropy.  On the other hand, for certain proteins inclusion of trans-pQTL can increase statistical power.


The figure below from Donoghue et al.[@donoghue2025integration] shows the proportion of variance explained in various asthma-relevant proteins by cis vs trans pQTL.  For certain proteins, it is necessary to use trans pQTLs to explain an adequate proportion of the variance.

![pqtl_cis_vs_trans](https://github.com/user-attachments/assets/0d3e4597-7e60-49e2-be07-bc01442fc3d2)


### GWAS Summary Statistics

In addition to tables of pQTLs, the PPP has also released full GWAS summary statistics for each protein measured, which can be downloaded from [Synapse](https://www.synapse.org/Synapse:syn51365303). These summary statistics can be used for colocalization analysis.  For example, after one finds MR evidence that levels of a protein causally affect a phenotype, one can use a high probability of colocalization between pQTLs for protein levels and GWAS hits of the phenotype to provide additional support for the causal connection[@zuber2022combining].


## Links

[Genetics Podcast Interview About PPP](https://www.youtube.com/watch?v=nSJJjU61JAk)
# UKBB PPP

The UK Biobank Pharma Proteomics Project (PPP)[@sun2023plasma] is an effort to characterize the plasma proteomes of a sample of UK Biobank participants.  The PPP was funded and organized by a consortium consisting of the research arms of major pharmaceutical companies.  The aim was to produce datasets suitable for drug-discovery research.

## Motivation

Often, GWAS hits are mechanistically mysterious: it is unclear what biological pathway links a given genetic polymorphism with an associated phenotype.  At a high level the PPP is motivated by the observation that since proteins are the fundamental biological effector molecules, they can help explain GWAS hits.  If we observe that a GWAS hit for a phenotype is also associated with changed plasma levels of a protein, we may reason that the protein lies along the biological pathway leading to the GWAS hit.

## Public Data
 
While access to individual-level PPP data requires UKBB accreditation, there are a number of powerful summary PPP datasets that have been publicly released.

### pQTLs

A plasma pQTL for a protein is a genetic variation that has likely-causal association with variations in plasma levels of that protein. One can find such pQTLs by running a GWAS with the plasma level of the target protein as the phenotype, then fine-mapping the results.   The original Sun et al. paper[@sun2023plasma] provides pQTLs for 2525 different proteins (see Supplementary Tables 9 and 10).  These pQTLs are divided into cis-pQTLs (within 1Mb of the gene encoding the protein they affect), and trans-pQTLs.  

pQTLs can be used in  [Mendelian randomization](../Bioinformatics_Concepts/Mendelian_Randomization.md) (MR). In pQTL-based MR, the pQTLs for a protein are instruments, the level of the protein is the exposure, and a phenotype is the outcome.  The goal is causal inference about the effect of the protein on the phenotype.  In a pQTL MR study, one must decide whether to exclusively use cis-pQTLs, or to use both trans- and cis- pQTLs.  On the one hand, using only cis-pQTLs reduces the risk of violations of the MR assumptions due to pleiotropy, since the pathway from the genetic variant to the protein level is more straightforward for cis-pQTLs.  On the other hand, for certain proteins, inclusion of trans-pQTL can increase statistical power.


The figure below from Donoghue et al.[@donoghue2025integration] shows the proportion of variance in levels of asthma-relevant proteins explained by cis vs trans pQTL.  For certain proteins, it is necessary to use trans pQTLs to explain an adequate  of the variance.

![pqtl_cis_vs_trans](https://github.com/user-attachments/assets/0d3e4597-7e60-49e2-be07-bc01442fc3d2)


### GWAS Summary Statistics

In addition to tables of pQTLs, the PPP has also released full GWAS summary statistics for each protein measured. There are available at [Synapse](https://www.synapse.org/Synapse:syn51365303). These summary statistics are useful for colocalization analysis.  For example, after one finds MR evidence suggesting that levels of a protein causally affect a phenotype, one can use colocalization to further support or refute a causal connection[@zuber2022combining].


## Links

[Genetics Podcast Interview About PPP](https://www.youtube.com/watch?v=nSJJjU61JAk)
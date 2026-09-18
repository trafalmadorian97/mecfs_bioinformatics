# Sumstats Standardization

It is often useful to convert GWAS summary statistics to a standard form before analyzing them. The use of a standard form facilitates comparisons between GWAS, and is expected by certain post-GWAS analysis tools.  Summary-statistics standardization a complex, multi-step process.

## Strand Standardization

In the DNA double helix, the two strands are complementary: adenine (A) always pairs with thymine (T), and cytosine (C) always pairs with guanine (G)[^judson_note].  This means that the same variant can be represented equivalently by describing its base sequence on either strand.

For instance, suppose that a GWAS summary statistics file contains the following data: 

|   CHROM |   GENPOS | ALLELE0   | ALLELE1   |    A1FREQ |         BETA |     LOG10P |
|--------:|---------:|:----------|:----------|----------:|-------------:|-----------:|
|       5 | 90000126 | G         | A         | 0.0154907 |  0.0388003   | 0.346247   |
|       5 | 90002604 | T         | TG        | 0.16074   | -0.000160003 | 0.00323221 |
|       5 | 90002963 | G         | C         | 0.0109188 | -0.0368269   | 0.257525   |
|       5 | 90003830 | G         | A         | 0.160737  | -0.00019273  | 0.0038963  |

This data can be represented equivalently on the other strand as 


|   CHROM |   GENPOS | ALLELE0   | ALLELE1 |    A1FREQ |         BETA |     LOG10P |
|--------:|---------:|:----------|:--------|----------:|-------------:|-----------:|
|       5 | 90000126 | C         | T       | 0.0154907 |  0.0388003   | 0.346247   |
|       5 | 90002604 | A         | CA      | 0.16074   | -0.000160003 | 0.00323221 |
|       5 | 90002963 | C         | G       | 0.0109188 | -0.0368269   | 0.257525   |
|       5 | 90003830 | C         | T       | 0.160737  | -0.00019273  | 0.0038963  |

(Note that at multi-base alleles like the indel T/TG, we flip the order of the bases.)

While GENPOS always refers to genomic coordinates in the $5'\to 3'$ direction along the forward strand, summary statistics files have been known to report alleles on either the forward or the reverse strand[@hartwig2016two]. This ambiguity is problematic when one seeks to jointly analyze two GWAS, or look up variants from a GWAS in a reference.

In most cases, the ambiguity can be resolved by consulting a reference genome sequence [FASTA file](https://en.wikipedia.org/wiki/FASTA_format).  For each genetic variant, one compares the alleles to the reference sequence at the given genomic coordinate. If there is a match, we know the alleles are correctly represented on the forward strand.  If instead there is match with the complement of one the alleles, we know the alleles are represented on the reverse strand, and we can standardize the genetic variant by substituting the alleles for their complements.  If there is no match between the alleles and either the forward or reverse strand, this indicates an error.


### Palindromic variants

Consider the C/G variant in the table above.  At position 90002963, the reference genome (GRCh38) has the base G.  There are two possible interpretations:

- ALLELE 1 is the reference allele, while ALLELE 0 is the alternate allele, and the alleles refer to the forward strand.
- ALLELE 0 is the reference allele, while ALLELE 1 is the alternate allele, and the alleles refer to the reverse strand.


It is impossible to distinguish these cases using only CHROM/GENPOS/ALLELE0/ALLELE1. This is an example what is called a _palindromic variant_.


If one has access to a reference database of standardized forward-strand genetic variants with allele frequencies, one can sometime disambiguate palindromic variants using these frequencies. Suppose as above that we have a palindromic variant and ALLELE1 matches the reference allele in the database. The procedure is:


- If the allele frequency of ALLELE1 is close to 0.5, we can't disambiguate. Otherwise:
- Compare the allele frequency of ALLELE1 to the frequency of the reference allele in the database.  If there is a close match, we can assume that the allele pair in the summary statistics is on the forward strand.
- Compare the allele frequency of ALLELE1 to one minus the frequency of the reference allele in the database.  If there is a close match, we can assume the allele pair in the summary statistics is on the reverse strand.


## Ambiguous Indels


Consider the following summary statistics row



| CHROM |          GENPOS | ALLELE0 | ALLELE1 |        A1FREQ | 
|------:|----------------:|:--------|:--------|--------------:|
|    15 |        54698192 | TG      | T       |         0.845 |


Suppose that due to incomplete documentation (which is not uncommon for summary statistics), we do not know whether ALLELE0 is the reference allele (so that the variant is a deletion), or whether ALLELE1 is the reference allele (so that the variant is an insertion).

Unfortunately, just comparing against a FASTA file containing the hg19 reference sequence cannot resolve this.  To see why, note that at hg19 genomic position 54698192 on chromosome 15, the reference sequence has the value TGGGGGG. Thus both T and TG match the reference sequence at genetic position 15:54698192. As in the case of palindromic variants, we can only confidently resolve ambiguous indels using allele frequency information from a population database for a matched population.




## Left Normalization

Even if we restrict ourselves to the forward strand, and even if we fix which allele is the reference allele, the same genetic variant can still be represented in multiple ways.

This principle is demonstrated by the following examples from the [Center for Statistical Genetics at the University of Michigan:](https://genome.sph.umich.edu/w/index.php?title=Variant_Normalization)

![trim_example_normalization](https://github.com/user-attachments/assets/2f87c2be-4e00-481b-b130-1df7554029cf)


![example_2_normalization](https://github.com/user-attachments/assets/962344dd-4ce7-42a1-8912-27448c96cfc7)


What is needed is a unique canonical form to which all equivalent variants can be normalized.  Fortunately, there is an algorithm[@tan2015unified] called *left-normalization* that converts variant to such a canonical form. While many modern GWAS pre-left-normalize their alleles, it is often advisable to run the left-normalization algorithm, just to be sure.

### Algorithm details


todo



[^judson_note]:For a discussion of the structure of DNA, see Chapter 4 of _Molecular Biology of the Cell_[@alberts2022molecular].  For a historical account of the discovery of this structure, see Judson's book[@judson1996eighth].

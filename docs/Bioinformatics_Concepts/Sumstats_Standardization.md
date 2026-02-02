# Sumstats Standardization

Before analyzing GWAS summary statistics, it is often desirable to bring them to a standard form.  This is a complex, multi-step process.

## Step 1: Strand Standardization

In the double-stranded DNA molecule, the two strands are complementary: adenine (A) always pairs with thymine (T), and cytosine (C) always pairs with guanine (G).  This means that the same variant can be represented equivalently by describing its base sequences on either strand.

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

(Note that at multi-base alleles like the indel T/TG, we flip the order of the bases in addition to computing their complements)

While GENPOS always refers to genomic coordinate in the $5'\to 3'$ direction along the forward strand, summary statistics files have been known report alleles on either the forward or the reverse strand[@hartwig2016two]. This ambiguity is problematic when one seeks to jointly analyze two GWAS, or look up variants from a GWAS in a reference.

In most cases, the ambiguity can be resolved by consulting a reference genome sequence [FASTA file](https://en.wikipedia.org/wiki/FASTA_format).  For each genetic variant, one compares the alleles to the reference sequence at the given genomic coordinate. If there is a match, we know the alleles are correctly represented on the forward strand.  If instead there is match with the complement of one the alleles, we know the alleles are represented on the reverse strand.  We can standardize the genetic variant by substituting the alleles for their complements.  If there is no match between the alleles and either the forward or reverse strand, this likely indicates an error.


### Palindromic variants

Consider the C/G variant in the table above.  At position 90002963, the reference genome (GRCh38) has the base G.  There are two possible interpretations:

- ALLELE 1 is the reference allele, while ALLELE 0 is the alternate allele, and the alleles refer to the forward strand.
- ALLELE 0 is the reference allele, while ALLELE 1 is the alternate allele, and the alleles refer to the reverse strand.


It is impossible to distinguish these cases using only CHROM/GENPOS/ALLELE0/ALLELE1.


If one has access to a reference database of standardized genetic variants with allele frequencies, one can sometime disambiguate using these frequencies. Suppose as above that we have a palindromic variant and ALLELE 1 matches the reference allele in the database. The procedure is:

- Compare the allele frequency of ALLELE1 to the frequency of the reference allele in the database.  If there is a close match, we can assume that the allele pair in the summary statistics are on the forward strand.
- Compare the allele frequency of ALLELE1 to one minus the frequency of the reference allele in the database.  If there is a close match, we can assume the allele pair in the summary statistics are on the reverse strand.
- If the allele frequency of ALLELE1 is close to 0.5, we can't disambiguate.




## Step 2: Left Normalization

Even if we restrict ourselves to the forward strand, the same genetic variant can still be represented in multiple ways.

This principle is demonstrated by the following examples from the [Center for Statistical Genetics at the University of Michigan:](https://genome.sph.umich.edu/w/index.php?title=Variant_Normalization)

![trim_example_normalization](https://github.com/user-attachments/assets/2f87c2be-4e00-481b-b130-1df7554029cf)


![example_2_normalization](https://github.com/user-attachments/assets/962344dd-4ce7-42a1-8912-27448c96cfc7)


Clearly, what is needed is a unique canonical form to which all equivalent variants can be normalized.  Fortunately, there is a left-normalization algorithm that does just this. While many modern GWAS pre-left-normalize their alleles, it is often advisable to run left-normalization, just to be sure.

### Algorithm details


todo


## Step 3: Harmonization between GWAS, or between a GWAS and an LD Matrix

todo


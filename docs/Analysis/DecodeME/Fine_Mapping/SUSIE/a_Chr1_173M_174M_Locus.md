---
hide:
- navigation
- toc
---
# Chr1 173.5M-174.5M


## Methodology

In an attempt to narrow the DecodeME[@genetics2025initial] GWAS-1 signal, I [fine-mapped](../../../../Bioinformatics_Concepts/Fine_Mapping.md) the hit on chromosome 1 using [SUSIE](https://stephenslab.github.io/susieR/)[@wang2020simple].

As a linkage disequilibrium reference, I used a UK Biobank LD matrix hosted on [AWS Open Data](https://registry.opendata.aws/ukbb-ld/).  Because this LD reference uses GRCh37 coordinates, I used GWASLab liftover to translate the DecodeME summary statistics to GRCh37.


As a sensitivity analysis, I ran SUSIE 4 times:

- Once with $L=10$
- Once with $L=2$
- Once with $L=1$
- Once with $L=10$ and strict variant filtering

Here $L$ refers to the maximum number of credible sets that can found by SUSIE.  A lower value of $L$ corresponds to increased regularization, since it decreases the ability of SUSIE to use extra credible sets to fit noise.  Weissbrod et al.[@weissbrod2020functionally] observe that setting $L$ to 1 protects against instability or inaccuracy due to a mispatch between the LD reference population and the GWAS population, because when $L=1$, SUSIE no longer depends on the LD matrix.   They also point out that when $L=2$, even though SUSIE still depends on the LD reference population, empirically tends to be robust to moderate levels of population mismatch. I thus used the $L=1$ and $L=2$ runs to evaluate whether LD mismatch issues could be influencing SUSIE's results.


"Variant filtering" refers to removal of variants with a high z-score and a high likelihood ratio according to a [Kriging](https://en.wikipedia.org/wiki/Kriging)-based likelihood ratio test.  Zou et al.[zou2022fine] propose this filtering strategy to mitigate instability in SUSIE due to mismatch between the LD and GWAS populations.  In the first three runs above, I filter variants with $LR\ge 2$ and $|z|\ge 2$, [consistent with the SUSIE documentation](https://stephenslab.github.io/susieR/reference/kriging_rss.html).  In the final run I instead filter variants with $LR\ge 2$ and $|z|\ge 1$, to evaluate the sensitivity of the results to the filtering threshold.

In my SUSIE runs, I retained palindromic SNPs whose orientation GWASLAB was able to determine from allele frequencies in the Thousand Genomes Project reference VCF file.


## Results

In all 4 runs, SUSIE found a single diffuse credible set.  Moreover, this credible set contained the same 86 variants in all four runs, as illustrated in the UpSet plot below:


![upset_chrom_1](https://github.com/user-attachments/assets/0a4ff4a0-0925-403d-ad3b-2a512f406c46)


The next figure illustrates the SUSIE results for $L=10$. It is representative of the other runs.

![chr1_stackplot](https://github.com/user-attachments/assets/78104e68-7402-46b9-8ddf-70d83d9e4b60)

- The top panel is a heatmap in which pixel (i,j) is colored according to the squared correlation between variants i and j.  The heatmap thus reveals the local linkage disequilibrium (LD) structure in the vicinity of the GWAS hit, which is a determinant of SUSIE's results when $L>1$.

- The second panel shows a local Manhattan plot.

- The third panel shows SUSIE posterior inclusion probability (PIP).

- The bottom panel shows genes in the region of the GWAS hit. 

Overall, SUSIE has returned a diffuse signal in a region with a number of plausible genes.


The table below lists the full detailed SUSIE results for the $L=10$ case

??? info "Variant List"

    | cs   |   CHR |       POS | EA    | NEA   |      alpha |         mu |        PIP |
    |:-----|------:|----------:|:------|:------|-----------:|-----------:|-----------:|
    | L1   |     1 | 173815290 | C     | T     | 0.0349532  | -0.022181  | 0.0349532  |
    | L1   |     1 | 173853127 | C     | T     | 0.0334384  | -0.0221481 | 0.0334384  |
    | L1   |     1 | 173865586 | T     | C     | 0.0332454  | -0.0221439 | 0.0332454  |
    | L1   |     1 | 173815111 | C     | T     | 0.0332354  | -0.0221436 | 0.0332354  |
    | L1   |     1 | 173866074 | A     | T     | 0.0332311  | -0.0221435 | 0.0332311  |
    | L1   |     1 | 173878862 | C     | T     | 0.0328582  | -0.0221352 | 0.0328582  |
    | L1   |     1 | 173812639 | A     | C     | 0.0327682  | -0.0221331 | 0.0327682  |
    | L1   |     1 | 173851310 | A     | G     | 0.0313245  | -0.0220996 | 0.0313245  |
    | L1   |     1 | 173859100 | G     | A     | 0.0296411  | -0.0220585 | 0.0296411  |
    | L1   |     1 | 173863209 | A     | G     | 0.0296014  | -0.0220575 | 0.0296014  |
    | L1   |     1 | 173863569 | A     | T     | 0.0292311  | -0.0220482 | 0.0292311  |
    | L1   |     1 | 173863567 | T     | G     | 0.0292311  | -0.0220482 | 0.0292311  |
    | L1   |     1 | 173846590 | G     | T     | 0.0276321  | -0.0220062 | 0.0276321  |
    | L1   |     1 | 173863568 | T     | A     | 0.0275657  | -0.0220044 | 0.0275657  |
    | L1   |     1 | 173832336 | T     | C     | 0.0272003  | -0.0219944 | 0.0272003  |
    | L1   |     1 | 173878832 | C     | T     | 0.0261331  | -0.0219645 | 0.0261331  |
    | L1   |     1 | 173838788 | T     | TG    | 0.0261228  | -0.0219642 | 0.0261228  |
    | L1   |     1 | 173855298 | T     | A     | 0.0255997  | -0.0219491 | 0.0255997  |
    | L1   |     1 | 173846110 | A     | G     | 0.025353   | -0.0219418 | 0.025353   |
    | L1   |     1 | 173857283 | G     | A     | 0.0243066  | -0.0219102 | 0.0243066  |
    | L1   |     1 | 173848009 | G     | A     | 0.023978   | -0.0219    | 0.023978   |
    | L1   |     1 | 173824813 | T     | C     | 0.0234489  | -0.0218833 | 0.0234489  |
    | L1   |     1 | 173842467 | G     | A     | 0.0230966  | -0.0218719 | 0.0230966  |
    | L1   |     1 | 173870321 | G     | GTAC  | 0.0230538  | -0.0218705 | 0.0230538  |
    | L1   |     1 | 173857037 | T     | C     | 0.0228147  | -0.0218627 | 0.0228147  |
    | L1   |     1 | 173881871 | T     | C     | 0.0203149  | -0.0217753 | 0.0203149  |
    | L1   |     1 | 173844051 | T     | A     | 0.0184477  | -0.0217024 | 0.0184477  |
    | L1   |     1 | 173820365 | C     | T     | 0.0145082  | -0.0215198 | 0.0145082  |
    | L1   |     1 | 173832772 | CA    | C     | 0.0144368  | -0.021516  | 0.0144368  |
    | L1   |     1 | 173878471 | G     | A     | 0.0118415  | -0.0213641 | 0.0118415  |
    | L1   |     1 | 173831882 | G     | A     | 0.00912575 | -0.0211628 | 0.00912575 |
    | L1   |     1 | 173743879 | CAAAA | C     | 0.00696391 | -0.0209519 | 0.00696391 |
    | L1   |     1 | 173717200 | ACT   | A     | 0.00675049 | -0.0209275 | 0.00675049 |
    | L1   |     1 | 173767443 | T     | A     | 0.00660937 | -0.0209109 | 0.00660937 |
    | L1   |     1 | 173783493 | C     | T     | 0.00650478 | -0.0208983 | 0.00650478 |
    | L1   |     1 | 173699007 | G     | A     | 0.00592378 | -0.0208246 | 0.00592378 |
    | L1   |     1 | 173698510 | T     | C     | 0.00585287 | -0.0208151 | 0.00585287 |
    | L1   |     1 | 173709616 | G     | T     | 0.0056909  | -0.020793  | 0.0056909  |
    | L1   |     1 | 173734270 | CAACA | C     | 0.0056717  | -0.0207903 | 0.0056717  |
    | L1   |     1 | 173683954 | T     | C     | 0.00516697 | -0.0207165 | 0.00516697 |
    | L1   |     1 | 173755936 | TGAAG | T     | 0.00476139 | -0.0206516 | 0.00476139 |
    | L1   |     1 | 174210076 | T     | TTG   | 0.00289739 | -0.0202525 | 0.00289739 |
    | L1   |     1 | 174128994 | G     | A     | 0.00239407 | -0.0200971 | 0.00239407 |
    | L1   |     1 | 174158856 | C     | T     | 0.00235282 | -0.0200829 | 0.00235282 |
    | L1   |     1 | 174111115 | T     | C     | 0.00229859 | -0.0200638 | 0.00229859 |
    | L1   |     1 | 174066947 | T     | C     | 0.00227798 | -0.0200564 | 0.00227798 |
    | L1   |     1 | 174069469 | CT    | C     | 0.00227531 | -0.0200555 | 0.00227531 |
    | L1   |     1 | 174084104 | A     | G     | 0.0022655  | -0.0200519 | 0.0022655  |
    | L1   |     1 | 174152688 | G     | A     | 0.00225032 | -0.0200464 | 0.00225032 |
    | L1   |     1 | 174146656 | C     | A     | 0.00222885 | -0.0200385 | 0.00222885 |
    | L1   |     1 | 174085043 | TACA  | T     | 0.00219654 | -0.0200266 | 0.00219654 |
    | L1   |     1 | 174076864 | C     | G     | 0.00218651 | -0.0200228 | 0.00218651 |
    | L1   |     1 | 174064481 | A     | T     | 0.0021761  | -0.0200189 | 0.0021761  |
    | L1   |     1 | 174068049 | A     | G     | 0.00215682 | -0.0200116 | 0.00215682 |
    | L1   |     1 | 174069981 | C     | T     | 0.00215307 | -0.0200102 | 0.00215307 |
    | L1   |     1 | 174191694 | C     | T     | 0.00214725 | -0.0200079 | 0.00214725 |
    | L1   |     1 | 174062911 | C     | T     | 0.00214307 | -0.0200063 | 0.00214307 |
    | L1   |     1 | 174057626 | T     | C     | 0.00214136 | -0.0200057 | 0.00214136 |
    | L1   |     1 | 174241165 | C     | G     | 0.0021385  | -0.0200046 | 0.0021385  |
    | L1   |     1 | 174093837 | G     | T     | 0.00213651 | -0.0200038 | 0.00213651 |
    | L1   |     1 | 174157834 | A     | G     | 0.00213405 | -0.0200029 | 0.00213405 |
    | L1   |     1 | 174170341 | C     | T     | 0.00213381 | -0.0200028 | 0.00213381 |
    | L1   |     1 | 174086292 | T     | C     | 0.00211439 | -0.0199953 | 0.00211439 |
    | L1   |     1 | 174050667 | G     | A     | 0.00210878 | -0.0199931 | 0.00210878 |
    | L1   |     1 | 174068929 | A     | G     | 0.00207848 | -0.0199812 | 0.00207848 |
    | L1   |     1 | 174079700 | G     | C     | 0.00207189 | -0.0199786 | 0.00207189 |
    | L1   |     1 | 174118271 | C     | G     | 0.00206586 | -0.0199762 | 0.00206586 |
    | L1   |     1 | 174075924 | G     | T     | 0.00206197 | -0.0199746 | 0.00206197 |
    | L1   |     1 | 174075925 | C     | T     | 0.00206197 | -0.0199746 | 0.00206197 |
    | L1   |     1 | 174059522 | C     | G     | 0.002049   | -0.0199694 | 0.002049   |
    | L1   |     1 | 174161430 | A     | G     | 0.00203663 | -0.0199645 | 0.00203663 |
    | L1   |     1 | 174072203 | G     | A     | 0.00203588 | -0.0199642 | 0.00203588 |
    | L1   |     1 | 174078417 | A     | G     | 0.00202943 | -0.0199615 | 0.00202943 |
    | L1   |     1 | 174079584 | T     | A     | 0.00202752 | -0.0199608 | 0.00202752 |
    | L1   |     1 | 174079585 | G     | A     | 0.00202752 | -0.0199608 | 0.00202752 |
    | L1   |     1 | 174069421 | T     | C     | 0.00191843 | -0.0199152 | 0.00191843 |
    | L1   |     1 | 174176363 | G     | A     | 0.00190403 | -0.019909  | 0.00190403 |
    | L1   |     1 | 174144293 | CACAT | C     | 0.00189623 | -0.0199056 | 0.00189623 |
    | L1   |     1 | 174197408 | G     | A     | 0.00187952 | -0.0198983 | 0.00187952 |
    | L1   |     1 | 174078860 | C     | T     | 0.00187677 | -0.0198971 | 0.00187677 |
    | L1   |     1 | 174113408 | G     | A     | 0.0018707  | -0.0198944 | 0.0018707  |
    | L1   |     1 | 174152476 | C     | T     | 0.0018538  | -0.0198869 | 0.0018538  |
    | L1   |     1 | 174150850 | C     | T     | 0.00182508 | -0.019874  | 0.00182508 |
    | L1   |     1 | 174060062 | G     | A     | 0.00180084 | -0.0198629 | 0.00180084 |
    | L1   |     1 | 174326702 | G     | A     | 0.00161599 | -0.0197731 | 0.00161599 |
    | L1   |     1 | 174288602 | G     | A     | 0.00158912 | -0.0197592 | 0.00158912 |⏎

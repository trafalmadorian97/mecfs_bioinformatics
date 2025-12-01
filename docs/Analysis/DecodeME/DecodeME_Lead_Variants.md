# DecodeME Lead Variants

As an initial analysis step, we apply [GWASLab's procedure for extracting lead variants](https://cloufield.github.io/gwaslab/ExtractLead/) to the DecodeME GWAS 1 data.  This procedure groups together significant genetic variants using a sliding-window approach, then reports the most significant variant in each region.

## Table of Variants from GWASLAB

Here is the table of lead variants produces by GWASLab:

|    | SNPID            |   CHR |       POS | EA   | NEA   |      EAF |       BETA |        SE |   CHISQ |   MLOG10P |      N | GENE    |
|---:|:-----------------|------:|----------:|:-----|:------|---------:|-----------:|----------:|--------:|----------:|-------:|:--------|
|   | 1:173846152:T:C  |     1 | 173846152 | C    | T     | 0.325279 | -0.0759185 | 0.0136323 | 31.0138 |   7.59142 | 275488 | DARS2   |
|   | 6:26239176:A:G   |     6 |  26239176 | G    | A     | 0.261233 |  0.0825251 | 0.0140356 | 34.5711 |   8.3862  | 275488 | H4C6    |
|   | 6:97984426:C:CA  |     6 |  97984426 | CA   | C     | 0.546314 | -0.068408  | 0.0125368 | 29.7742 |   7.3139  | 275488 | MMS22L  |
|   | 15:54866724:A:G  |    15 |  54866724 | G    | A     | 0.311707 |  0.0785482 | 0.0135977 | 33.369  |   8.11788 | 275488 | UNC13C  |
|   | 17:52183006:C:T  |    17 |  52183006 | T    | C     | 0.329679 |  0.0805953 | 0.0134577 | 35.8658 |   8.67492 | 275488 | CA10    |
|   | 20:48914387:T:TA |    20 |  48914387 | TA   | T     | 0.633808 |  0.0909054 | 0.0133414 | 46.4275 |  11.0219  | 275488 | ARFGEF2 |

Note that the genomic coordinates in the table above refer to genome build 38.

## Comparison with lead variants reported in DecodeME
The lead variants reported above by GWASLAB all agree with the variants reported by the DecodeME preprint.  However, in some cases, GWASLAB assigns these lead variants to different genes than GWASLAB.

In particular


-  17:52183006:C:T is assigned to CA10 by both GWASLAB and the DecodeME preprint.
-  20:48914387:T:TA is assigned to both ARFGEF2 and CSE1L by the DecodeME preprint.  DecodeME finds ARFGEF2.
-  The other variants are assigned to different genes by GWASLAB and the DecodeME preprint.  This reflects ambiguity about how to assign GWAS signals to genes.

## Reproducing
To reproduce this analysis, run the [initial DecodeME analysis script][mecfs_bio.analysis.decode_me_initial_analysis]



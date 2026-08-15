# Cloud Computation


Certain bioinformatics operations require much greater computational resources than are typically available on a local development machine.  Fitting a Bayesian model to GWAS summary statistics for the purposes of polygenic prediction or genome-wide fine mapping is an example of such an expensive bioinformation operation[^sbayesrc_note]. 

Since it running these computationally costly operations locally would be impractical, `mecfs-bio` Tasks instead run them in cloud.


## AWS




[^sbayesrc_note]: For instance, the SBayesRC fine-mapping paper[wu2026genome] reports that running SBayesRC on summary statistics with 13M bariants required *"150Gb of RAM and 13h of computation using 24 CPU cores"*
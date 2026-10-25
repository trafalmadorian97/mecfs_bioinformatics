"""
Region plot of DecodeME GWAS-1 around SLC6A2, the gene encoding the norepinephrine transporter (NET).

Centered on the midpoint of the SLC6A2 gene body (GRCh37 chr16:55,689,516-55,740,102).
Almost no variant in this window reaches P < 0.01, so all variants are drawn (skip=0) rather than only those with P < 0.01.
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_annovar_37_sumstats import (
    DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_region_plot_task_arbitrary_locus import (
    GwasLabRegionPlotTargetLocusTask,
)

DECODE_ME_SLC6A2_REGION_PLOT_37 = GwasLabRegionPlotTargetLocusTask.create(
    asset_id="decode_me_region_plot_slc6a2_locus",
    sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
    vcf_name_for_ld="1kg_eur_hg19",
    chrom=16,
    pos=55_714_809,
    buffer=250_000,
    skip=0,
    title="DecodeME GWAS-1 (ME/CFS): SLC6A2 (norepinephrine transporter) region",
)

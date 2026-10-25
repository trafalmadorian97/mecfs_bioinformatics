"""
Region plot of the Keaton et al. diastolic blood pressure GWAS around SLC6A2, the gene encoding the norepinephrine transporter (NET).

Centered on the midpoint of the SLC6A2 gene body (GRCh37 chr16:55,689,516-55,740,102).
"""

from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_standard_analysis import (
    KEATON_DBP_STANDARD_ANALYSIS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_region_plot_task_arbitrary_locus import (
    GwasLabRegionPlotTargetLocusTask,
)

_KEATON_DBP_SUMSTATS = KEATON_DBP_STANDARD_ANALYSIS.magma_tasks.sumstats_task
assert isinstance(_KEATON_DBP_SUMSTATS, GWASLabCreateSumstatsTask)

KEATON_DBP_SLC6A2_REGION_PLOT_37 = GwasLabRegionPlotTargetLocusTask.create(
    asset_id="keaton_dbp_region_plot_slc6a2_locus",
    sumstats_task=_KEATON_DBP_SUMSTATS,
    vcf_name_for_ld="1kg_eur_hg19",
    chrom=16,
    pos=55_714_809,
    buffer=250_000,
    title="Diastolic blood pressure (Keaton et al. 2024): SLC6A2 (norepinephrine transporter) region",
)

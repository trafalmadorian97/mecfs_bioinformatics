"""
Task to create a stacked plot illustrating the fine-mapping of the DecodeME Chromosome 1 locus
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.without_palindromes.susie_finemap_decode_me_37_chr1_173_locus import (
    DECODE_ME_GWAS_1_SUSIE_FINEMAP_CHR1_173_000_001_LOCUS,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.raw.magma_ensembl_gene_location_reference_data_build_37 import (
    MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    HeatmapOptions,
    RegionSelectDefault,
    SusieStackPlotTask,
)

DECODE_ME_GWAS_1_SUSIE_FINEMAP_CHR1_173_000_001_LOCUS_STACKPLOT = (
    SusieStackPlotTask.create(
        asset_id="stackplot_decode_me_gwas_1_37_susie_finemap_chr1_173000001_locus",
        susie_task=DECODE_ME_GWAS_1_SUSIE_FINEMAP_CHR1_173_000_001_LOCUS,
        gene_info_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
        gene_info_pipe=IdentityPipe(),
        region_mode=RegionSelectDefault(),
        heatmap_options=HeatmapOptions(
            heatmap_bin_options=None, mode="ld2", cmap="plasma"
        ),
    )
)

from mecfs_bio.asset_generator.genetic_correlation_asset_generator import (
    genetic_corr_by_ct_ldsc_asset_generator,
)
from mecfs_bio.assets.gwas.ldl.willer_et_al.analysis.willer_ldl_standard_analysis import (
    WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.myocardial_infarction.analysis.mi_standard_analysis import (
    MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
    QuantPhenotype,
    SumstatsSource,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_min_in_col import (
    FilterRowsByMinInCol,
)
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.build_system.task.pipes.to_polars_pipe import ToPolarsPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

MI_LDL_WILLER_CORRELATION = genetic_corr_by_ct_ldsc_asset_generator(
    base_name="willer_ldl_mi_corr",
    sources=[
        SumstatsSource(
            task=WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="willer_LDL",
            sample_info=QuantPhenotype(),
            pipe=SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN,
                constant=89138,  # source: data file
            ),
        ),
        SumstatsSource(
            MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="MI",
            sample_info=BinaryPhenotypeSampleInfo(
                sample_prevalence=0.09,  # 39,074/(39,074+392,979)
                estimated_population_prevalence=0.04,  # https://www.ncbi.nlm.nih.gov/books/NBK83160/
            ),
            pipe=CompositePipe(
                [
                    ComputeBetaIfNeededPipe(),
                    ComputeSEPipe(),
                    FilterRowsByMinInCol(1e-15, col=GWASLAB_SE_COL),
                    ToPolarsPipe(),
                    SetColToConstantPipe(GWASLAB_SAMPLE_SIZE_COLUMN, constant=432053),
                ]
            ),
        ),
    ],
)

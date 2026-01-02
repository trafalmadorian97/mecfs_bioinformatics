from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.raw.raw_sch_pgc2022 import (
    PGC_2022_SCH_RAW,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)

SCH_PGC_2022_STANDARD_ANALYSIS = concrete_standard_analysis_generator_assume_already_has_rsid(
    base_name="pgc_2022_sch",
    raw_gwas_data_task=PGC_2022_SCH_RAW,
    fmt=GWASLabColumnSpecifiers(
        rsid="ID",
        chrom="CHROM",
        pos="POS",
        ea="A1",  # since a comment says ##INFO=<ID=BETA,Number=1,Type=Float,Description="beta or ln(OR) of A1">
        nea="A2",
        info="IMPINFO",
        beta="BETA",
        se="SE",
        p="PVAL",
        n="NEFF",  # use effective sample siz
        snpid=None,
        OR=None,
    ),
    sample_size=58749,  # from summary statistics file
    include_master_gene_lists=False,
)

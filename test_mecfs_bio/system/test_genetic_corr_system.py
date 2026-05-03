"""
System test for `GeneticCorrelationByCTLDSCTask`.

Computes genetic correlation between two real GWAS summary statistic datasets
(Lee et al. 2018 educational attainment + Bentham et al. 2015 systemic lupus
erythematosus) using the production 1000 Genomes EUR LD-score reference.

These two studies were chosen because they are the smallest pair already wired
into the production CT-LDSC asset generator
(`mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_initial_asset_generator`)
and they exercise both phenotype-info paths (quantitative + binary).

This test is run weekly by `.github/workflows/system_test_genetic_corr_by_ct_ldsc.yml`.
implemented by Claude
"""

from pathlib import Path

import pandas as pd

from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.processed_gwas_data.lee_et_al_magma_task_generator import (
    LEE_ET_AL_2018_COMBINED_MAGMA_TASKS,
)
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.analysis_results.bentham_2015_standard_analysis import (
    BENTHAM_LUPUS_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.auxiliary.prevalence_info import (
    BENTHAM_LUPUS_PREVALENCE_INFO,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    GeneticCorrelationByCTLDSCTask,
    QuantPhenotype,
    SumstatsSource,
)
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN
from test_mecfs_bio.system.util import log_on_error

CT_LDSC_SYSTEM_TEST_TASK = GeneticCorrelationByCTLDSCTask.create(
    asset_id="ct_ldsc_system_test_lee_bentham",
    sources=[
        SumstatsSource(
            task=LEE_ET_AL_2018_COMBINED_MAGMA_TASKS.sumstats_task,
            alias="Educational_attainment",
            pipe=SetColToConstantPipe(GWASLAB_SAMPLE_SIZE_COLUMN, 257841),
            sample_info=QuantPhenotype(),
        ),
        SumstatsSource(
            task=BENTHAM_LUPUS_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="Lupus",
            pipe=SetColToConstantPipe(GWASLAB_SAMPLE_SIZE_COLUMN, 14267),
            sample_info=BENTHAM_LUPUS_PREVALENCE_INFO,
        ),
    ],
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    build="19",
)


def test_genetic_corr_by_ct_ldsc(tmp_path: Path):
    info_store = tmp_path / "info_store.yaml"
    asset_root = tmp_path / "asset_store"
    with log_on_error(info_store):
        asset_root.mkdir(parents=True, exist_ok=True)
        runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = runner.run([CT_LDSC_SYSTEM_TEST_TASK], incremental_save=True)
        assert result is not None

        rg_asset = result[CT_LDSC_SYSTEM_TEST_TASK.asset_id]
        assert isinstance(rg_asset, FileAsset)
        df = pd.read_csv(rg_asset.path)
        assert len(df) == 1
        for col in ("p1", "p2", "rg", "se", "z", "p"):
            assert col in df.columns, f"missing column {col} in {list(df.columns)}"
        assert df["rg"].notna().all()
        # LDSC point estimates can drift slightly outside [-1, 1] when h2 is small
        assert (df["rg"].abs() <= 1.5).all()

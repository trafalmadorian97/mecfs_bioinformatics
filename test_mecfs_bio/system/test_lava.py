"""
System test for LavaTask: local genetic correlation analysis using LAVA.

This test downloads reference data (g1000 EUR LD reference, locus definitions)
and two GWAS summary statistics datasets (DecodeME ME/CFS and Liu et al. IBD),
runs the LavaTask, and verifies the output files are produced.

NOTE: This test downloads large reference files and may take several minutes.
It is intended to run as a daily/weekly system test, not on every pull request.
"""

from pathlib import Path

import pandas as pd

from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_liftover_to_37_with_rsid import (
    LIU_ET_AL_2023_IBD_LIFTOVER_TO_37_WITH_RSID,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_build_37_with_rsid import (
    DECODE_ME_GWAS_1_LIFTOVER_TO_37_WITH_RSID,
)
from mecfs_bio.assets.reference_data.lava_ld_reference.g1000_eur.processed.lava_thousand_geomes_eur_ld_ref_extracted import (
    LAVA_G100_EUR_LD_REF_EXTRACTED,
)
from mecfs_bio.assets.reference_data.lava_locus_file.default.raw.default_lava_locus_file import (
    DEFAULT_LAVA_LOCUS_FILE,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.lava_task import (
    BIVAR_RESULTS_FILENAME,
    UNIV_RESULTS_FILENAME,
    LavaPhenotypeDataSource,
    LavaTask,
    LDReferenceInfo,
)
from test_mecfs_bio.system.util import log_on_error

LAVA_SYSTEM_TEST_TASK = LavaTask.create(
    asset_id="lava_system_test",
    sources=[
        LavaPhenotypeDataSource(
            task=DECODE_ME_GWAS_1_LIFTOVER_TO_37_WITH_RSID,
            alias="mecfs",
        ),
        LavaPhenotypeDataSource(
            task=LIU_ET_AL_2023_IBD_LIFTOVER_TO_37_WITH_RSID,
            alias="ibd",
        ),
    ],
    ld_reference_info=LDReferenceInfo(
        ld_ref_task=LAVA_G100_EUR_LD_REF_EXTRACTED,
        filename_prefix="g1000_eur",
    ),
    lava_locus_definitions_task=DEFAULT_LAVA_LOCUS_FILE,
)


def test_lava_local_genetic_correlation(tmp_path: Path):
    """
    Test that the LavaTask can run end-to-end with real data:
    - Downloads g1000 EUR LD reference
    - Downloads default locus definitions
    - Downloads and processes DecodeME ME/CFS and Liu et al. IBD GWAS data
    - Runs LAVA univariate heritability and bivariate genetic correlation
    - Verifies output CSV files are produced with expected columns
    """
    info_store = tmp_path / "info_store.yaml"
    asset_root = tmp_path / "asset_store"

    with log_on_error(info_store):
        asset_root.mkdir(parents=True, exist_ok=True)
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run(
            [LAVA_SYSTEM_TEST_TASK],
            incremental_save=True,
        )
        assert result is not None

        lava_output = result[LAVA_SYSTEM_TEST_TASK.asset_id]
        output_dir = lava_output.path

        # Verify univariate results file exists and has expected columns
        univ_path = output_dir / UNIV_RESULTS_FILENAME
        assert univ_path.exists(), f"Univariate results not found at {univ_path}"
        univ_df = pd.read_csv(univ_path)
        assert len(univ_df) > 0, "Univariate results should not be empty"
        assert "phen" in univ_df.columns
        assert "h2.obs" in univ_df.columns
        assert "p" in univ_df.columns
        assert "locus" in univ_df.columns
        assert set(univ_df["phen"].unique()) == {"mecfs", "ibd"}

        # Verify bivariate results file exists
        bivar_path = output_dir / BIVAR_RESULTS_FILENAME
        assert bivar_path.exists(), f"Bivariate results not found at {bivar_path}"
        bivar_df = pd.read_csv(bivar_path)
        if len(bivar_df) > 0:
            assert "phen1" in bivar_df.columns
            assert "phen2" in bivar_df.columns
            assert "rho" in bivar_df.columns
            assert "p" in bivar_df.columns
            assert "locus" in bivar_df.columns

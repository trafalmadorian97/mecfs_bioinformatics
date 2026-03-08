"""
System tests for LavaTask: local genetic correlation analysis using LAVA.


Implemented by Github copilot
Two system tests are provided:

1. test_lava_local_genetic_correlation:
   Uses full GWAS data (DecodeME ME/CFS and Liu et al. IBD) with g1000 EUR LD
   reference. Downloads large files and is intended for daily/weekly runs.

2. test_lava_vignette_data:
   Uses the small sample data bundled with the LAVA R package vignettes.
   Downloads only a few MB and runs much faster.
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
from mecfs_bio.assets.reference_data.lava_vignette_data.lava_vignette_data import (
    LAVA_VIGNETTE_BMI_SUMSTATS,
    LAVA_VIGNETTE_LD_REF,
    LAVA_VIGNETTE_LOCI,
    LAVA_VIGNETTE_NEURO_SUMSTATS,
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
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from test_mecfs_bio.system.util import log_on_error

# --- Full data test (slow) ---

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

# --- Vignette data test (fast) ---

# Pipes to rename vignette BMI columns to gwaslab format
# BMI columns: CHR BP SNPID_UKB A1 A2 MAF BETA SE STAT P NMISS INFO_UKB
_BMI_RENAME_PIPE = CompositePipe(
    pipes=[
        RenameColPipe(old_name="SNPID_UKB", new_name="rsID"),
        RenameColPipe(old_name="A1", new_name="EA"),
        RenameColPipe(old_name="A2", new_name="NEA"),
        RenameColPipe(old_name="NMISS", new_name="N"),
    ]
)

# Pipes to rename vignette neuro columns to gwaslab format
# Neuro columns: CHR POS SNP A1 A2 EAF BETA SE Z P N_analyzed INFO
_NEURO_RENAME_PIPE = CompositePipe(
    pipes=[
        RenameColPipe(old_name="SNP", new_name="rsID"),
        RenameColPipe(old_name="A1", new_name="EA"),
        RenameColPipe(old_name="A2", new_name="NEA"),
        RenameColPipe(old_name="N_analyzed", new_name="N"),
    ]
)

LAVA_VIGNETTE_TEST_TASK = LavaTask.create(
    asset_id="lava_vignette_test",
    sources=[
        LavaPhenotypeDataSource(
            task=LAVA_VIGNETTE_BMI_SUMSTATS,
            alias="bmi",
            pipe=_BMI_RENAME_PIPE,
        ),
        LavaPhenotypeDataSource(
            task=LAVA_VIGNETTE_NEURO_SUMSTATS,
            alias="neuro",
            pipe=_NEURO_RENAME_PIPE,
        ),
    ],
    ld_reference_info=LDReferenceInfo(
        ld_ref_task=LAVA_VIGNETTE_LD_REF,
        filename_prefix="g1000_test",
    ),
    lava_locus_definitions_task=LAVA_VIGNETTE_LOCI,
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


def test_lava_vignette_data(tmp_path: Path):
    """
    Fast system test using the small sample data from the LAVA vignettes.

    Uses the vignette's g1000_test LD reference (~3MB) and two continuous
    phenotypes (bmi, neuro) with 42 test loci.
    Column names in the vignette data are renamed to gwaslab format via pipes.

    This test verifies that:
    - LavaTask runs end-to-end with vignette data
    - Univariate results contain both phenotypes
    - Bivariate results are produced for loci with significant signal in both
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
            [LAVA_VIGNETTE_TEST_TASK],
            incremental_save=True,
        )
        assert result is not None

        lava_output = result[LAVA_VIGNETTE_TEST_TASK.asset_id]
        output_dir = lava_output.path

        # Verify univariate results
        univ_path = output_dir / UNIV_RESULTS_FILENAME
        assert univ_path.exists(), f"Univariate results not found at {univ_path}"
        univ_df = pd.read_csv(univ_path)
        assert len(univ_df) > 0, "Univariate results should not be empty"
        assert "phen" in univ_df.columns
        assert "h2.obs" in univ_df.columns
        assert "p" in univ_df.columns
        assert "locus" in univ_df.columns
        assert set(univ_df["phen"].unique()) == {"bmi", "neuro"}

        # Verify bivariate results
        bivar_path = output_dir / BIVAR_RESULTS_FILENAME
        assert bivar_path.exists(), f"Bivariate results not found at {bivar_path}"
        bivar_df = pd.read_csv(bivar_path)
        # The LAVA vignette data has strong signal for both bmi and neuro
        # at several loci, so we expect bivariate results
        assert len(bivar_df) > 0, "Bivariate results should not be empty"
        assert "phen1" in bivar_df.columns
        assert "phen2" in bivar_df.columns
        assert "rho" in bivar_df.columns
        assert "p" in bivar_df.columns
        assert "locus" in bivar_df.columns

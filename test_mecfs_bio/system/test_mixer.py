"""
Added by Claude

System tests for MixerTask: univariate and bivariate MiXeR analysis.

Uses the small hello-world data from themixer repository
(chr21 and chr22 only, ~35K SNPs). Downloads a few MB and runs in
under 2 minutes with fast-run arguments in the univariate case and under 7 minutes in the bivariate case.
"""

import json
from pathlib import Path

from mecfs_bio.assets.reference_data.mixer.hello_world.mixer_hello_world_data import (
    MIXER_HELLO_WORLD_PREPARED,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.mixer.bivariate_mixer_task import (
    MIXER_BIVARIATE_FIT_JSON_PATTERN,
    MIXER_BIVARIATE_TEST_JSON_PATTERN,
    BivariateMixerTask,
)
from mecfs_bio.build_system.task.mixer.mixer_task import (
    MixerTask,
    PreformattedMixerDataSource,
)
from mecfs_bio.build_system.task.mixer.mixer_univariate_combine import (
    COMBINED_FIT_FILENAME_PREFIX,
    COMBINED_TEST_FILENAME_PREFIX,
    MixerRunSource,
    MixerUnivariateCombine,
)
from mecfs_bio.build_system.task.mixer.mixer_univariate_results import (
    TEST_OUTPUT_PREFIX,
    MixerUnivariateSummarizeResultsTask,
)
from test_mecfs_bio.system.util import log_on_error


def test_mixer_univariate_hello_world(tmp_path: Path):
    """
    Fast system test: runs univariate fit1 + test1 on trait1 from the
    cross-trait MiXeR hello-world example.

    Uses the hello-world PLINK data (chr21-22) with generated .ld files.
    Skips the --extract flag for simplicity.
    Verifies that fit1 and test1 produce valid JSON output files.
    """
    info_store = tmp_path / "info_store.yaml"
    asset_root = tmp_path / "asset_store"

    mixer_task = MixerTask.create(
        asset_id="mixer_univariate_hello_world_test",
        trait_1_source=PreformattedMixerDataSource(
            task=MIXER_HELLO_WORLD_PREPARED,
            filename="trait1.sumstats.gz",
            alias="trait1",
        ),
        ref_data_directory_task=MIXER_HELLO_WORLD_PREPARED,
        ld_file_pattern="g1000_eur_hm3_chr@.ld",
        bim_file_pattern="g1000_eur_hm3_chr@.bim",
        extract_file_pattern_gen=None,
        chr_args="21-22",
        extra_args=(
            "--fit-sequence",
            "diffevo-fast",
            "neldermead-fast",
            "--diffevo-fast-repeats",
            "2",
            "--seed",
            "123",
        ),
        reps_to_perform=[1],
    )
    combine_task = MixerUnivariateCombine.create(
        asset_id="mixer_univariate_hello_world_test_combine",
        mixer_source_runs=[MixerRunSource(task=mixer_task, rep=1)],
        trait_name="dummy_trait",
    )
    plot_task = MixerUnivariateSummarizeResultsTask.create(
        asset_id="mixer_univariate_hello_world_test_plot",
        combine_task=combine_task,
        trait_name="dummy_trait",
    )

    with log_on_error(info_store):
        asset_root.mkdir(parents=True, exist_ok=True)
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run(
            [mixer_task, combine_task, plot_task],
            incremental_save=True,
        )
        assert result is not None

        output = result[mixer_task.asset_id]
        assert isinstance(output, DirectoryAsset)
        output_dir = output.path

        # Verify fit1 output
        fit_json_path = output_dir / "trait1.fit.1.json"
        assert fit_json_path.exists(), f"fit1 output not found at {fit_json_path}"
        with open(fit_json_path) as f:
            fit_data = json.load(f)
        assert "pi" in fit_data or "params" in fit_data, (
            f"fit1 JSON missing expected keys, got: {list(fit_data.keys())}"
        )
        # Verify test1 output
        test_json_path = output_dir / "trait1.test.1.json"
        assert test_json_path.exists(), f"test1 output not found at {test_json_path}"
        with open(test_json_path) as f:
            test_data = json.load(f)
        assert len(test_data) > 0, "test1 JSON should not be empty"

        combine_output = result[combine_task.asset_id]
        assert isinstance(combine_output, DirectoryAsset)
        combine_output_dir = combine_output.path
        combined_fit_file_path = combine_output_dir / (
            COMBINED_FIT_FILENAME_PREFIX + ".json"
        )
        combined_test_file_path = combine_output_dir / (
            COMBINED_TEST_FILENAME_PREFIX + ".json"
        )
        assert combined_fit_file_path.exists(), (
            f"combined fit file not found at {combined_fit_file_path}"
        )
        assert combined_test_file_path.exists(), (
            f"combined test file not found at {combined_test_file_path}"
        )

        plot_output = result[plot_task.asset_id]
        assert isinstance(plot_output, DirectoryAsset)
        plot_output_dir = plot_output.path
        assert (plot_output_dir / f"{TEST_OUTPUT_PREFIX}.power.png").exists()


FAST_RUN_EXTRA_ARGS = (
    "--fit-sequence",
    "diffevo-fast",
    "neldermead-fast",
    "--diffevo-fast-repeats",
    "2",
    "--seed",
    "123",
)


def _make_univariate_task(
    asset_id: str,
    trait_filename: str,
    trait_alias: str,
) -> MixerTask:
    """Helper: create a single-rep univariate MixerTask for a hello-world trait."""
    return MixerTask.create(
        asset_id=asset_id,
        trait_1_source=PreformattedMixerDataSource(
            task=MIXER_HELLO_WORLD_PREPARED,
            filename=trait_filename,
            alias=trait_alias,
        ),
        ref_data_directory_task=MIXER_HELLO_WORLD_PREPARED,
        ld_file_pattern="g1000_eur_hm3_chr@.ld",
        bim_file_pattern="g1000_eur_hm3_chr@.bim",
        extract_file_pattern_gen=None,
        chr_args="21-22",
        extra_args=FAST_RUN_EXTRA_ARGS,
        reps_to_perform=[1],
    )


def test_mixer_bivariate_hello_world(tmp_path: Path):
    """
    Fast system test: runs bivariate (cross-trait) MiXeR on trait1 and trait2
    from the hello-world example.

    Steps: univariate fit1+test1 for each trait, then bivariate fit2+test2.
    Verifies that fit2 and test2 produce valid JSON output files.
    """
    info_store = tmp_path / "info_store.yaml"
    asset_root = tmp_path / "asset_store"

    trait1_univariate = _make_univariate_task(
        asset_id="bivar_test_trait1_univariate",
        trait_filename="trait1.sumstats.gz",
        trait_alias="trait1",
    )
    trait2_univariate = _make_univariate_task(
        asset_id="bivar_test_trait2_univariate",
        trait_filename="trait2.sumstats.gz",
        trait_alias="trait2",
    )

    bivariate_task = BivariateMixerTask.create(
        asset_id="bivar_test_bivariate",
        trait_1_source=PreformattedMixerDataSource(
            task=MIXER_HELLO_WORLD_PREPARED,
            filename="trait1.sumstats.gz",
            alias="trait1",
        ),
        trait_2_source=PreformattedMixerDataSource(
            task=MIXER_HELLO_WORLD_PREPARED,
            filename="trait2.sumstats.gz",
            alias="trait2",
        ),
        ref_data_directory_task=MIXER_HELLO_WORLD_PREPARED,
        trait_1_univariate_task=trait1_univariate,
        trait_2_univariate_task=trait2_univariate,
        extract_file_pattern_gen=None,
        chr_args="21-22",
        extra_args=FAST_RUN_EXTRA_ARGS,
        ld_file_pattern="g1000_eur_hm3_chr@.ld",
        bim_file_pattern="g1000_eur_hm3_chr@.bim",
    )

    with log_on_error(info_store):
        asset_root.mkdir(parents=True, exist_ok=True)
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run(
            [trait1_univariate, trait2_univariate, bivariate_task],
            incremental_save=True,
        )
        assert result is not None

        # Verify bivariate fit2 output
        bivar_output = result[bivariate_task.asset_id]
        assert isinstance(bivar_output, DirectoryAsset)
        bivar_dir = bivar_output.path

        fit_json_name = MIXER_BIVARIATE_FIT_JSON_PATTERN.replace("@", "1")
        fit_json_path = bivar_dir / fit_json_name
        assert fit_json_path.exists(), f"fit2 output not found at {fit_json_path}"
        with open(fit_json_path) as f:
            fit_data = json.load(f)
        assert len(fit_data) > 0, "fit2 JSON should not be empty"

        # Verify bivariate test2 output
        test_json_name = MIXER_BIVARIATE_TEST_JSON_PATTERN.replace("@", "1")
        test_json_path = bivar_dir / test_json_name
        assert test_json_path.exists(), f"test2 output not found at {test_json_path}"
        with open(test_json_path) as f:
            test_data = json.load(f)
        assert len(test_data) > 0, "test2 JSON should not be empty"

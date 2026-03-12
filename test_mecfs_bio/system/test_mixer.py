"""
System tests for MixerTask: univariate and bivariate MiXeR analysis.

Uses the small hello-world data from the comorment/mixer repository
(chr21 and chr22 only, ~35K SNPs). Downloads a few MB and runs in
under 2 minutes with fast-run arguments.
"""

import json
from pathlib import Path

from mecfs_bio.assets.reference_data.mixer.hello_world.mixer_hello_world_data import (
    MIXER_HELLO_WORLD_PREPARED,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.mixer_bivariate_task import (
    MixerTask,
    PreformattedMixerDataSource,
    UnivariateMode,
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
        mixer_mode=UnivariateMode(),
        ce_data_directory_task=MIXER_HELLO_WORLD_PREPARED,
        ld_file_pattern="g1000_eur_hm3_chr@.ld",
        bim_file_pattern="g1000_eur_hm3_chr@.bim",
        extract_file_pattern_gen=None,
        chr_args="21-22",
        extra_args=(
            "--fit-sequence", "diffevo-fast", "neldermead-fast",
            "--diffevo-fast-repeats", "2",
            "--seed", "123",
        ),
        reps_to_perform=[1],
    )

    with log_on_error(info_store):
        asset_root.mkdir(parents=True, exist_ok=True)
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run(
            [mixer_task],
            incremental_save=True,
        )
        assert result is not None

        output = result[mixer_task.asset_id]
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

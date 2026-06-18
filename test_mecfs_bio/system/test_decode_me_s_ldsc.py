from pathlib import Path

import pandas as pd
import pytest

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.s_ldsc.decode_me_sldsc import (
    DECODE_ME_S_LDSC,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from test_mecfs_bio.system.util import log_on_error

# Expected cell-type S-LDSC coefficient p-values for the GTEx/Franke multi-tissue gene-expression
# dataset, from a clean from-scratch build on gwaslab 4.1.7 (CI run 26691178699, 2026-05-30).
# This guards the documented DecodeME S-LDSC results against changes -- e.g. a gwaslab upgrade --
# that would materially shift them: the gwaslab 3.6.x -> 4.1.7 upgrade changed the upstream
# liftover/harmonization and moved these p-values ~5-10%, which the old "did it crash?" test could
# not catch. If a dependency/pipeline change legitimately alters these, update both the values
# here AND docs/Analysis/ME_CFS/DecodeME/j_S-LDSC_DecodeME_Analysis.md.
EXPECTED_MULTITISSUE_PVALUES = {
    "A08.186.211.Brain": 2.334044e-07,
    "A08.186.211.730.885.287.500.571.735.Visual.Cortex": 1.507441e-06,
    "A08.186.211.464.405.Hippocampus": 2.035723e-06,
    "Brain_Cortex": 2.346594e-06,
}

# Relative tolerance for the p-value comparison. S-LDSC is deterministic given fixed inputs, so
# this is loose enough to absorb trivial floating-point noise but tight enough to flag the kind of
# material shift (~5-10%) that the gwaslab major-version upgrade produced.
P_VALUE_REL_TOL = 1e-2


def test_run_s_lsdc(tmp_path: Path):
    """
    Test that we can run stratified linkage disequilibrium score regression on the DECODE ME data,
    and that the headline cell-type p-values still match the documented results.
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
        # Emit the full cts result tables to the logs (disable pandas truncation), so every
        # tissue's p-value is recoverable from the CI logs rather than just the head/tail of the
        # large GTEx/Franke and Roadmap-chromatin tables. The cts task logs the result via an
        # f-string repr, which respects these process-wide display options.
        with pd.option_context(
            "display.max_rows",
            None,
            "display.max_columns",
            None,
            "display.width",
            1000,
            "display.max_colwidth",
            None,
        ):
            store = test_runner.run(list(DECODE_ME_S_LDSC.get_terminal_tasks()))
        assert store is not None
        _assert_documented_multitissue_pvalues(store)


def _assert_documented_multitissue_pvalues(store) -> None:
    """
    Check the GTEx/Franke cell-type coefficient p-values against the documented values.

    The cts result table is a transitive dependency of the S-LDSC terminal tasks, so it is present
    in the runner's returned asset store.
    """
    cell_task = DECODE_ME_S_LDSC.partitioned_tasks[
        "multi_tissue_gene_expression"
    ].cell_analysis_task
    asset = store[cell_task.asset_id]
    assert isinstance(asset, FileAsset)
    results = pd.read_csv(asset.path)
    pvalues = dict(zip(results["Name"], results["Coefficient_P_value"]))

    for name, expected in EXPECTED_MULTITISSUE_PVALUES.items():
        assert name in pvalues, f"tissue {name!r} missing from S-LDSC results"
        assert pvalues[name] == pytest.approx(expected, rel=P_VALUE_REL_TOL), (
            f"S-LDSC coefficient p-value for {name!r} is {pvalues[name]:.4e}, "
            f"expected ~{expected:.4e} (>{P_VALUE_REL_TOL:.0%} change). A dependency or pipeline "
            "change may have materially altered the documented DecodeME S-LDSC results; review, "
            "and if intended, update the expected values here and the documentation."
        )

    # The strongest CNS hit should remain the most significant tissue.
    assert min(pvalues, key=lambda tissue: pvalues[tissue]) == "A08.186.211.Brain"

import tempfile
from pathlib import Path

import pandas as pd

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_lead_variants import (
    DECODE_ME_GWAS_1_LEAD_VARIANTS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_manhattan import (
    DECODE_ME_GWAS_1_MANHATTAN_PLOT,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_manhattan_and_qq import (
    DECODE_ME_GWAS_1_MANHATTAN_AND_QQ_PLOT,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SNPID_COL

expected_vars = {
    "1:173846152:T:C",
    "6:26239176:A:G",
    "6:97984426:C:CA",
    "15:54866724:A:G",
    "17:52183006:C:T",
    "20:48914387:T:TA",
}


def test_run_initial_analysis():
    """
    Test that we can run the initial DecodeME analysis:
    - Download data
    - Filter lead variants
    - Produce Manhattan plots

    """
    with tempfile.TemporaryDirectory() as tempdirname:
        tempdir = Path(tempdirname)
        info_store = tempdir / "info_store.yaml"
        asset_root = tempdir / "asset_store"
        asset_root.mkdir(parents=True, exist_ok=True)
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run(
            [
                DECODE_ME_GWAS_1_MANHATTAN_AND_QQ_PLOT,
                DECODE_ME_GWAS_1_MANHATTAN_PLOT,
                DECODE_ME_GWAS_1_LEAD_VARIANTS,
            ]
        )
        variants_asset = result[DECODE_ME_GWAS_1_LEAD_VARIANTS.asset_id]
        assert isinstance(variants_asset, FileAsset)
        df: pd.DataFrame = pd.read_csv(variants_asset.path, index_col=None)
        lead_vars = set(df[GWASLAB_SNPID_COL].tolist())
        assert expected_vars == lead_vars

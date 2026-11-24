from pathlib import Path

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.task.consolidate_ld_scores_task import (
    LD_SCORE_CHROM_COL,
    LD_SCORE_LD_SCORE_COL,
    LD_SCORE_M_5_50_COL,
    LD_SCORE_RSID_COL,
)
from mecfs_bio.build_system.task.decode_database.build_decode_somascan_protein_task import (
    DECODE_SLIM_FILENAME,
    BuildDecodeSomascanProteinTask,
)
from mecfs_bio.build_system.task.decode_database.decode_somascan_protein import (
    DecodeSomascanProtein,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_h2 import batched_h2
from mecfs_bio.build_system.task.ppp_ldsc.decode_protein_heritability_task import (
    DecodeProteinHeritabilityTask,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_ENTREZ_ID_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.build_system.wf.s3_access import S3Endpoint, S3ObjectRef
from mecfs_bio.constants.csf_database_constants import EntrezGeneId, SeqId, UniProtId
from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_H2_N_P05_COL,
    DECODE_H2_N_P95_COL,
    DECODE_H2_PROTEIN_ID_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.ppp_database_constants import PPP_INDEX_IS_STRAND_AMBIGUOUS_COL
from mecfs_bio.constants.ppp_ldsc_constants import (
    PPP_H2_H2_COL,
    PPP_H2_N_BAR_COL,
    PPP_H2_VARIANT_SET_COL,
    PPP_VARIANT_SET_ALL,
    PPP_VARIANT_SET_CIS_EXCLUDED,
)

_N_SNPS = 3_000
_M = 100_000.0
_N_CENTER = 33_000.0
_H2 = 0.3
_ENDPOINT = S3Endpoint(url="https://fake.example", region="r", credentials_profile="p")


def _parquet_meta(asset_id: str) -> ResultTableMeta:
    return ResultTableMeta(
        id=AssetId(asset_id),
        trait="t",
        project="p",
        sub_dir=Path("x"),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )


def _protein_task(
    index_task: FakeTask, seq_id: str, entrez_id: str | None, symbol: str
) -> BuildDecodeSomascanProteinTask:
    protein = DecodeSomascanProtein(
        seq_id=SeqId(seq_id),
        gene_symbol=symbol,
        target_name=symbol,
        uniprot=UniProtId("P00000"),
        entrez_gene_id=None if entrez_id is None else EntrezGeneId(entrez_id),
        ref=S3ObjectRef(endpoint=_ENDPOINT, bucket="b", key=f"{seq_id}.txt.gz"),
        size=1,
        md5="0" * 32,
    )
    return BuildDecodeSomascanProteinTask.create(
        index_task=index_task,
        protein=protein,
        asset_id=f"protein_{seq_id}",
        index_name="test",
        group="smp",
    )


def _slim(rng: np.random.Generator, ld: np.ndarray) -> pl.DataFrame:
    """SE = 1 and BETA = Z with Z^2 following the LDSC expectation; N varies by +-8%."""
    n = _N_CENTER * (1 + rng.uniform(-0.08, 0.08, size=_N_SNPS))
    chi2 = 1 + _N_CENTER * _H2 / _M * ld + rng.normal(0, 0.1, size=_N_SNPS)
    return pl.DataFrame(
        {
            GWASLAB_BETA_COL: np.sqrt(np.maximum(chi2, 0.0)),
            GWASLAB_SE_COL: np.ones(_N_SNPS),
            GWASLAB_SAMPLE_SIZE_COLUMN: n,
        }
    ).cast(pl.Float32)


def test_decode_h2_uses_median_n_and_matches_ppp_kernel(tmp_path: Path) -> None:
    rng = np.random.default_rng(0)
    pos = 10_000 * np.arange(1, _N_SNPS + 1)
    rsids = [f"rs{i}" for i in range(_N_SNPS)]
    ld = rng.uniform(1.0, 50.0, size=_N_SNPS)
    index_task = FakeTask(_parquet_meta("index"))
    ld_task = FakeTask(_parquet_meta("ld"))
    magma_task = FakeTask(_parquet_meta("magma"))
    with_coords = _protein_task(index_task, "1-1", "348", "APOE")
    without_coords = _protein_task(index_task, "2-1", None, "NOTAGENE")
    assets: dict[AssetId, Asset] = {}

    def put(asset_id: AssetId, frame: pl.DataFrame) -> None:
        path = tmp_path / f"{asset_id}.parquet"
        frame.write_parquet(path)
        assets[asset_id] = FileAsset(path)

    put(
        index_task.asset_id,
        pl.DataFrame(
            {
                GWASLAB_CHROM_COL: np.ones(_N_SNPS, dtype=np.int64),
                GWASLAB_POS_COL: pos,
                GWASLAB_RSID_COL: rsids,
                PPP_INDEX_IS_STRAND_AMBIGUOUS_COL: np.zeros(_N_SNPS, dtype=bool),
            }
        ),
    )
    put(
        ld_task.asset_id,
        pl.DataFrame(
            {
                LD_SCORE_RSID_COL: rsids,
                LD_SCORE_CHROM_COL: np.ones(_N_SNPS, dtype=np.int64),
                LD_SCORE_LD_SCORE_COL: ld,
                LD_SCORE_M_5_50_COL: np.full(_N_SNPS, _M),
            }
        ),
    )
    put(
        magma_task.asset_id,
        pl.DataFrame(
            {
                GENE_INFO_ENTREZ_ID_COL: [348],
                GENE_INFO_CHROM_COL: ["1"],
                GENE_INFO_START_COL: [15_000_000],
                GENE_INFO_END_COL: [15_010_000],
                GENE_INFO_STRAND_COL: ["+"],
                GENE_INFO_NAME_COL: ["APOE"],
            }
        ),
    )
    slims = {}
    for task in (with_coords, without_coords):
        directory = tmp_path / task.asset_id
        directory.mkdir()
        slims[task.asset_id] = _slim(rng, ld)
        slims[task.asset_id].write_parquet(directory / DECODE_SLIM_FILENAME)
        assets[task.asset_id] = DirectoryAsset(directory)

    task = DecodeProteinHeritabilityTask.create(
        asset_id="decode_h2",
        protein_tasks=(with_coords, without_coords),
        index_task=index_task,
        ld_scores_task=ld_task,
        gene_coords_task=magma_task,
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(
        scratch, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assert isinstance(result, FileAsset)
    table = pl.read_parquet(result.path)

    for protein_task in (with_coords, without_coords):
        slim = slims[protein_task.asset_id]
        n = slim[GWASLAB_SAMPLE_SIZE_COLUMN].to_numpy().astype(float)
        chi2 = (
            slim[GWASLAB_BETA_COL].to_numpy().astype(float)
            / slim[GWASLAB_SE_COL].to_numpy().astype(float)
        ) ** 2
        expected = batched_h2(chi2[:, None], ld, np.array([np.median(n)]), _M)
        rows = table.filter(
            pl.col(DECODE_H2_PROTEIN_ID_COL) == str(protein_task.protein.seq_id)
        )
        all_row = rows.filter(pl.col(PPP_H2_VARIANT_SET_COL) == PPP_VARIANT_SET_ALL)
        assert all_row[PPP_H2_N_BAR_COL].item() == pytest.approx(np.median(n))
        assert all_row[DECODE_H2_N_P05_COL].item() == pytest.approx(
            np.quantile(n, 0.05)
        )
        assert all_row[DECODE_H2_N_P95_COL].item() == pytest.approx(
            np.quantile(n, 0.95)
        )
        assert all_row[PPP_H2_H2_COL].item() == pytest.approx(expected.h2[0], rel=1e-9)

    variant_sets = dict(
        table.group_by(DECODE_H2_PROTEIN_ID_COL)
        .agg(pl.col(PPP_H2_VARIANT_SET_COL).sort())
        .rows()
    )
    assert variant_sets == {
        "1-1": [PPP_VARIANT_SET_ALL, PPP_VARIANT_SET_CIS_EXCLUDED],
        "2-1": [PPP_VARIANT_SET_ALL],
    }

from pathlib import Path, PurePath

import numpy as np
import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.protein_coloc.coloc_display_table_task import (
    DISPLAY_SCHEMA,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_locus_stack_plot_task import (
    STACK_PLOT_INDEX_FILENAME,
    ColocLocusStackPlotTask,
    _LocusSpec,
    _trait_track,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.genetic_map_constants import (
    GMAP_CM_COL,
    GMAP_POS_COL,
    GMAP_RATE_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    DISPLAY_CIS_STATUS_COL,
    DISPLAY_GENE_COL,
    DISPLAY_LOCUS_COL,
    DISPLAY_PP_H4_COL,
    DISPLAY_PROTEIN_COL,
    LOCUS_ID_COL,
    N_PROTEINS_COL,
    PLOT_PNG_FILENAME_COL,
    PLOT_SVG_FILENAME_COL,
)
from test_mecfs_bio.unit.build_system.task.protein_coloc.coloc_result_fixtures import (
    COLOC_ID,
    coloc_task,
    write_coloc_dir,
)

_CHR1_LOCUS = "chr1:1000-3000"
_CHR2_LOCUS = "chr2:1000-3000"

# (gene, protein, locus, pp_h4): tested pairs of the fixture coloc dir, so their windows exist.
_DISPLAY_ROWS = [
    ("AAA", "OID1", _CHR1_LOCUS, 0.95),
    ("CCC", "OID3", _CHR2_LOCUS, 0.90),
    ("BBB", "OID2", _CHR1_LOCUS, 0.85),
]


def _table_task(asset_id: str) -> FakeTask:
    return FakeTask(
        ResultTableMeta(
            id=AssetId(asset_id),
            trait="synthetic",
            project="test",
            sub_dir=PurePath("analysis"),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
    )


def _write_display(path: Path, rows: list[tuple[str, str, str, float]]) -> None:
    pl.DataFrame(
        {
            DISPLAY_GENE_COL: [r[0] for r in rows],
            DISPLAY_PROTEIN_COL: [r[1] for r in rows],
            DISPLAY_LOCUS_COL: [r[2] for r in rows],
            DISPLAY_CIS_STATUS_COL: ["trans"] * len(rows),
            DISPLAY_PP_H4_COL: [r[3] for r in rows],
        },
        schema={
            DISPLAY_GENE_COL: pl.String,
            DISPLAY_PROTEIN_COL: pl.String,
            DISPLAY_LOCUS_COL: pl.String,
            DISPLAY_CIS_STATUS_COL: pl.String,
            DISPLAY_PP_H4_COL: pl.Float64,
        },
    ).select(
        [
            pl.col(c)
            if c
            in (
                DISPLAY_GENE_COL,
                DISPLAY_PROTEIN_COL,
                DISPLAY_LOCUS_COL,
                DISPLAY_CIS_STATUS_COL,
                DISPLAY_PP_H4_COL,
            )
            else pl.lit(None).alias(c)
            for c in DISPLAY_SCHEMA.names()
        ]
    ).cast(DISPLAY_SCHEMA).write_parquet(path)


def _write_trait(path: Path) -> None:
    rows = []
    for chrom in (1, 2):
        pos = np.arange(900, 3101, 10)
        z = np.zeros(len(pos))
        z[len(pos) // 2] = 8.0
        rows.append(
            pl.DataFrame(
                {
                    GWASLAB_CHROM_COL: chrom,
                    GWASLAB_POS_COL: pos,
                    GWASLAB_EFFECT_ALLELE_COL: "A",
                    GWASLAB_NON_EFFECT_ALLELE_COL: "G",
                    GWASLAB_BETA_COL: z * 0.01,
                    GWASLAB_SE_COL: 0.01,
                }
            )
        )
    pl.concat(rows).write_parquet(path)


def _write_genetic_map(path: Path) -> None:
    pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [1, 1, 2, 2],
            GMAP_POS_COL: [1000, 3000, 1000, 3000],
            GMAP_RATE_COL: [0.5, 2.0, 1.0, 1.0],
            GMAP_CM_COL: [0.0, 0.1, 0.0, 0.1],
        }
    ).write_parquet(path)


def _write_genes(path: Path) -> None:
    pl.DataFrame(
        {
            GENE_INFO_CHROM_COL: ["1", "2"],
            GENE_INFO_START_COL: [1500, 1500],
            GENE_INFO_END_COL: [2500, 2500],
            GENE_INFO_STRAND_COL: ["+", "-"],
            GENE_INFO_NAME_COL: ["AAA", "CCC"],
        }
    ).write_parquet(path)


def _stack_plots(
    tmp_path: Path, display_rows: list[tuple[str, str, str, float]], max_loci: int
) -> Path:
    coloc_dir = write_coloc_dir(tmp_path / "coloc")
    paths = {
        name: tmp_path / f"{name}.parquet"
        for name in ("display", "trait", "gmap", "genes")
    }
    _write_display(paths["display"], display_rows)
    _write_trait(paths["trait"])
    _write_genetic_map(paths["gmap"])
    _write_genes(paths["genes"])
    task = ColocLocusStackPlotTask.create(
        asset_id="stack_plots",
        coloc_task=coloc_task(),
        trait_task=_table_task("trait"),
        display_task=_table_task("display"),
        gene_info_task=_table_task("genes"),
        genetic_map_task=_table_task("gmap"),
        max_loci=max_loci,
    )
    assets: dict[AssetId, Asset] = {
        AssetId(COLOC_ID): DirectoryAsset(coloc_dir),
        **{AssetId(name): FileAsset(path) for name, path in paths.items()},
    }
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(
        scratch, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assert isinstance(result, DirectoryAsset)
    return result.path


def test_one_figure_per_display_locus(tmp_path: Path) -> None:
    out = _stack_plots(tmp_path, _DISPLAY_ROWS, max_loci=10)
    index = pl.read_parquet(out / STACK_PLOT_INDEX_FILENAME)
    assert index[LOCUS_ID_COL].to_list() == [_CHR1_LOCUS, _CHR2_LOCUS]
    assert index[N_PROTEINS_COL].to_list() == [2, 1]
    for filename in (
        index[PLOT_PNG_FILENAME_COL].to_list() + index[PLOT_SVG_FILENAME_COL].to_list()
    ):
        assert (out / filename).is_file()


def test_max_loci_truncates_by_best_pp_h4(tmp_path: Path) -> None:
    out = _stack_plots(tmp_path, _DISPLAY_ROWS, max_loci=1)
    index = pl.read_parquet(out / STACK_PLOT_INDEX_FILENAME)
    assert index[LOCUS_ID_COL].to_list() == [_CHR1_LOCUS]
    assert len(list(out.glob("*.png"))) == 1


def test_empty_display_table_gives_empty_index(tmp_path: Path) -> None:
    out = _stack_plots(tmp_path, [], max_loci=10)
    index = pl.read_parquet(out / STACK_PLOT_INDEX_FILENAME)
    assert index.height == 0
    assert list(out.glob("*.png")) == []


def test_trait_track_drops_unusable_variants() -> None:
    """An SE = 0 row (infinite z) or a non-finite beta must not become the plotted trait lead."""
    trait = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [1, 1, 1, 1],
            GWASLAB_POS_COL: [1100, 1500, 2000, 2500],
            GWASLAB_BETA_COL: [0.5, float("nan"), 0.08, 0.01],
            GWASLAB_SE_COL: [0.0, 0.01, 0.01, 0.01],
        }
    ).lazy()
    spec = _LocusSpec(locus_id=_CHR1_LOCUS, chrom=1, start=1000, end=3000, proteins=())
    track = _trait_track(trait, spec)
    assert track.pos.tolist() == [2000, 2500]

"""A hand-built TraitProteinColocTask result directory for the display and plot tests."""

from pathlib import Path

import numpy as np
import polars as pl

from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    BETA_PROTEIN_COL,
    BETA_TRAIT_COL,
    CIS_STATUS_COL,
    GENE_COL,
    LOCUS_CHR_COL,
    LOCUS_END_COL,
    LOCUS_ID_COL,
    LOCUS_LEAD_EA_COL,
    LOCUS_LEAD_NEA_COL,
    LOCUS_LEAD_POS_COL,
    LOCUS_START_COL,
    PP_H4_COL,
    PROTEIN_ID_COL,
    RAW_RESULTS_FILENAME,
    SE_PROTEIN_COL,
    SE_TRAIT_COL,
    STATUS_COL,
    STATUS_EXCLUDED_MHC,
    STATUS_NOT_SELECTED_BUDGET,
    STATUS_TESTED,
    WINDOWS_DIRNAME,
    pair_file_stem,
)

COLOC_ID = "coloc"

# (LOCUS_ID, PROTEIN_ID, GENE, PP_H4, STATUS)
RESULT_ROWS = [
    ("chr1:1000-3000", "OID1", "AAA", 0.95, STATUS_TESTED),
    ("chr1:1000-3000", "OID2", "BBB", 0.85, STATUS_TESTED),
    ("chr2:1000-3000", "OID3", "CCC", 0.30, STATUS_TESTED),
    ("chr2:1000-3000", "OID4", "DDD", None, STATUS_NOT_SELECTED_BUDGET),
    ("chr6:1000-3000", None, None, None, STATUS_EXCLUDED_MHC),
]


def coloc_task() -> FakeTask:
    return FakeTask(ResultDirectoryMeta(id=COLOC_ID, trait="synthetic", project="test"))


def write_coloc_dir(directory: Path) -> Path:
    windows = directory / WINDOWS_DIRNAME
    windows.mkdir(parents=True)
    rows = []
    for locus_id, protein_id, gene, pp_h4, status in RESULT_ROWS:
        chrom = int(locus_id.split(":")[0].removeprefix("chr"))
        rows.append(
            {
                LOCUS_ID_COL: locus_id,
                LOCUS_CHR_COL: chrom,
                LOCUS_START_COL: 1000,
                LOCUS_END_COL: 3000,
                LOCUS_LEAD_POS_COL: 2000,
                LOCUS_LEAD_EA_COL: "A",
                LOCUS_LEAD_NEA_COL: "G",
                PROTEIN_ID_COL: protein_id,
                GENE_COL: gene,
                CIS_STATUS_COL: "cis",
                PP_H4_COL: pp_h4,
                STATUS_COL: status,
            }
        )
        if status == STATUS_TESTED:
            assert protein_id is not None
            _window(chrom).write_parquet(
                windows / f"{pair_file_stem(locus_id, protein_id)}.parquet"
            )
    pl.DataFrame(rows).cast(
        {
            LOCUS_CHR_COL: pl.Int32,
            LOCUS_START_COL: pl.Int32,
            LOCUS_END_COL: pl.Int32,
            LOCUS_LEAD_POS_COL: pl.Int32,
        }
    ).write_parquet(directory / RAW_RESULTS_FILENAME)
    return directory


def _window(chrom: int) -> pl.DataFrame:
    rng = np.random.default_rng(chrom)
    pos = np.arange(1000, 3001, 20)
    z = rng.normal(0, 1, len(pos))
    z[len(pos) // 2] = 10.0
    return pl.DataFrame(
        {
            GWASLAB_CHROM_COL: chrom,
            GWASLAB_POS_COL: pos,
            GWASLAB_EFFECT_ALLELE_COL: "A",
            GWASLAB_NON_EFFECT_ALLELE_COL: "G",
            BETA_TRAIT_COL: z * 0.01,
            SE_TRAIT_COL: 0.01,
            BETA_PROTEIN_COL: z * 0.02,
            SE_PROTEIN_COL: 0.02,
        }
    )

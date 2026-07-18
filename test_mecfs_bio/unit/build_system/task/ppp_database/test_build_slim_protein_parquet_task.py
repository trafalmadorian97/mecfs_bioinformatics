import gzip
import math
from pathlib import Path

import polars as pl
import pyarrow.parquet as pq
import pytest

from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    PppProteinFile,
    align_protein_to_index,
    aligned_columns,
    write_slim_aligned_parquet,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

_REGENIE_HEADER = (
    "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA"
)


def test_ppp_protein_file_tar_filename():
    # The tar name is derived from the structured fields, not stored; panels may contain an
    # underscore (Inflammation_II) and it must round-trip the manifest's filename column.
    protein = PppProteinFile(
        gene="A1BG",
        uniprot="P04217",
        oid="OID30771",
        version="v1",
        panel="Inflammation_II",
        synid="syn52359450",
    )
    assert protein.tar_filename == "A1BG_P04217_OID30771_v1_Inflammation_II.tar"


def _write_chr_gz(path: Path, rows: list[str]) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(_REGENIE_HEADER + "\n")
        handle.write("\n".join(rows) + "\n")


def test_align_protein_to_index():
    # Index deliberately NOT position-sorted, to prove output follows index row order.
    index = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [2, 1, 1],
            GWASLAB_POS_COL: [300, 100, 200],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["T", "G", "T"],
            GWASLAB_EFFECT_ALLELE_COL: ["A", "A", "C"],
        }
    )
    protein = pl.DataFrame(
        {
            "CHROM": [1, 1, 5],
            "GENPOS": [100, 200, 999],
            "ALLELE0": ["G", "C", "A"],
            "ALLELE1": ["A", "T", "G"],
            "BETA": [0.5, 0.3, 0.9],
            "SE": [0.1, 0.2, 0.4],
        }
    )

    out = align_protein_to_index(index, protein, include_sample_size=False)
    assert out.columns == aligned_columns(False)
    beta = out[GWASLAB_BETA_COL].to_list()
    se = out[GWASLAB_SE_COL].to_list()
    assert math.isnan(beta[0]) and math.isnan(se[0])  # chr2:300 absent
    assert beta[1] == pytest.approx(0.5, abs=1e-6)  # same orientation
    assert beta[2] == pytest.approx(-0.3, abs=1e-6)  # swapped -> flipped


def test_align_protein_to_index_with_sample_size():
    # N is allele-independent: present variants carry the (constant) protein N, absent
    # variants NaN, regardless of orientation.
    index = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [2, 1, 1],
            GWASLAB_POS_COL: [300, 100, 200],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["T", "G", "T"],
            GWASLAB_EFFECT_ALLELE_COL: ["A", "A", "C"],
        }
    )
    protein = pl.DataFrame(
        {
            "CHROM": [1, 1, 5],
            "GENPOS": [100, 200, 999],
            "ALLELE0": ["G", "C", "A"],
            "ALLELE1": ["A", "T", "G"],
            "BETA": [0.5, 0.3, 0.9],
            "SE": [0.1, 0.2, 0.4],
            "N": [42000, 42000, 42000],
        }
    )

    out = align_protein_to_index(index, protein, include_sample_size=True)
    assert out.columns == aligned_columns(True)
    assert GWASLAB_SAMPLE_SIZE_COLUMN in out.columns
    n = out[GWASLAB_SAMPLE_SIZE_COLUMN].to_list()
    assert math.isnan(n[0])  # chr2:300 absent
    assert n[1] == pytest.approx(42000.0)  # same orientation
    assert n[2] == pytest.approx(42000.0)  # flipped orientation, N unaffected


def test_write_slim_aligned_parquet(tmp_path: Path):
    extracted = tmp_path / "extracted" / "PROTEIN"
    extracted.mkdir(parents=True)
    # chr1: two variants (one same orientation, one swapped); chr2: one variant.
    _write_chr_gz(
        extracted / "discovery_chr1_x.gz",
        [
            "1 100 1:100:G:A:imp:v1 G A 0.2 0.9 100 ADD 0.5 0.1 1 1 NA",
            "1 200 1:200:C:T:imp:v1 C T 0.3 0.9 100 ADD 0.3 0.2 1 1 NA",
        ],
    )
    _write_chr_gz(
        extracted / "discovery_chr2_x.gz",
        ["2 300 2:300:T:A:imp:v1 T A 0.1 0.9 100 ADD 0.7 0.3 1 1 NA"],
    )
    # A chr1 variant the protein does not carry -> a normal, expected NaN.
    index = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [1, 1, 1, 2],
            GWASLAB_POS_COL: [100, 200, 250, 300],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G", "T", "C", "T"],
            GWASLAB_EFFECT_ALLELE_COL: ["A", "C", "G", "A"],
        }
    )

    out_path = tmp_path / "slim.parquet"
    write_slim_aligned_parquet(
        tmp_path / "extracted", index, out_path, include_sample_size=False
    )

    out = pl.read_parquet(out_path)
    assert out.columns == aligned_columns(False)
    assert out.height == 4
    beta = out[GWASLAB_BETA_COL].to_list()
    se = out[GWASLAB_SE_COL].to_list()
    assert beta[0] == pytest.approx(0.5, abs=1e-6)
    assert beta[1] == pytest.approx(-0.3, abs=1e-6)  # swapped
    assert math.isnan(beta[2]) and math.isnan(se[2])  # variant absent from protein
    assert beta[3] == pytest.approx(0.7, abs=1e-6)
    assert [se[0], se[1], se[3]] == pytest.approx([0.1, 0.2, 0.3], abs=1e-6)

    # The output must actually use byte-stream-split (the baked-in format).
    encodings = pq.ParquetFile(out_path).metadata.row_group(0).column(0).encodings
    assert "BYTE_STREAM_SPLIT" in encodings


def test_write_slim_aligned_parquet_finds_x_chromosome_file(tmp_path: Path):
    # PPP names the X file with the letter (discovery_chrX_...) but codes the
    # chromosome numerically inside it (CHROM=23), and the index follows the numeric
    # coding. Globbing for chr23 would therefore find nothing and, before this was
    # handled, would have failed every protein once X entered the index.
    extracted = tmp_path / "extracted" / "PROTEIN"
    extracted.mkdir(parents=True)
    _write_chr_gz(
        extracted / "discovery_chrX_PROTEIN.gz",
        ["23 2781514 X:2699555:C:A:imp:v1 C A 0.4 1.2 33084 ADD -0.5 0.2 1 1 NA"],
    )
    index = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [23],
            GWASLAB_POS_COL: [2781514],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["C"],
            GWASLAB_EFFECT_ALLELE_COL: ["A"],
        }
    )

    out_path = tmp_path / "slim.parquet"
    write_slim_aligned_parquet(
        tmp_path / "extracted", index, out_path, include_sample_size=False
    )

    out = pl.read_parquet(out_path)
    assert out.height == 1
    assert out[GWASLAB_BETA_COL][0] == pytest.approx(-0.5, abs=1e-6)
    assert out[GWASLAB_SE_COL][0] == pytest.approx(0.2, abs=1e-6)


def test_write_slim_aligned_parquet_with_sample_size(tmp_path: Path):
    # With include_sample_size, the slim file gains a constant-per-protein N column
    # (float32, NaN where the protein lacks a variant) and still uses byte-stream-split.
    extracted = tmp_path / "extracted" / "PROTEIN"
    extracted.mkdir(parents=True)
    _write_chr_gz(
        extracted / "discovery_chr1_x.gz",
        [
            "1 100 1:100:G:A:imp:v1 G A 0.2 0.9 42000 ADD 0.5 0.1 1 1 NA",
            "1 200 1:200:C:T:imp:v1 C T 0.3 0.9 42000 ADD 0.3 0.2 1 1 NA",
        ],
    )
    index = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [1, 1, 1],
            GWASLAB_POS_COL: [100, 200, 250],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G", "T", "C"],
            GWASLAB_EFFECT_ALLELE_COL: ["A", "C", "G"],
        }
    )

    out_path = tmp_path / "slim.parquet"
    write_slim_aligned_parquet(
        tmp_path / "extracted", index, out_path, include_sample_size=True
    )

    out = pl.read_parquet(out_path)
    assert out.columns == aligned_columns(True)
    n = out[GWASLAB_SAMPLE_SIZE_COLUMN].to_list()
    assert n[0] == pytest.approx(42000.0) and n[1] == pytest.approx(42000.0)
    assert math.isnan(n[2])  # index variant absent from protein
    # N must not break the baked-in byte-stream-split format.
    encodings = pq.ParquetFile(out_path).metadata.row_group(0).column(2).encodings
    assert "BYTE_STREAM_SPLIT" in encodings


def test_write_slim_aligned_parquet_index_not_chromosome_sorted_raises(tmp_path: Path):
    # Processing a chromosome at a time only reproduces index row order because the
    # index is chromosome-sorted. Nothing about an index makes that true by
    # construction, so if the index sort order ever changes, this must fail loudly
    # rather than silently emit beta/se against the wrong variants -- row order is the
    # contract every per-protein file is written against.
    extracted = tmp_path / "extracted" / "PROTEIN"
    extracted.mkdir(parents=True)
    _write_chr_gz(
        extracted / "discovery_chr1_x.gz",
        ["1 100 1:100:G:A:imp:v1 G A 0.2 0.9 100 ADD 0.5 0.1 1 1 NA"],
    )
    _write_chr_gz(
        extracted / "discovery_chr2_x.gz",
        ["2 300 2:300:T:A:imp:v1 T A 0.1 0.9 100 ADD 0.7 0.3 1 1 NA"],
    )
    # Sorted by position first, so the chromosomes interleave.
    index = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [2, 1],
            GWASLAB_POS_COL: [300, 100],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["T", "G"],
            GWASLAB_EFFECT_ALLELE_COL: ["A", "A"],
        }
    )

    with pytest.raises(AssertionError, match="chromosome-sorted"):
        write_slim_aligned_parquet(
            tmp_path / "extracted",
            index,
            tmp_path / "slim.parquet",
            include_sample_size=False,
        )


def test_write_slim_aligned_parquet_missing_chromosome_file_raises(tmp_path: Path):
    # An index chromosome with no regenie file must fail loudly rather than write a
    # whole block of NaNs: the download is far too expensive to discover afterwards
    # that the output is silently useless.
    extracted = tmp_path / "extracted" / "PROTEIN"
    extracted.mkdir(parents=True)
    _write_chr_gz(
        extracted / "discovery_chr1_x.gz",
        ["1 100 1:100:G:A:imp:v1 G A 0.2 0.9 100 ADD 0.5 0.1 1 1 NA"],
    )
    # chr2 is in the index but has no file (e.g. a chromosome-naming mismatch).
    index = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [1, 2],
            GWASLAB_POS_COL: [100, 300],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G", "T"],
            GWASLAB_EFFECT_ALLELE_COL: ["A", "A"],
        }
    )

    with pytest.raises(AssertionError, match="chr2"):
        write_slim_aligned_parquet(
            tmp_path / "extracted",
            index,
            tmp_path / "slim.parquet",
            include_sample_size=False,
        )

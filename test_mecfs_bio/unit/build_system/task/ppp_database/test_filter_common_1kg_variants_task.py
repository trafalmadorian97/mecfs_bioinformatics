from pathlib import Path

import polars as pl
import pytest

from mecfs_bio.build_system.task.ppp_database.filter_common_1kg_variants_task import (
    FILTERED_COLUMNS,
    ONEKG_AF_COL,
    ONEKG_ALT_COL,
    ONEKG_REF_COL,
    filter_vcf_to_frame,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL

# Minimal VCF: a common SNV, an indel, a rare SNV, another common SNV.
_VCF = (
    "##fileformat=VCFv4.2\n"
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
    "##contig=<ID=chr1>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tG\t.\t.\tAF=0.2\n"  # common SNV -> kept
    "chr1\t200\t.\tA\tAT\t.\t.\tAF=0.3\n"  # indel -> dropped (-v snps)
    "chr1\t300\t.\tC\tT\t.\t.\tAF=0.001\n"  # rare SNV -> dropped (MAF)
    "chr1\t400\t.\tG\tA\t.\t.\tAF=0.5\n"  # common SNV -> kept
)


def test_filter_vcf_to_frame(tmp_path: Path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)

    out = filter_vcf_to_frame(vcf, maf_threshold=0.01, work_dir=tmp_path)

    assert out.columns == FILTERED_COLUMNS
    assert out.height == 2
    assert out[GWASLAB_CHROM_COL].to_list() == [1, 1]
    assert out[GWASLAB_POS_COL].to_list() == [100, 400]
    assert out[ONEKG_REF_COL].to_list() == ["A", "G"]
    assert out[ONEKG_ALT_COL].to_list() == ["G", "A"]
    assert out[ONEKG_AF_COL].to_list() == pytest.approx([0.2, 0.5], abs=1e-6)
    assert out.schema[GWASLAB_CHROM_COL] == pl.Int32

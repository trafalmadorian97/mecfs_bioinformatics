"""
Unit tests for `GeneticCorrelationByCTLDSCTask` and its helpers.

These tests focus on the pure / lightweight components of
`mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task`
that we can exercise without running LDSC end-to-end (which would require
a 22-chromosome LD-score reference and genome-wide HapMap3 coverage).

For end-to-end coverage of `GeneticCorrelationByCTLDSCTask.execute`, see
`test_mecfs_bio/system/test_genetic_corr_system.py`.

Implemented by Claude
"""

from pathlib import Path

import gwaslab as gl
import pandas as pd
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    FilterSettings,
    QuantPhenotype,
    SumstatsSource,
    filter_sumstats,
    get_compatible_snps_polars,
    load_and_preprocess_sumstats,
)


@frozen
class HapmapSNP:
    """A single HapMap3 variant. Hardcoded entries below are real HapMap3 SNPs."""

    rsid: str
    a1: str
    a2: str
    chrom: int
    pos: int


# A small subset of real HapMap3 SNPs (build hg38), copied verbatim from the
# HapMap3 reference list. HapMap3 is a frozen reference panel that will not
# change, so it is safe to hardcode these here rather than read them from
# gwaslab's internal directory layout (which can change between versions).
#
# All entries are non-palindromic SNPs (no A/T or C/G pairs) and non-indels,
# so they survive `filter_sumstats(FilterSettings())`.
_HAPMAP3_SAMPLE: list[HapmapSNP] = [
    HapmapSNP(rsid="rs3094315", a1="G", a2="A", chrom=1, pos=817186),
    HapmapSNP(rsid="rs3131972", a1="A", a2="G", chrom=1, pos=817341),
    HapmapSNP(rsid="rs3131969", a1="A", a2="G", chrom=1, pos=818802),
    HapmapSNP(rsid="rs1048488", a1="C", a2="T", chrom=1, pos=825532),
    HapmapSNP(rsid="rs3115850", a1="T", a2="C", chrom=1, pos=825767),
    HapmapSNP(rsid="rs2286139", a1="C", a2="T", chrom=1, pos=826352),
    HapmapSNP(rsid="rs12562034", a1="A", a2="G", chrom=1, pos=833068),
    HapmapSNP(rsid="rs4040617", a1="G", a2="A", chrom=1, pos=843942),
    HapmapSNP(rsid="rs2980300", a1="T", a2="C", chrom=1, pos=850609),
    HapmapSNP(rsid="rs2519031", a1="A", a2="G", chrom=1, pos=854250),
]


def _hapmap3_snps(n: int) -> list[HapmapSNP]:
    assert n <= len(_HAPMAP3_SAMPLE), (
        f"only {len(_HAPMAP3_SAMPLE)} hardcoded HapMap3 SNPs available"
    )
    return _HAPMAP3_SAMPLE[:n]


def _make_sumstats(
    rsids: list[str],
    chroms: list[int],
    positions: list[int],
    eas: list[str],
    neas: list[str],
    study: str = "trait1",
) -> gl.Sumstats:
    n = len(rsids)
    df = pd.DataFrame(
        {
            "rsID": rsids,
            "CHR": chroms,
            "POS": positions,
            "EA": eas,
            "NEA": neas,
            "BETA": [0.01] * n,
            "SE": [0.05] * n,
            "P": [0.5] * n,
            "N": [10000] * n,
        }
    )
    return gl.Sumstats(
        df,
        rsid="rsID",
        chrom="CHR",
        pos="POS",
        ea="EA",
        nea="NEA",
        beta="BETA",
        se="SE",
        p="P",
        n="N",
        study=study,
        build="38",
        verbose=False,
    )


def _sumstats_from_snps(snps: list[HapmapSNP], study: str = "trait1") -> gl.Sumstats:
    return _make_sumstats(
        rsids=[s.rsid for s in snps],
        chroms=[s.chrom for s in snps],
        positions=[s.pos for s in snps],
        eas=[s.a1 for s in snps],
        neas=[s.a2 for s in snps],
        study=study,
    )


def test_filter_sumstats_keeps_only_hapmap3():
    real = _hapmap3_snps(5)
    fake = [
        HapmapSNP(rsid="rs999999991", a1="A", a2="G", chrom=1, pos=10_000_001),
        HapmapSNP(rsid="rs999999992", a1="A", a2="G", chrom=1, pos=10_000_002),
        HapmapSNP(rsid="rs999999993", a1="A", a2="G", chrom=1, pos=10_000_003),
    ]
    sumstats = _sumstats_from_snps(real + fake)

    filter_sumstats(sumstats, FilterSettings(), build="38")

    assert set(sumstats.data["rsID"]) == {s.rsid for s in real}


def test_filter_sumstats_removes_indels_and_palindromes():
    snps = _hapmap3_snps(4)

    # Row 0: leave as a normal SNP (A/G)
    # Row 1: indel (AT/A)
    # Row 2: palindromic (A/T)
    # Row 3: leave as a normal SNP (different alleles to row 0)
    sumstats = _make_sumstats(
        rsids=[s.rsid for s in snps],
        chroms=[s.chrom for s in snps],
        positions=[s.pos for s in snps],
        eas=["A", "AT", "A", "C"],
        neas=["G", "A", "T", "T"],
    )

    filter_sumstats(
        sumstats,
        FilterSettings(keep_only_hapmap=False),
        build="38",
    )

    surviving = set(sumstats.data["rsID"])
    assert surviving == {snps[0].rsid, snps[3].rsid}


def test_get_compatible_snps_polars():
    # i: rs1 A/G, rs2 A/G,        rs3 C/G,         rs4 A/T,             rs5 A/G
    # j: rs1 A/G (match),
    #    rs2 G/A (flipped match),
    #    rs3 G/C (other-strand flipped match — complement-reversed of i's "C/G" is "C/G", flipped gives "G/C"),
    #    rs4 A/C (incompatible),
    #    rs5 missing entirely (still produces no row in the inner-join)
    df_i = pd.DataFrame(
        {
            "rsID": ["rs1", "rs2", "rs3", "rs4", "rs5"],
            "EA": ["A", "A", "C", "A", "A"],
            "NEA": ["G", "G", "G", "T", "G"],
        }
    )
    df_j = pd.DataFrame(
        {
            "rsID": ["rs1", "rs2", "rs3", "rs4"],
            "EA": ["A", "G", "G", "A"],
            "NEA": ["G", "A", "C", "C"],
        }
    )

    result = get_compatible_snps_polars(df_i, df_j)

    assert set(result["rsID"]) == {"rs1", "rs2", "rs3"}


def test_load_and_preprocess_sumstats(tmp_path: Path):
    real = _hapmap3_snps(5)
    fake = [HapmapSNP(rsid="rs999999991", a1="A", a2="G", chrom=1, pos=10_000_001)]
    sumstats = _sumstats_from_snps(real + fake, study="trait_a")

    pickle_path = tmp_path / "sumstats.pickle"
    gl.dump_pickle(sumstats, path=pickle_path)

    source_id = AssetId("trait_a_sumstats")
    fake_task = FakeTask(
        meta=GWASLabSumStatsMeta(
            id=source_id, trait="dummy_trait", project="dummy_project"
        )
    )

    quant = QuantPhenotype()
    source = SumstatsSource(task=fake_task, alias="trait_a", sample_info=quant)

    def fetch(asset_id: AssetId) -> Asset:
        assert asset_id == source_id
        return FileAsset(pickle_path)

    out_sumstats, _, _ = load_and_preprocess_sumstats(
        source=source, fetch=fetch, settings=FilterSettings(), build="38"
    )

    assert set(out_sumstats.data["rsID"]) == {s.rsid for s in real}

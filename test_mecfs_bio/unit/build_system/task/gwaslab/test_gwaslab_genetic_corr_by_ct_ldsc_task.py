"""
Unit tests for `GeneticCorrelationByCTLDSCTask` and its helpers.

Most tests here focus on the pure / lightweight components of
`mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task`.

`test_execute_runs_ldsc_and_returns_numeric_rg` additionally drives
`GeneticCorrelationByCTLDSCTask.execute` all the way through gwaslab's vendored
cross-trait LDSC code on tiny synthetic data, using a hand-built single-chromosome
LD-score reference and a rigged positive-heritability fit. It is a fast guard for
the gwaslab / numpy-version incompatibility that otherwise only the (weekly)
`test_mecfs_bio/system/test_genetic_corr_system.py` would catch.

Implemented by Claude
"""

import gzip
import math
from pathlib import Path

import gwaslab as gl
import pandas as pd
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    FilterSettings,
    GeneticCorrelationByCTLDSCTask,
    QuantPhenotype,
    SumstatsSource,
    filter_sumstats,
    get_compatible_snps_polars,
    load_and_preprocess_sumstats,
)
from mecfs_bio.build_system.wf.base_wf import make_wf


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
    ses: list[float] | None = None,
    betas: list[float] | None = None,
    ns: list[int] | None = None,
) -> gl.Sumstats:
    n = len(rsids)
    df = pd.DataFrame(
        {
            "rsID": rsids,
            "CHR": chroms,
            "POS": positions,
            "EA": eas,
            "NEA": neas,
            "BETA": [0.01] * n if betas is None else betas,
            "SE": [0.05] * n if ses is None else ses,
            "P": [0.5] * n,
            "N": [10000] * n if ns is None else ns,
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
    gl.dump_pickle(sumstats, path=str(pickle_path))

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


def test_load_and_preprocess_sumstats_drops_degenerate_z(tmp_path: Path):
    """Variants with SE == 0 give an infinite LDSC Z (BETA/SE) and abort the SVD
    inside the cross-trait regression, so they must not reach LDSC.

    Harmonised GWAS Catalog files do contain such variants: the Kerrebijn
    fibromyalgia sumstats report BETA as the smallest normal double with SE == 0.
    """
    snps = _hapmap3_snps(4)
    degenerate = {snps[1].rsid, snps[2].rsid}
    sumstats = _make_sumstats(
        rsids=[s.rsid for s in snps],
        chroms=[s.chrom for s in snps],
        positions=[s.pos for s in snps],
        eas=[s.a1 for s in snps],
        neas=[s.a2 for s in snps],
        study="trait_a",
        ses=[0.05, 0.0, 0.0, 0.05],
    )

    pickle_path = tmp_path / "sumstats.pickle"
    gl.dump_pickle(sumstats, path=str(pickle_path))

    source_id = AssetId("trait_a_sumstats")
    fake_task = FakeTask(
        meta=GWASLabSumStatsMeta(
            id=source_id, trait="dummy_trait", project="dummy_project"
        )
    )
    source = SumstatsSource(
        task=fake_task, alias="trait_a", sample_info=QuantPhenotype()
    )

    def fetch(asset_id: AssetId) -> Asset:
        assert asset_id == source_id
        return FileAsset(pickle_path)

    out_sumstats, _, _ = load_and_preprocess_sumstats(
        source=source, fetch=fetch, settings=FilterSettings(), build="38"
    )

    assert set(out_sumstats.data["rsID"]) == {s.rsid for s in snps} - degenerate


# --- End-to-end LDSC regression test -----------------------------------------
#
# This is the one test that drives GeneticCorrelationByCTLDSCTask.execute all the
# way through gwaslab's vendored cross-trait LDSC code,
#
# reaching the right lines requires gwaslab to estimate a POSITIVE heritability for
# both traits (otherwise it takes the _negative_hsq branch and never calls
# float()). We arrange that by rigging Z^2 = 1 + (N * h2 / M) * L2 with h2 > 0,
# so the LDSC regression of chi-square on LD score has a positive slope.

_LDSC_N = 20_000
_LDSC_M_5_50 = 1_000_000
_LDSC_H2 = 0.4
# Distinct LD score per SNP so the regression has spread.
_LDSC_L2 = [40.0 + 12.0 * i for i in range(len(_HAPMAP3_SAMPLE))]


def _rigged_betas(l2s: list[float]) -> list[float]:
    # With SE == 1, BETA == Z, so BETA = sqrt(E[Z^2]) gives a clean positive-h2 fit.
    return [math.sqrt(1.0 + (_LDSC_N * _LDSC_H2 / _LDSC_M_5_50) * l2) for l2 in l2s]


def _ldsc_sumstats(snps: list[HapmapSNP], l2s: list[float], study: str) -> gl.Sumstats:
    return _make_sumstats(
        rsids=[s.rsid for s in snps],
        chroms=[s.chrom for s in snps],
        positions=[s.pos for s in snps],
        eas=[s.a1 for s in snps],
        neas=[s.a2 for s in snps],
        study=study,
        betas=_rigged_betas(l2s),
        ses=[1.0] * len(snps),
        ns=[_LDSC_N] * len(snps),
    )


def _write_ld_reference(
    directory: Path, snps: list[HapmapSNP], l2s: list[float]
) -> None:
    """Write a synthetic single-chromosome LDSC reference (all SNPs are chr1).

    The reference is padded with SNPs absent from the sumstats so that it differs
    in length from them. gwaslab's smart_merge has a fast path (for identically
    shaped/ordered frames) that calls the removed pandas API `df.drop('SNP', 1)`;
    real references never trigger it because they dwarf the sumstats, so the
    padding keeps this test on the same pd.merge path as production.
    """
    rows = [(s.chrom, s.rsid, s.pos, l2) for s, l2 in zip(snps, l2s)]
    for i in range(5):
        rows.append((1, f"rs_pad{i}", 900_000 + i, 55.0 + i))
    with gzip.open(directory / "LDscore.1.l2.ldscore.gz", "wt") as f:
        f.write("CHR\tSNP\tBP\tL2\n")
        for chrom, snp, bp, l2 in rows:
            f.write(f"{chrom}\t{snp}\t{bp}\t{l2}\n")
    (directory / "LDscore.1.l2.M_5_50").write_text(f"{_LDSC_M_5_50}\n")
    (directory / "LDscore.1.l2.M").write_text(f"{_LDSC_M_5_50}\n")


def test_execute_runs_ldsc_and_returns_numeric_rg(tmp_path: Path):
    snps = _hapmap3_snps(len(_HAPMAP3_SAMPLE))

    ld_dir = tmp_path / "ld_ref"
    ld_dir.mkdir()
    _write_ld_reference(ld_dir, snps, _LDSC_L2)

    # trait2 is stored in reversed SNP order so the two sumstats differ in SNP
    # sequence at merge time; this, too, keeps gwaslab off the buggy smart_merge
    # fast path (which fires only when the frames match exactly).
    trait1 = _ldsc_sumstats(snps, _LDSC_L2, study="trait1")
    trait2 = _ldsc_sumstats(snps[::-1], _LDSC_L2[::-1], study="trait2")

    ids = {}
    sources = []
    for name, sumstats in (("trait1", trait1), ("trait2", trait2)):
        pickle_path = tmp_path / f"{name}.pickle"
        gl.dump_pickle(sumstats, path=str(pickle_path))
        source_id = AssetId(f"{name}_sumstats")
        ids[source_id] = FileAsset(pickle_path)
        fake_task = FakeTask(
            meta=GWASLabSumStatsMeta(id=source_id, trait=name, project="test")
        )
        sources.append(
            SumstatsSource(task=fake_task, alias=name, sample_info=QuantPhenotype())
        )

    ld_ref_id = AssetId("ld_ref")
    ld_ref_task = FakeTask(meta=SimpleDirectoryMeta(id=ld_ref_id))

    task = GeneticCorrelationByCTLDSCTask.create(
        asset_id="ct_ldsc_synthetic",
        sources=sources,
        ld_ref_task=ld_ref_task,
        build="38",
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == ld_ref_id:
            return DirectoryAsset(ld_dir)
        return ids[asset_id]

    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()
    result = task.execute(scratch_dir=scratch_dir, fetch=fetch, wf=make_wf())

    assert isinstance(result, FileAsset)
    df = pd.read_csv(result.path)
    assert len(df) == 1
    # A numeric (non-NA) rg proves gwaslab reached RG.__init__'s float() lines;
    # under numpy >= 2.4 that raises and this value would be NA / the run crashes.
    assert df["rg"].notna().all()
    assert df["rg"].abs().iloc[0] <= 1.5

from pathlib import Path

import polars as pl

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr1_174_128_548 import (
    POLYFUN_EXPLAIN_CHR1_174,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    COMBINED_CS_FILENAME,
    PIP_COLUMN,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)
from test_mecfs_bio.system.util import log_on_error

# The run config whose explainability assets this test builds. The L=10 run is the
# one whose PolyFun credible set is characterized in detail in the documentation.
_L10_LABEL = "l10"

# The documented top PolyFun SUSIE variant at the DecodeME chr1:173.5M-174.5M locus,
# from docs/Analysis/ME_CFS/DecodeME/Fine_Mapping/SUSIE-PolyFun_(External Prior)/
# a_Polyfun_Chr1_173M_174M_Locus.md. POS is hg19 (the Broad UKBB LD panel's build);
# note the gwaslab SNPID keeps the original hg38 position, so the lead variant is
# identified here by chromosome, hg19 position, and allele pair rather than SNPID.
# This guards the documented lead variant against a dependency or pipeline change
# that would move it. If such a change legitimately relocates the lead variant,
# update both the values here AND the documentation page above.
_EXPECTED_LEAD_CHROM = 1
_EXPECTED_LEAD_POS = 173_855_298
_EXPECTED_LEAD_ALLELES = frozenset({"A", "T"})


def test_decode_me_polyfun_susie_lead_variant(tmp_path: Path):
    """
    Test that we can run the DecodeME chr1 PolyFun-SUSIE explainability workflow for the L=10 run
    config -- the full uniform/polyfun SUSIE pair plus its contrast, explainability plot, and
    display tables -- and that the top PolyFun credible-set variant still matches the documented
    lead variant.
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
        group = POLYFUN_EXPLAIN_CHR1_174.groups_by_label[_L10_LABEL]
        # Build every explainability asset hanging off the L=10 run config. These transitively
        # build the UKBB LD download, the harmonized sumstats, and both SUSIE runs (uniform +
        # polyfun), so running them exercises the whole L=10 fine-mapping path end to end.
        store = test_runner.run(
            [
                group.plot_png,
                group.plot_svg,
                group.top_line_table,
                group.detailed_table,
                group.per_variant_annotation_table,
            ]
        )
        assert store is not None
        _assert_documented_lead_variant(store, group.susie_polyfun)


def _assert_documented_lead_variant(store, polyfun_susie_task) -> None:
    """
    Check the highest-PIP variant of the L=10 PolyFun SUSIE run against the documented lead variant.

    The PolyFun SUSIE run's combined credible-set table is a transitive dependency of the L=10
    explainability assets, so it is present in the runner's returned asset store.
    """
    asset = store[polyfun_susie_task.asset_id]
    assert isinstance(asset, DirectoryAsset)
    combined_cs = pl.read_parquet(asset.path / COMBINED_CS_FILENAME)
    top = combined_cs.sort(PIP_COLUMN, descending=True).row(0, named=True)

    chrom = int(top[GWASLAB_CHROM_COL])
    pos = int(top[GWASLAB_POS_COL])
    alleles = frozenset(
        {top[GWASLAB_EFFECT_ALLELE_COL], top[GWASLAB_NON_EFFECT_ALLELE_COL]}
    )

    assert (chrom, pos, alleles) == (
        _EXPECTED_LEAD_CHROM,
        _EXPECTED_LEAD_POS,
        _EXPECTED_LEAD_ALLELES,
    ), (
        f"Top PolyFun SUSIE variant is chr{chrom}:{pos} {sorted(alleles)}, expected "
        f"chr{_EXPECTED_LEAD_CHROM}:{_EXPECTED_LEAD_POS} {sorted(_EXPECTED_LEAD_ALLELES)}. "
        "A dependency or pipeline change may have moved the documented DecodeME chr1 PolyFun "
        "lead variant; review, and if intended, update the expected values here and the "
        "documentation."
    )

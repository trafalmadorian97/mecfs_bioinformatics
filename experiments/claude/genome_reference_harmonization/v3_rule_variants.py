"""
V3 follow-up: how candidate tightenings of the stringent ambiguous-indel rules change
kept and wrong counts on build-38 DecodeME.

Variants, each applied before decide_ambiguous_indels:
- baseline: the panel as built;
- drop_monomorphic: panel records whose EUR AF is 0 or 1 are treated as absent;
- drop_rare_records_<t>: panel records with min(AF, 1 - AF) < t are treated as absent;
- require_both_records: rows lacking either reading's record are dropped.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.v3_rule_variants \
    2>&1 | tee experiments/claude/genome_reference_harmonization/v3_rule_variants.log
"""

import attrs
import polars as pl

from experiments.claude.genome_reference_harmonization.tune_ambiguous_indel_rules import (
    BASE_OPTIONS,
    BUILD_38_TABLE,
    GRID_DISTANCES,
    GRID_MARGINS,
    ChromosomeAmbiguousIndels,
    _ambiguous_indels,
)
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.ambiguous_indels import (
    INDEL_ACTION_COL,
    decide_ambiguous_indels,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    chromosomes_to_harmonize,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    ACTION_KEEP,
    ACTION_SWAP,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

RARE_RECORD_THRESHOLDS = (0.001, 0.01)
VARIANT_COL = "variant"
DISTANCE_COL = "distance"
MARGIN_COL = "margin"
KEPT_COL = "kept"
WRONG_COL = "wrong"
_HAS_KEEP_COL = "_has_keep"
_HAS_FLIP_COL = "_has_flip"


def _drop_records(panel: pl.DataFrame, min_maf: float) -> pl.DataFrame:
    af = pl.col(PANEL_AF_COL)
    return panel.filter(pl.min_horizontal(af, 1 - af) >= min_maf)


def _drop_monomorphic(panel: pl.DataFrame) -> pl.DataFrame:
    af = pl.col(PANEL_AF_COL)
    return panel.filter((af > 0) & (af < 1))


def _both_records_mask(item: ChromosomeAmbiguousIndels) -> pl.Series:
    keys = [GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL]
    keep_keys = item.panel.select(
        GWASLAB_POS_COL,
        pl.col(PANEL_REF_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
        pl.col(PANEL_ALT_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
    ).with_columns(pl.lit(True).alias(_HAS_KEEP_COL))
    flip_keys = item.panel.select(
        GWASLAB_POS_COL,
        pl.col(PANEL_REF_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.col(PANEL_ALT_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
    ).with_columns(pl.lit(True).alias(_HAS_FLIP_COL))
    joined = (
        item.rows.select(keys)
        .join(keep_keys, on=keys, how="left", maintain_order="left")
        .join(flip_keys, on=keys, how="left", maintain_order="left")
    )
    return joined[_HAS_KEEP_COL].fill_null(False) & joined[_HAS_FLIP_COL].fill_null(
        False
    )


def main() -> None:
    assets = DEFAULT_RUNNER.run(
        [
            BUILD_38_TABLE,
            UCSC_HG38_INDEXED_FASTA,
            THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
        ]
    )
    fasta_asset = assets[UCSC_HG38_INDEXED_FASTA.asset_id]
    panel_asset = assets[THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset) and isinstance(
        panel_asset, FileAsset
    )
    fasta = IndexedFasta.open(fasta_asset.path)
    sumstats = scan_sumstats_as_polars(
        assets[BUILD_38_TABLE.asset_id], BUILD_38_TABLE.meta, IdentityPipe()
    )
    per_chromosome = [
        _ambiguous_indels(sumstats, chrom, fasta, panel_asset.path)
        for chrom in chromosomes_to_harmonize(sumstats, fasta, BASE_OPTIONS)
    ]
    variants: dict[str, list[ChromosomeAmbiguousIndels]] = {
        "baseline": per_chromosome,
        "drop_monomorphic": [
            attrs.evolve(item, panel=_drop_monomorphic(item.panel))
            for item in per_chromosome
        ],
        **{
            f"drop_rare_records_{threshold}": [
                attrs.evolve(item, panel=_drop_records(item.panel, threshold))
                for item in per_chromosome
            ]
            for threshold in RARE_RECORD_THRESHOLDS
        },
        "require_both_records": [
            ChromosomeAmbiguousIndels(
                rows=item.rows.filter(_both_records_mask(item)), panel=item.panel
            )
            for item in per_chromosome
        ],
    }
    results = []
    for name, items in variants.items():
        for distance in GRID_DISTANCES:
            for margin in GRID_MARGINS:
                options = attrs.evolve(
                    BASE_OPTIONS,
                    indel_max_af_distance=distance,
                    indel_min_af_margin=margin,
                )
                decisions = pl.concat(
                    [
                        decide_ambiguous_indels(item.rows, item.panel, options)
                        for item in items
                    ]
                )
                results.append(
                    {
                        VARIANT_COL: name,
                        DISTANCE_COL: distance,
                        MARGIN_COL: margin,
                        KEPT_COL: int(
                            (decisions[INDEL_ACTION_COL] == ACTION_KEEP).sum()
                        ),
                        WRONG_COL: int(
                            (decisions[INDEL_ACTION_COL] == ACTION_SWAP).sum()
                        ),
                    }
                )
    table = pl.DataFrame(results)
    with pl.Config(tbl_rows=200, tbl_cols=-1, tbl_width_chars=200):
        print(table)
        print("\nbest kept per variant among cells with the fewest wrong choices:")
        print(
            table.sort(
                [VARIANT_COL, WRONG_COL, KEPT_COL], descending=[False, False, True]
            )
            .group_by(VARIANT_COL, maintain_order=True)
            .head(3)
        )


if __name__ == "__main__":
    main()

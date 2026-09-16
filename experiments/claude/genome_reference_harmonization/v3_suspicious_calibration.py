"""
Revised V3, step 3: calibrate max_suspicious_indel_fraction under the production trust code.

For each testable trusted dataset (DecodeME build 38 as the truth set, plus the three GWAS
Catalog harmonised files), run the production count_trust_evidence_genome_wide against the
hg38 FASTA and 1000 Genomes EUR hg38 panel, at the recommended untrusted-path resolution
(indel_max_af_distance 0.02, indel_min_af_margin 0.3). This exercises the monomorphic-panel
filter and the checkable denominator, so the reported suspicious fraction is exactly what
decide_trust sees. We want DecodeME trusted and the harmonised files untrusted at
max_suspicious_indel_fraction 1e-4.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.v3_suspicious_calibration \
    2>&1 | tee experiments/claude/genome_reference_harmonization/v3_suspicious_calibration.log
"""

import tempfile
from pathlib import Path

import polars as pl

from experiments.claude.genome_reference_harmonization.tune_ambiguous_indel_rules import (
    BUILD_38_TABLE,
)
from experiments.claude.genome_reference_harmonization.v3_other_truth_sets import (
    Build,
    _load_table,
    _open,
)
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    chromosomes_to_harmonize,
    count_trust_evidence_genome_wide,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.trust import decide_trust
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL

# The recommended untrusted-path resolution (V3 grid: fewest wrong swaps). The suspicious
# gate reuses these, and max_suspicious_indel_fraction / min_checkable_ambiguous_indels keep
# their proposed defaults so this run confirms them.
CALIB_OPTIONS = GenomeReferenceHarmonizationOptions(
    indel_max_af_distance=0.02, indel_min_af_margin=0.3
)
_HG38_PICKLE_ROOT = Path("assets/base_asset_store/gwas")
_PICKLES = {
    "MVP myocardial infarction": _HG38_PICKLE_ROOT
    / "myocardial_infaction/million_veterans/gwaslab_sumstats/million_veterans_mi_raw_sumstats_task_no_liftover.pickle",
    "Bellenguez Alzheimer's": _HG38_PICKLE_ROOT
    / "alzheimers/bellenguez_et_al/gwaslab_sumstats/bellenguez_et_al_alz_raw_sumstats_task_no_liftover.pickle",
    "Kerrebijn fibromyalgia": _HG38_PICKLE_ROOT
    / "fibromyalgia/kerrebijn_et_al/gwaslab_sumstats/kerrebijin_fibro_raw_sumstats_task_no_liftover.pickle",
}


def _report(label: str, sumstats: pl.LazyFrame, fasta: IndexedFasta, panel_path: Path) -> None:
    sumstats = sumstats.filter(pl.col(GWASLAB_CHROM_COL).is_in(list(fasta.entries)))
    chromosomes = chromosomes_to_harmonize(sumstats, fasta, CALIB_OPTIONS)
    evidence = count_trust_evidence_genome_wide(
        sumstats, chromosomes, fasta, panel_path, CALIB_OPTIONS
    )
    trusted = decide_trust(evidence, CALIB_OPTIONS)
    suspicious = evidence.suspicious
    print(f"\n######## {label}")
    print(f"  counts: {evidence.counts}")
    print(
        f"  suspicious={suspicious.suspicious:,} checkable={suspicious.checkable:,} "
        f"fraction={suspicious.fraction:.3e}"
    )
    print(
        f"  decide_trust (max_suspicious_indel_fraction={CALIB_OPTIONS.max_suspicious_indel_fraction:.0e}, "
        f"min_checkable={CALIB_OPTIONS.min_checkable_ambiguous_indels}): trusted={trusted}"
    )


def main() -> None:
    opened = _open(
        Build(
            "hg38",
            UCSC_HG38_INDEXED_FASTA,
            THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
        )
    )
    fasta, panel_path = opened.fasta, opened.panel_path
    decode_asset = DEFAULT_RUNNER.run([BUILD_38_TABLE])[BUILD_38_TABLE.asset_id]
    decode_sumstats = scan_sumstats_as_polars(
        decode_asset, BUILD_38_TABLE.meta, IdentityPipe()
    )
    _report("DecodeME build 38 (truth set)", decode_sumstats, fasta, panel_path)
    with tempfile.TemporaryDirectory() as scratch:
        for label, pickle_path in _PICKLES.items():
            parquet_path = Path(scratch) / "table.parquet"
            has_eaf = _load_table(pickle_path, parquet_path)
            if not has_eaf:
                print(f"\n######## {label}\n  no EAF column: skipped")
                continue
            _report(label, pl.scan_parquet(parquet_path), fasta, panel_path)


if __name__ == "__main__":
    main()

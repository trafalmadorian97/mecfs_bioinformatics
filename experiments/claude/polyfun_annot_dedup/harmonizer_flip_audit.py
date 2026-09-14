"""Does HarmonizeGWASWithReferenceViaAlleles flip the allele order of DecodeME variants?

For each DecodeME polyfun-explain locus, run the CURRENT harmonizer Task code (the cached outputs
could predate code changes, since the build cache ignores code) against the locally cached inputs,
then match each output row back to its input row by SNPID and classify it:
  - kept: output EA/NEA == input EA/NEA
  - flipped: output EA/NEA == input NEA/EA
For flipped rows, also check that BETA was negated, EAF became 1 - EAF, whether the cases/controls
allele-frequency columns were (not) flipped, and which orientation matches the hg19 FASTA.
Finally compare the freshly computed output with the cached asset.

Run:
  pixi r python experiments/claude/polyfun_annot_dedup/harmonizer_flip_audit.py \
    2>&1 | tee experiments/claude/polyfun_annot_dedup/harmonizer_flip_audit.log
"""

import tempfile
from pathlib import Path

import numpy as np
import polars as pl

from experiments.claude.polyfun_annot_dedup.benchmark_fasta_ref_check import HG19, read_fai, ref_match
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr1_174_128_548 import (
    POLYFUN_EXPLAIN_CHR1_174,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr6_26_215_000 import (
    POLYFUN_EXPLAIN_CHR6_26,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr6_97_505_620 import (
    POLYFUN_EXPLAIN_CHR6_97,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr15_54_925_638 import (
    POLYFUN_EXPLAIN_CHR15_54,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr17_50_237_377 import (
    POLYFUN_EXPLAIN_CHR17_50,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr20_47_653_230 import (
    POLYFUN_EXPLAIN_CHR20_47,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_chrom_pos_alleles import (
    HarmonizeGWASWithReferenceViaAlleles,
)
from mecfs_bio.build_system.wf.base_wf import make_wf

LOCI = {
    "chr1_174": POLYFUN_EXPLAIN_CHR1_174,
    "chr6_26": POLYFUN_EXPLAIN_CHR6_26,
    "chr6_97": POLYFUN_EXPLAIN_CHR6_97,
    "chr15_54": POLYFUN_EXPLAIN_CHR15_54,
    "chr17_50": POLYFUN_EXPLAIN_CHR17_50,
    "chr20_47": POLYFUN_EXPLAIN_CHR20_47,
}
FREQ_COLS = ["EAF", "A1FREQ_CASES", "A1FREQ_CONTROLS"]


def main() -> None:
    meta_to_path = DEFAULT_RUNNER.meta_to_path
    fai = read_fai(HG19)
    mm = np.memmap(HG19, dtype=np.uint8, mode="r")
    input_cache: dict[AssetId, pl.DataFrame] = {}
    for label, locus in LOCI.items():
        task = locus.groups[0].susie_polyfun.gwas_data_task
        assert isinstance(task, HarmonizeGWASWithReferenceViaAlleles)
        deps = {d.asset_id: d for d in task.deps}

        def fetch(asset_id: AssetId) -> Asset:
            return FileAsset(meta_to_path(deps[asset_id].meta))

        with tempfile.TemporaryDirectory() as scratch:
            out = task.execute(scratch_dir=Path(scratch), fetch=fetch, wf=make_wf())
            assert isinstance(out, FileAsset)
            fresh = pl.read_parquet(out.path)

        cached_path = meta_to_path(task.meta)
        cached = pl.read_parquet(cached_path)
        same_as_cache = fresh.sort("SNPID", "EA").equals(cached.sort("SNPID", "EA"))

        src_task = task.gwas_data_task
        if src_task.asset_id not in input_cache:
            input_cache[src_task.asset_id] = (
                scan_dataframe_asset(FileAsset(meta_to_path(src_task.meta)), src_task.meta)
                .collect()
                .to_polars()
                .with_columns(pl.col("EA").cast(pl.Utf8), pl.col("NEA").cast(pl.Utf8))
            )
        src = input_cache[src_task.asset_id]
        cr = task.chrom_range_filter
        assert cr is not None
        src_locus = src.filter(
            (pl.col("CHR") == cr.chrom) & pl.col("POS").is_between(cr.start, cr.end)
        )

        n_src_dup_snpid = src_locus.height - src_locus["SNPID"].n_unique()
        joined = fresh.join(
            src_locus.select(
                "SNPID",
                pl.col("EA").alias("in_EA"),
                pl.col("NEA").alias("in_NEA"),
                pl.col("BETA").alias("in_BETA"),
                *[pl.col(c).alias(f"in_{c}") for c in FREQ_COLS],
            ),
            on="SNPID",
            how="left",
        )
        kept = (joined["EA"] == joined["in_EA"]) & (joined["NEA"] == joined["in_NEA"])
        flipped = (joined["EA"] == joined["in_NEA"]) & (joined["NEA"] == joined["in_EA"])
        fl = joined.filter(flipped)
        is_indel = (fl["EA"].str.len_chars() > 1) | (fl["NEA"].str.len_chars() > 1)

        chrom = fl["CHR"].cast(pl.Utf8).to_numpy()
        pos = fl["POS"].to_numpy()
        in_nea_ref = ref_match(mm, fai, chrom, pos, fl["in_NEA"]) if fl.height else np.array([], bool)
        out_nea_ref = ref_match(mm, fai, chrom, pos, fl["NEA"]) if fl.height else np.array([], bool)

        print(
            f"\n=== {label} ({task.asset_id}) ===\n"
            f"  input rows in chrom range: {src_locus.height:,} (duplicate SNPIDs: {n_src_dup_snpid})\n"
            f"  output rows: {fresh.height:,} | kept orientation: {int(kept.sum()):,} | "
            f"flipped: {fl.height:,} (SNV {int((~is_indel).sum()):,}, indel {int(is_indel.sum()):,}) | "
            f"neither: {int((~kept & ~flipped).sum()):,}\n"
            f"  fresh output identical to cached asset: {same_as_cache}"
        )
        if fl.height:
            beta_ok = np.allclose(fl["BETA"].to_numpy(), -fl["in_BETA"].to_numpy())
            print(f"  flipped rows: BETA negated: {beta_ok}")
            for c in FREQ_COLS:
                if c in fl.columns:
                    became_complement = np.allclose(fl[c].to_numpy(), 1 - fl[f"in_{c}"].to_numpy(), atol=1e-5)
                    unchanged = np.allclose(fl[c].to_numpy(), fl[f"in_{c}"].to_numpy(), atol=1e-5)
                    print(f"  flipped rows: {c} -> 1-{c}: {became_complement}, unchanged: {unchanged}")
            print(
                f"  flipped rows: input NEA == hg19 ref: {int(in_nea_ref.sum())}/{fl.height}; "
                f"output NEA == hg19 ref: {int(out_nea_ref.sum())}/{fl.height}"
            )
            with pl.Config(tbl_cols=-1, tbl_width_chars=250):
                print(fl.select("SNPID", "CHR", "POS", "in_EA", "in_NEA", "EA", "NEA", "in_BETA", "BETA", "in_EAF", "EAF").head(8))


if __name__ == "__main__":
    main()

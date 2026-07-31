"""
Filter a 1000 Genomes EUR reference VCF down to common, biallelic SNVs.

This produces the reusable intermediate asset behind the common-1kg PPP index
mode: biallelic single-nucleotide variants with EUR MAF >= a threshold (default
0.01, ~8.8M variants for the 30x EUR hg38 panel). Columns: CHR, POS[hg38], REF, ALT,
AF. rsIDs are attached downstream (the VCF ID field is chr:pos:ref:alt, not rsIDs).

chrX is included. CHR follows gwaslab's numeric coding, so X arrives here spelled
chrX and leaves as 23.

The VCF arrives as a task dependency rather than a gwaslab download. gwaslab serves
its panels from a URL it can republish in place, and it did: the EUR hg38 panel
gained chrX at some point after October 2025, which silently changed the variant
set. That is unacceptable here, because the PPP variant index built from this asset
fixes the row order that every per-protein beta/se file is written against, so a
source that changes underneath the index invalidates thousands of downstream files.
Depending on a pinned, checksummed copy keeps the index reproducible.

bcftools is invoked through the pixi environment (a pinned dependency), not a system
binary.
"""

from pathlib import Path, PurePath

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.ppp_database.byte_stream_split_parquet import (
    write_byte_stream_split_parquet,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_CODE_FOR_NAME,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
)
from mecfs_bio.util.subproc.run_command import execute_command

# Column names of the 1000 Genomes reference VCF (biallelic reference/alt alleles).
ONEKG_REF_COL = "REF"
ONEKG_ALT_COL = "ALT"
ONEKG_AF_COL = "AF"

FILTERED_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    ONEKG_REF_COL,
    ONEKG_ALT_COL,
    ONEKG_AF_COL,
]

DEFAULT_MAF_THRESHOLD = 0.01

# The VCF spells the non-autosomes with letters (chrX); CHR is an integer column, so
# they are mapped to gwaslab's numeric codes before the cast. Anything not covered
# here reaches the cast as a letter and fails it, which is the intent: an unexpected
# contig should stop the build rather than be dropped or coerced.
_CHROM_NAME_TO_CODE_STR = {
    name: str(code) for name, code in GWASLAB_CHROM_CODE_FOR_NAME.items()
}


def filter_vcf_to_frame(
    vcf_path: Path, maf_threshold: float, work_dir: Path
) -> pl.DataFrame:
    """Run bcftools to select biallelic SNVs with AF in
    [maf_threshold, 1 - maf_threshold], returning CHR, POS, REF, ALT, AF.

    No chromosome is excluded; CHR comes back numerically coded (X as 23).

    bcftools is resolved from the pixi environment (PATH inside the build), not a
    system install.
    """
    tsv_path = work_dir / "filtered_1kg_common.tsv"
    maf_expr = f"INFO/AF>={maf_threshold:.6f} && INFO/AF<={1.0 - maf_threshold:.6f}"
    cmd = [
        "bcftools",
        "view",
        "-v",
        "snps",
        "-m2",
        "-M2",
        "-i",
        f"'{maf_expr}'",
        "-Ou",
        str(vcf_path),
        "|",
        "bcftools",
        "query",
        "-f",
        r"'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n'",
        "-o",
        str(tsv_path),
    ]
    execute_command(cmd)
    return (
        pl.read_csv(
            tsv_path,
            separator="\t",
            has_header=False,
            new_columns=[
                "chrom_str",
                GWASLAB_POS_COL,
                ONEKG_REF_COL,
                ONEKG_ALT_COL,
                ONEKG_AF_COL,
            ],
            schema_overrides={GWASLAB_POS_COL: pl.Int32, ONEKG_AF_COL: pl.Float32},
        )
        # Strip the 'chr' contig prefix, then map X/Y/MT onto their numeric codes.
        .with_columns(
            pl.col("chrom_str")
            .str.replace("^chr", "")
            .replace(_CHROM_NAME_TO_CODE_STR)
            .cast(pl.Int32)
            .alias(GWASLAB_CHROM_COL)
        )
        .select(FILTERED_COLUMNS)
    )


@frozen
class FilterCommon1kgVariantsTask(GeneratingTask):
    """
    Filter a 1000 Genomes EUR reference VCF to common biallelic SNVs.

    thousand_genomes: a task producing the EUR reference VCF as a FileAsset.
    maf_threshold: minor-allele-frequency cutoff (variant kept if AF is within
        [maf_threshold, 1 - maf_threshold]).
    """

    meta: ReferenceFileMeta
    thousand_genomes: Task
    maf_threshold: float

    @property
    def deps(self) -> list[Task]:
        return [self.thousand_genomes]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        vcf_asset = fetch(self.thousand_genomes.asset_id)
        assert isinstance(vcf_asset, FileAsset), (
            f"expected {self.thousand_genomes.asset_id} to be a FileAsset, "
            f"got {type(vcf_asset).__name__}"
        )
        filtered = filter_vcf_to_frame(vcf_asset.path, self.maf_threshold, scratch_dir)
        out_path = scratch_dir / "common_1kg_variants.parquet"
        write_byte_stream_split_parquet(
            filtered, out_path, float_columns=[ONEKG_AF_COL]
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        thousand_genomes: Task,
        asset_id: str,
        maf_threshold: float = DEFAULT_MAF_THRESHOLD,
    ) -> "FilterCommon1kgVariantsTask":
        # The panel this was filtered from is what distinguishes one of these assets
        # from another, so it comes from the source's meta rather than an argument.
        source_meta = thousand_genomes.meta
        assert isinstance(source_meta, ReferenceFileMeta), (
            f"expected {asset_id} source to carry a ReferenceFileMeta, "
            f"got {type(source_meta).__name__}"
        )
        return cls(
            meta=ReferenceFileMeta(
                group="onekg_common_variants",
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename="common_1kg_variants",
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            thousand_genomes=thousand_genomes,
            maf_threshold=maf_threshold,
        )

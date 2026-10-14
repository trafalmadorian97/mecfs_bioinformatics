import narwhals as nw
from pprint import pprint

import polars as pl

from mecfs_bio.build_system.reference.schemas.hg19_snp151_schema_valid_choms import HG19_SNP151_VALID_CHROMS


def get_unique_chroms():
    scan = pl.scan_parquet("assets/base_asset_store/reference_data/genome_annotations/build_37/raw/genome_annotation_database_build_37_as_parquet.parquet.zstd")
    chroms = scan.select(pl.col("chrom")).group_by(pl.col("chrom")).len()
    print(chroms)
    import pdb; pdb.set_trace()
    print("yo")

def get_unique_chroms_duck():
    scan = nw.scan_parquet("assets/base_asset_store/reference_data/genome_annotations/build_37/raw/genome_annotation_database_build_37_as_parquet.parquet.zstd", backend="ibis")
    chroms = scan.select(nw.col("chrom")).group_by(nw.col("chrom")).agg(nw.col("chrom").count().alias("count"))
    print(chroms)
    result = chroms.collect().to_pandas().sort_values("count", ascending=False)
    import pdb; pdb.set_trace()
    print("yo")

def measure_snp_post_dif():
    scan = nw.scan_parquet("assets/base_asset_store/reference_data/genome_annotations/build_37/raw/genome_annotation_database_build_37_as_parquet.parquet.zstd", backend="ibis")
    scan = scan.filter(nw.col("class") == "single").filter(nw.col("chrom").is_in(HG19_SNP151_VALID_CHROMS))
    scan = scan.with_columns((nw.col("chromEnd_zero_based")-nw.col("chromStart_zero_based") ).alias("diff"))
    scan.select(nw.col("diff")).group_by(nw.col("diff")).agg(nw.col("diff").count().alias("count"))
    print(scan)
    result = scan.collect().to_pandas().sort_values("count", ascending=False)
    print(result)
    import pdb; pdb.set_trace()
    print("yo")

def study_result():
    result = nw.scan_parquet(
        "assets/base_asset_store/gwas/ME_CFS/DecodeME/processed/decode_me_gwas_1_liftover_to_build_37_with_rsid.parquet",
        backend="ibis",
    )
    counts = result.select(nw.col("rsID")).group_by(nw.col("rsID")).agg(nw.col("rsID").count().alias("count"))
    counts_filtered = counts.filter(nw.col("count")>1)
    joined = result.join(counts_filtered, on="rsID").collect()
    #observed a few dups here.  I would assume this is a consequence of liftover.
    head = result.head().collect()
    print(head)
    import pdb; pdb.set_trace()
    print("yo")


def what_does_observed_column_look_like_for_snps():
    scan = nw.scan_parquet("assets/base_asset_store/reference_data/genome_annotations/build_37/raw/genome_annotation_database_build_37_as_parquet.parquet.zstd", backend="ibis")
    scan = scan.filter(nw.col("class") == "single")
    observed = scan.select(nw.col("observed")).group_by(nw.col("observed")).agg(nw.col("observed").count().alias("count"))
    print(observed)
    result = observed.collect().to_pandas().sort_values("count", ascending=False)
    import pdb; pdb.set_trace()
    print("yo")


def look_at_dups_in_original():
    pass
    # scan =nw.scan_parquet("assets/base_asset_store/gwas/ME_CFS/DecodeME/processed/decode_me_gwas_1_liftover_to_37_parquet_file.parquet")
    # scan.filter(nw.col(""))



if __name__ == "__main__":
    # get_unique_chroms_duck()
    # measure_snp_post_dif()
    study_result()
    # what_does_observed_column_look_like_for_snps()
    # get_unique_chroms()


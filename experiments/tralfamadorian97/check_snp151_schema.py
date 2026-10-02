import polars as pl
import polars_streaming_csv_decompression

from mecfs_bio.build_system.reference.schemas.hg19_sn151_schema import HG19_SNP151_SCHEMA


def go():
    print("running.. ")
    frame =polars_streaming_csv_decompression.streaming_csv(
        "assets/base_asset_store/reference_data/genome_annotations/build_37/raw/snp151.txt.gz",
        separator="\t"
    )

    print("collecting scheam.. ")
    schema = frame.collect_schema()
    print(schema)
    import pdb; pdb.set_trace()
    print("yo")


def gogo():
    df=pl.read_csv("dummy_snp151.txt", separator="\t", has_header=False,
                new_columns=HG19_SNP151_SCHEMA
                )
    df.write_parquet("dummy_snp_parquet.parquet")


if __name__ == '__main__':
    gogo()
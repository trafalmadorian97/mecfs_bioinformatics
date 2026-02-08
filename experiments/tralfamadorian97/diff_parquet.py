import narwhals as nw
import duckdb
import polars as pl
from pathlib import Path


def diff_parquet(pth1: Path, pth2: Path, out_dir: Path, columns:list[str]):
    # 1. Scan the files (this does not load them into memory)
    out_dir.mkdir(parents=True, exist_ok=True)
    df1 = pl.scan_parquet(pth1)
    df2 = pl.scan_parquet(pth2)

    # 2. Find rows in df1 but not in df2
    added = df1.join(df2, on=columns, how="anti")

    # 3. Find rows in df2 but not in df1
    removed = df2.join(df1, on=columns, how="anti")

    # 4. Combine and stream to disk
    # Use .sink_parquet() for out-of-core streaming
    added_path = out_dir/"added.parquet"
    removed_path = out_dir/"removed.parquet"
    added.sink_parquet(added_path)
    removed.sink_parquet(removed_path)
    print(f"Added rows: {len(pl.scan_parquet(added_path))}")
    print(f"Removed rows: {len(pl.scan_parquet(removed_path))}")
    print("done")

def duckdb_diff_parquet(
pth1: Path, pth2: Path, out_dir: Path
):
    print("running duckdb diff`")
    out_dir.mkdir(parents=True, exist_ok=True)
    # Initialize DuckDB (using a persistent file for "spilling" if needed)
    con = duckdb.connect(database='diff_temp.db')

    # Define the diff query
    # This finds rows in file1 that are not in file2 and vice versa
    query = f"""
        (SELECT * FROM read_parquet('{pth1}') EXCEPT SELECT * FROM read_parquet('{pth2}'))
        UNION ALL
        (SELECT * FROM read_parquet('{pth2}') EXCEPT SELECT * FROM read_parquet('{pth1}'))
    """

    # Stream the result directly to a new parquet file
    query = f"COPY ({query}) TO '{out_dir/"diff.parquet"}' (FORMAT PARQUET);"
    print(f"Running:\n {query}")
    con.execute(query)

    print("done")



def diff_narwhals(
        pth1: Path, pth2: Path, out_dir: Path
):

    out_dir.mkdir(parents=True, exist_ok=True)
    print("running narwhals diff`")

    duckdb.sql(f"SET memory_limit='{2}GB'")
    df1= nw.scan_parquet(pth1,backend="ibis")
    df2 = nw.scan_parquet(pth2,backend="ibis")
    df1_minus_df2 = df1.join(df2, on=df1.columns, how="anti")
    df2_minus_df1 = df2.join(df1, on=df2.columns, how="anti")

    df1_minus_df2.sink_parquet(out_dir/"df1_minus_df2.parquet")
    # df2_minus_df1.sink_parquet(out_dir/"df2_minus_df1.parquet")
    print("done")

def compare_proc_parquet():
    diff_parquet(
        Path("/tmp/pytest-of-paiforsyth/pytest-118/test_run_hba_magma0/asset_store/reference_data/db_snp_reference_data/build_37/annovar/db_snp150_annovar_proc_parquet.parquet.zstd"),
Path("assets/base_asset_store/reference_data/db_snp_reference_data/build_37/annovar/db_snp150_annovar_proc_parquet.parquet.zstd"  ),
out_dir=Path("data")/"proc_parquet_comparison",
       columns=[ "rsid"]
    )

if __name__ == "__main__":
    compare_proc_parquet()
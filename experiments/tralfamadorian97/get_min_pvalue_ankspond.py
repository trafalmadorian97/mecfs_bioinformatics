import polars as pl
def go():
    import duckdb

    # Define the path to your gzipped TSV file
    file_path = 'assets/base_asset_store/gwas/ankylosing_spondylitis/uk_biobank/raw/GCST90474065.h.tsv.gz'

    # Connect to an in-memory DuckDB database
    con = duckdb.connect()

    # Query the minimum p-value using read_csv_auto
    # DuckDB natively handles gzip decompression and detects headers/delimiters
    query = f"""
         SELECT MIN(p_value) AS min_p
         FROM read_csv_auto('{file_path}', delim='\t', header=True)
     """

    # Execute the query and fetch the single row result
    result = con.execute(query).fetchone()
    min_p_value = result[0]

    print(f"The minimum p-value is: {min_p_value}")

if __name__=="__main__":
    go()
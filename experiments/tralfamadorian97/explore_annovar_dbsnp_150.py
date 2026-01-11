import duckdb





"""
Goal: determine why there are apparently duplicate variants in dbsnp150
"""



def go():
    df= duckdb.sql(
        f"""
        select int_chrom, POS, ALT, REF,  COUNT(*)
        from 'assets/base_asset_store/reference_data/db_snp_reference_data/build_37/processed/db_snp150_annovar_proc_parquet_rename.parquet'
        group by int_chrom, POS, ALT, REF
        order by COUNT(*) desc
        limit 10
        """
    ).df()
    print(df)
    import pdb; pdb.set_trace()
    print("yo")
    # Here is the result:
    """
    0         21   21356260     CATA   -            17
1         21   21356260    GCATA   G            17
2          9  107080214    GATAC   G            17
3          9  107080214     ATAC   -            17
4          7  102371702     AAGG   -            16
5         10   32588468  CCATCAT   C            16
6         10   32588468   CATCAT   -            16
7          5  106976655     GGTA   -            16
8          7  102371702    AAAGG   A            16
9          5  106976655    TGGTA   T            16
    """

def look_at_duplicated():
    df = duckdb.sql(
        """
        select *
        from 'assets/base_asset_store/reference_data/db_snp_reference_data/build_37/processed/db_snp150_annovar_proc_parquet_rename.parquet'
        where int_chrom=21 and POS=21356260 and ALT='GCATA' and ref='G'
        """
    ).df()
    print(df)
    import pdb; pdb.set_trace()
    print("yo")


def count_dupliated():
    df= duckdb.sql(
        f"""
        select num_dups, COUNT(*)
        from
        (
        select int_chrom, POS, ALT, REF,  COUNT(*) as num_dups
        from 'assets/base_asset_store/reference_data/db_snp_reference_data/build_37/processed/db_snp150_annovar_proc_parquet_rename.parquet'
        group by int_chrom, POS, ALT, REF
        )
        group by num_dups
        order by COUNT(*)
        """
    ).df()
    print(df)
    import pdb; pdb.set_trace()
    print("yo")
    """
    Result:
        num_dups  count_star()
0         15             2
1         17             4
2         14             6
3         16             6
4         12             8
5         11            12
6         10            22
7          9            26
8          8           102
9          7           186
10         6           768
11         5          6650
12         4         62380
13         3        440502
14         2       2477667
15         1     350160447
    """


def count_dupliated_ibd():
    df= duckdb.sql(
        f"""
        select num_dups, COUNT(*)
        from
        (
        select int_chrom, POS, ALT, REF,  COUNT(*) as num_dups
        from 'assets/base_asset_store/gwas/inflammatory_bowl_disease/liu_et_al_2023_ibd_meta/processed/liu_et_al_2023_ibd_eur_harmonize_assign_rsids_via_snp150_annovar_with_dups.parquet'
        group by int_chrom, POS, ALT, REF
        )
        group by num_dups
        order by COUNT(*)
        """
    ).df()
    print(df)
    import pdb; pdb.set_trace()
    print("yo")
    """
       num_dups  count_star()
0        10             1
1         9             2
2         8             7
3         7            15
4         6            39
5         5           129
6         4          1171
7         3         12458
8         2        133877
9         1       9403846
So on the order of 1% of rsids here actually refer to duplicates of the same genetic variant
    """

if __name__ == "__main__":
    count_dupliated_ibd()
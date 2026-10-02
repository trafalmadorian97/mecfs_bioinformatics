import pandas as pd


def compar_ensemble_ids():
   df_ga = pd.read_csv("assets/base_asset_store/gwas/ME_CFS/DecodeME/processed/magma/decode_me_gwas_1_build_37_magma_ensemble_gene_analysis/gene_analysis_output.genes.out",
                    delim_whitespace=True,)
   df_ref = pd.read_csv(
       "assets/base_asset_store/reference_data/rna_seq_data/gtex/processed/gtex_v10_rna_seq_median_tissue_expression_prep_for_magma.tsv",
       sep="\t",
   )
   genes_ga = set(df_ga["GENE"].tolist())
   df_ref[["GENE","ver"]] = df_ref["Name"].str.split(".", expand=True)
   genes_ref = set(df_ref["GENE"].tolist())
   common = genes_ref & genes_ga
   print(f"number of genes in common is:{ len(common)}")
   import pdb; pdb.set_trace()
   print("yo")


if __name__ == "__main__":
    compar_ensemble_ids()

"""
Download and preprocess data from the UK Biobank Pharma Proteomics Project GWAS of plasma BTN1A1
"""

from mecfs_bio.asset_generator.ukbb_ppp_single_gwas_generator import ubbb_ppp_gwas_prep

BTN1A1_UKBB_PPP_GWAS_PROCESS = ubbb_ppp_gwas_prep(
    gene_name="btn1a1",
    syn_id="syn52362950",
    expected_filename="BTN1A1_Q13410_OID31363_v1_Oncology_II.tar",
    region_to_plot_chrom=6,
    region_to_plot_37_pos=26465384,
)

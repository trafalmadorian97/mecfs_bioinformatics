from mecfs_bio.build_system.task.gwaslab.gwaslab_cell_analysis_by_sldsc import (
    ldcts_processed_to_raw,
    ldcts_raw_to_processed,
)

_sample_ldcts_file_contents = """Adipose_Subcutaneous\tMulti_tissue_gene_expr_1000Gv3_ldscores/GTEx.1.,Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.control.
Adipose_Visceral_(Omentum)\tMulti_tissue_gene_expr_1000Gv3_ldscores/GTEx.2.,Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.control.
Adrenal_Gland\tMulti_tissue_gene_expr_1000Gv3_ldscores/GTEx.3.,Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.control.
Artery_Aorta\tMulti_tissue_gene_expr_1000Gv3_ldscores/GTEx.4.,Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.control."""


def test_round_trip_ldcts_file_contents():
    processed = ldcts_raw_to_processed(_sample_ldcts_file_contents)
    raw = ldcts_processed_to_raw(processed)
    assert raw == _sample_ldcts_file_contents

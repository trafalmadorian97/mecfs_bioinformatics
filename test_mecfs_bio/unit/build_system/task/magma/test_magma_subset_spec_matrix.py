import pandas as pd

from mecfs_bio.build_system.task.magma.magma_subset_specificity_matrix_using_top_labels import (
    get_spec_matrix_filtered,
)


def test_magma_subset_spec_matrix():
    dummy_spec_matrix = pd.DataFrame(
        {
            "GENE": [1, 2, 3],
            "CLUSTER1": [1, 1, 1],
            "CLUSTER2": [2, 1, 2],
            "CLUSTER3": [2, 1, 3],
        }
    )
    dummy_covar_result = pd.DataFrame(
        {
            "VARIABLE": ["CLUSTER1", "CLUSTER2", "CLUSTER3"],
            "TYPE": [
                "COVAR",
                "COVAR",
                "COVAR",
            ],
            "NGENES": [3, 3, 3],
            "BETA": [0.01, 0, -0.02],
            "P": [0.01, 0.1, 0.000001],
        }
    )
    spec_matrix_filtered = get_spec_matrix_filtered(
        covar_result=dummy_covar_result,
        spec_matrix=dummy_spec_matrix,
        nominal_sig_level=0.01,
    )
    assert spec_matrix_filtered.columns.tolist() == ["GENE", "CLUSTER3"]

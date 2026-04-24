from pathlib import Path

from mecfs_bio.constants.curated_gene_set_collections.curated_mecfs_gene_sets import \
    CURATED_POTENTIAL_MECFS_GENE_SETS_ALL
from mecfs_bio.util.gene_set.msigdb_lookup import gene_set_jaccard_index


def go():
    result = gene_set_jaccard_index(
        parquet_path=Path("assets/base_asset_store/reference_data/gene_set_data/msigdb/extracted/msigdb_human_gene_sets_table_parquet_from_sqllite.parquet.zstd"),
        specs=CURATED_POTENTIAL_MECFS_GENE_SETS_ALL
    )
    print(result)
    import pdb; pdb.set_trace()
    print(result)

"""
Top overlap:

                                      standard_name_1                                    standard_name_2  jaccard_index  n_intersection  n_union
0           GOBP_AUTONOMIC_NERVOUS_SYSTEM_DEVELOPMENT        GOBP_SYMPATHETIC_NERVOUS_SYSTEM_DEVELOPMENT       0.466667              21       45
1           GOBP_AUTONOMIC_NERVOUS_SYSTEM_DEVELOPMENT    GOBP_PARASYMPATHETIC_NERVOUS_SYSTEM_DEVELOPMENT       0.400000              18       45
2                  HALLMARK_INTERFERON_GAMMA_RESPONSE                 HALLMARK_INTERFERON_ALPHA_RESPONSE       0.325893              73      224
3           GOBP_AUTONOMIC_NERVOUS_SYSTEM_DEVELOPMENT            GOBP_ENTERIC_NERVOUS_SYSTEM_DEVELOPMENT       0.311111              14       45
4         GOBP_SYMPATHETIC_NERVOUS_SYSTEM_DEVELOPMENT    GOBP_PARASYMPATHETIC_NERVOUS_SYSTEM_DEVELOPMENT       0.258065               8       31
5                              REACTOME_TCR_SIGNALING      REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR       0.251185              53      211
6                              REACTOME_TCR_SIGNALING  REACTOME_REGULATION_OF_T_CELL_ACTIVATION_BY_CD...       0.233227              73      313
7                                    HALLMARK_HYPOXIA                                HALLMARK_GLYCOLYSIS       0.194030              65      335
8             HP_DISTAL_PERIPHERAL_SENSORY_NEUROPATHY         HP_MITOCHONDRIAL_RESPIRATORY_CHAIN_DEFECTS       0.181818               6       33
9       REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR  REACTOME_REGULATION_OF_T_CELL_ACTIVATION_BY_CD...       0.126316              48      380
10           GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION         REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION       0.122172              27      221
11                GOBP_REGULATION_OF_VASOCONSTRICTION                  GOBP_REGULATION_OF_BLOOD_PRESSURE       0.117647              28      238
12            KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION                   KEGG_ARACHIDONIC_ACID_METABOLISM       0.116129              18      155
13                KEGG_NEUROTROPHIN_SIGNALING_PATHWAY                               BIOCARTA_NGF_PATHWAY       0.114504              15      131
14                 HALLMARK_OXIDATIVE_PHOSPHORYLATION                       KEGG_CITRATE_CYCLE_TCA_CYCLE       0.104762              22      210
15                             GOBP_RESPONSE_TO_VIRUS                 HALLMARK_INTERFERON_GAMMA_RESPONSE       0.099656              58      582
16  GOBP_MICROGLIAL_CELL_ACTIVATION_INVOLVED_IN_IM...                          GOBP_ASTROCYTE_ACTIVATION       0.096774               3       31
17            GOBP_ENTERIC_NERVOUS_SYSTEM_DEVELOPMENT        GOBP_SYMPATHETIC_NERVOUS_SYSTEM_DEVELOPMENT       0.093750               3       32
18                               BIOCARTA_NGF_PATHWAY                                  PID_TRAIL_PATHWAY       0.090909               4       44
19     GOBP_NITRIC_OXIDE_MEDIATED_SIGNAL_TRANSDUCTION                         WP_EFFECTS_OF_NITRIC_OXIDE       0.088235               3       34

"""

if __name__ == '__main__':
    go()
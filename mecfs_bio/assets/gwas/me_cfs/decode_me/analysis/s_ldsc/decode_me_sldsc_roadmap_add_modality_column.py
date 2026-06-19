"""
Plan:
Split roadmap into different tables.  One per modality.  To do this, need a list of modalities

Then, add both Bonferroni and FDR to both tables


List of assays
>>> df = pd.read_csv("assets/base_asset_store/gwas/ME_CFS/DecodeME/analysis/decode_me_gwas_1_multi_tissue_chromatin_add_labels.csv")
>>> df
                                              Name   Coefficient  Coefficient_std_error  Coefficient_P_value  _Corrected P Value_  Reject Null                                  Cell Epigenetic_Assay                             Cell Type   Category
0                        Fetal_Brain_Female__DNase  1.022580e-07           1.681961e-08         6.020849e-10         2.944195e-07         True                    fetal_brain_female            DNase                    Fetal_Brain_Female        CNS
1    Brain_Dorsolateral_Prefrontal_Cortex__H3K27ac  5.080309e-08           8.568361e-09         1.522540e-09         3.508038e-07         True  brain_dorsolateral_prefrontal_cortex          H3K27ac  Brain_Dorsolateral_Prefrontal_Cortex        CNS
2                          Fetal_Brain_Male__DNase  9.543345e-08           1.625215e-08         2.152170e-09         3.508038e-07         True                      fetal_brain_male            DNase                      Fetal_Brain_Male        CNS
3    Brain_Dorsolateral_Prefrontal_Cortex__H3K4me3  1.706178e-07           2.991470e-08         5.869385e-09         7.175323e-07         True  brain_dorsolateral_prefrontal_cortex          H3K4me3  Brain_Dorsolateral_Prefrontal_Cortex        CNS
4                        Fetal_Brain_Male__H3K4me1  3.763501e-08           6.654997e-09         7.785455e-09         7.614175e-07         True                      fetal_brain_male          H3K4me1                      Fetal_Brain_Male        CNS
..                                             ...           ...                    ...                  ...                  ...          ...                                   ...              ...                                   ...        ...
484                        Stomach_Mucosa__H3K4me1 -2.157466e-08           5.730949e-09         9.999166e-01         9.999625e-01        False                        stomach_mucosa          H3K4me1                        Stomach_Mucosa  Digestive
485               Duodenum_Smooth_Muscle__H3K36me3 -3.772878e-08           9.995744e-09         9.999198e-01         9.999625e-01        False                duodenum_smooth_muscle         H3K36me3                Duodenum_Smooth_Muscle  Digestive
486                         Stomach_Mucosa__H3K9ac -5.798558e-08           1.515497e-08         9.999349e-01         9.999625e-01        False                        stomach_mucosa           H3K9ac                        Stomach_Mucosa  Digestive
487                 Rectal_Mucosa_Donor_29__H3K9ac -8.208379e-08           2.113282e-08         9.999487e-01         9.999625e-01        False                rectal_mucosa_donor_29           H3K9ac                Rectal_Mucosa_Donor_29  Digestive
488                  Rectal_Smooth_Muscle__H3K4me3 -1.034162e-07           2.611609e-08         9.999625e-01         9.999625e-01        False                  rectal_smooth_muscle          H3K4me3                  Rectal_Smooth_Muscle  Digestive

[489 rows x 10 columns]
>>> df["Epigenetic_Assay"]
0         DNase
1       H3K27ac
2         DNase
3       H3K4me3
4       H3K4me1
         ...
484     H3K4me1
485    H3K36me3
486      H3K9ac
487      H3K9ac
488     H3K4me3
Name: Epigenetic_Assay, Length: 489, dtype: object
>>> df["Epigenetic_Assay"].unique()
array(['DNase', 'H3K27ac', 'H3K4me3', 'H3K4me1', 'H3K9ac', 'H3K36me3'],
      dtype=object)
"""
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.s_ldsc.decode_me_sldsc import DECODE_ME_S_LDSC
from mecfs_bio.build_system.task.pipe_dataframe_task import PipeDataFrameTask, CSVOutFormat
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue

HYPOTHESIS_TESTING_COLUMNS = [
    "_Corrected P Value_",
    "Reject Null"
]

DROP_HP_COLS_PIPE = DropColPipe(
    cols_to_drop=HYPOTHESIS_TESTING_COLUMNS,
)

DECODE_ME_S_LDSC_ASSAY_SPECIFIC = [PipeDataFrameTask.create(
    source_task=DECODE_ME_S_LDSC.partitioned_tasks["multi_tissue_chromatin"].add_categories_task_unwrap,
    asset_id=f"decode_me_s_ldsc_result_multi_tissue_chromatin_{assay}_only",
    out_format=CSVOutFormat(sep=","),
    pipes=[
        FilterRowsByValue(
            target_column="Epigenetic_Assay",
            valid_values=[assay]
        ),
        DROP_HP_COLS_PIPE
    ]
) for assay in ["DNase","H3K27ac","H3K4me3","H3K4me1","H3K9ac","H3K36me3"]]



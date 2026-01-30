from mecfs_bio.assets.reference_data.rna_seq_data.gtex_v10_median_tissue_expression_rna_seq_raw import (
    GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    CSVOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.add_average_pipe import AddAveragePipe
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_min import FilterRowsByMin
from mecfs_bio.build_system.task.pipes.move_col_to_front_pipe import MoveColToFrontPipe
from mecfs_bio.build_system.task.pipes.shifted_log_pipe import ShiftedLogPipe
from mecfs_bio.build_system.task.pipes.str_split_exact_col import SplitExactColPipe
from mecfs_bio.build_system.task.pipes.winsorize_all import WinsorizeAllPipe

GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA_BRAIN_AVERAGE = PipeDataFrameTask.create(
    source_task=GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ,  # GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_EXTRACTED,
    asset_id="gtex_v10_rna_seq_median_tissue_expression_prep_for_magma_brain_average",
    out_format=CSVOutFormat(sep="\t"),
    pipes=[
        DropColPipe(cols_to_drop=["Description"]),
        FilterRowsByMin(min_value=1, exclude_columns=["Name"]),
        WinsorizeAllPipe(max_value=50, cols_to_exclude=["Name"]),
        SplitExactColPipe(
            col_to_split="Name", split_by=".", new_col_names=("Gene", "Version")
        ),
        DropColPipe(cols_to_drop=["Version", "Name"]),
        ShiftedLogPipe(cols_to_exclude=["Gene"], base=2, pseudocount=1),
        AddAveragePipe(
            cols_to_exclude=[
                "Gene",
                "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", 
                "Artery_Tibial", "Bladder", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", 
                "Cervix_Ectocervix", "Cervix_Endocervix", "Colon_Sigmoid", "Colon_Transverse", "Colon_Transverse_Mixed_Cell", 
                "Colon_Transverse_Mucosa", "Colon_Transverse_Muscularis", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", 
                "Esophagus_Muscularis", "Fallopian_Tube", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", 
                "Kidney_Medulla", "Liver", "Liver_Hepatocyte", "Liver_Mixed_Cell", "Liver_Portal_Tract", "Lung", "Minor_Salivary_Gland", 
                "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pancreas_Acini", "Pancreas_Islets", "Pancreas_Mixed_Cell", 
                "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", 
                "Small_Intestine_Terminal_Ileum_Lymphode_Aggregate", "Small_Intestine_Terminal_Ileum_Mixed_Cell", "Spleen", "Stomach", 
                "Stomach_Mixed_Cell", "Stomach_Mucosa", "Stomach_Muscularis", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood",
                ], 
            ),
        MoveColToFrontPipe(target_col="Gene"),
    ],
    backend="polars",
)

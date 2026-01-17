from mecfs_bio.build_system.task.get_uniprot_reference_data_task import (
    GetUniProtReferenceDataTask,
)

UNIPROT_LOOKUP = GetUniProtReferenceDataTask.create(
    "uniprot_lookup_table",
)

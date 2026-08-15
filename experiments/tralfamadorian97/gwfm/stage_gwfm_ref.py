from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.sbayes_rc_gwfm_ref.imputed_13m_v1.stage_reference_data_imputed_13m_v1 import \
    SBAYESRC_GWFM_IMPUTED13M_V1
from mecfs_bio.build_system.task.sbayesrc.stage_gwfm_reference_task import (
    StageGwfmReferenceTask,
)



def main() -> None:
    task =SBAYESRC_GWFM_IMPUTED13M_V1
    DEFAULT_RUNNER.run([task])


if __name__ == "__main__":
    main()

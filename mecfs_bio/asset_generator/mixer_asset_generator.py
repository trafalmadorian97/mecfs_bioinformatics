from typing import Mapping

from attrs import frozen

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.mixer_bivariate_task import MixerDataSource, MixerTask, UnivariateMode


@frozen
class UnivariateMixerTasks:
    tasks: Mapping[int, Task]
    def terminal_tasks(self)->list[Task]:
        return list(self.tasks.values())


def univariate_mixer_asset_generator(
        base_name:str,
        trait_1_source: MixerDataSource,
        reference_data_directory_task: Task,
):
    tasks= {}
    for rep in range(1, 21):
        tasks[rep] = MixerTask.create(
            asset_id=base_name+ f"_standard_mixer_rep_{rep}",
            trait_1_source=trait_1_source,
            mixer_mode=UnivariateMode(),
            ref_data_directory_task=reference_data_directory_task,
            reps_to_perform=[rep],
        )
    return UnivariateMixerTasks(tasks)
from typing import Mapping, Sequence

from attrs import frozen

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.mixer.mixer_task import MixerDataSource, MixerTask, UnivariateMode
from mecfs_bio.build_system.task.mixer.mixer_univariate_combine import MixerUnivariateCombine, MixerRunSource
from mecfs_bio.build_system.task.mixer.mixer_univariate_plot import MixerUnivariatePlot


@frozen
class UnivariateMixerTasks:
    run_tasks: Mapping[int, Task]
    combine_task: Task
    plot_task: Task
    def terminal_tasks(self)->list[Task]:
        # return list(self.run_tasks.values())
        return [self.plot_task]


def univariate_mixer_asset_generator(
        base_name:str,
        name_in_plot:str,
        trait_1_source: MixerDataSource,
        reference_data_directory_task: Task,
        reps:Sequence[int]=tuple(range(1,21)),
        threads: int=4
):
    tasks= {}
    for rep in reps:
        tasks[rep] = MixerTask.create(
            asset_id=base_name+ f"_standard_mixer_rep_{rep}",
            trait_1_source=trait_1_source,
            mixer_mode=UnivariateMode(),
            ref_data_directory_task=reference_data_directory_task,
            reps_to_perform=[rep],
            threads=threads,
        )
    combine_task = MixerUnivariateCombine.create(
        asset_id=base_name+"_univariate_mixer_combine",
        mixer_source_runs=[
            MixerRunSource(
               task=task,
                rep=num
            ) for num, task in tasks.items()
        ],
        trait_name=name_in_plot,
    )
    plot_task=MixerUnivariatePlot.create(
        asset_id=base_name+"_univariate_mixer_plot",
        combine_task=combine_task,
        trait_name=name_in_plot,
    )
    return UnivariateMixerTasks(tasks,
                                combine_task=combine_task,
                                plot_task=plot_task)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.mixer_bivariate_task import MixerBivariateTask, MixerDataSource


def go():
    MixerBivariateTask(
        meta=SimpleFileMeta("dummy"),
        trait_1_source=MixerDataSource(
            task=FakeTask("fake"),
            alias="task1",
            sample_info=None,
            pipe=None
        ),
        trait_2_source=MixerDataSource(
            task=FakeTask("fake"),
            alias="task1",
            sample_info=None,
            pipe=None
        )
    )
    import pdb; pdb.set_trace()
    print("yo")

if __name__ == '__main__':
    go()
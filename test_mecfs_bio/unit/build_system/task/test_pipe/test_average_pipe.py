import narwhals
import pandas as pd

from mecfs_bio.build_system.task.pipes.add_average_pipe import AddAveragePipe


def test_add_average_pipe():
    df = pd.DataFrame(
        {
            "nm": ["c", "d", "e"],
            "a": [0, 5, 4],
            "b": [50, 5, 6],
        }
    )
    nw_df = narwhals.from_native(df).lazy()
    pipe = AddAveragePipe(cols_to_exclude=["nm"])
    result = pipe.process(nw_df)
    result = result.collect().to_pandas()
    pd.testing.assert_frame_equal(
        result,
        pd.DataFrame(
            {
                "nm": ["c", "d", "e"],
                "a": [0, 5, 4],
                "b": [50, 5, 6],
                "Average": [25.0, 5.0, 5.0],
            }
        ),
    )

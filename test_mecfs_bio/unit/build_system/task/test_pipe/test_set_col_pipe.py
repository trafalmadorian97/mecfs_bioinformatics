import narwhals
import pandas as pd

from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe


def test_set_col_pipe():
    df = pd.DataFrame({"a": [1, 2, 3]})
    nw_frame = narwhals.from_native(df).lazy()
    pipe = SetColToConstantPipe(col_name="b", constant=5)
    result = pipe.process(nw_frame)
    pd_result = result.collect().to_pandas()
    assert (pd_result["b"] == 5).all()

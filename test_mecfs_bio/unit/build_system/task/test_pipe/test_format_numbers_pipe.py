import pandas as pd

from mecfs_bio.build_system.task.pipes.format_numbers_pipe import FormatFloatNumbersPipe


def test_format_numbers_pipe():
    df = pd.DataFrame(
        {
            "col1": ["a", "b", "c"],
            "col2": ["0.0012345678", "1234.5678", "blarg"],
        }
    )
    result = FormatFloatNumbersPipe(col="col2", format_str=".4g").process_pandas(df)
    assert result["col2"].tolist() == ["0.001235", "1235", "blarg"]

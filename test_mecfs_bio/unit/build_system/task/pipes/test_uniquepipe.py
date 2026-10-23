import narwhals
import pandas as pd

from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL, GWASLAB_RSID_COL


def test_unique_pipe_keep_none_drops_every_row_of_a_duplicated_key() -> None:
    frame = pd.DataFrame(
        {
            GWASLAB_RSID_COL: ["rs1", "rs2", "rs2", "rs3"],
            GWASLAB_POS_COL: [1, 2, 3, 4],
        }
    )
    result = (
        UniquePipe(by=[GWASLAB_RSID_COL], keep="none")
        .process(narwhals.from_native(frame).lazy())
        .collect()
        .to_polars()
    )
    assert sorted(result[GWASLAB_POS_COL].to_list()) == [1, 4]

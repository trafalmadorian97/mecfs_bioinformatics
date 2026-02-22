import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_INFO_SCORE_COL


@frozen
class FilterRowsByInfoScorePipe(DataProcessingPipe):
    min_score: float
    info_col: str = GWASLAB_INFO_SCORE_COL

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.filter(narwhals.col(self.info_col) >= self.min_score)

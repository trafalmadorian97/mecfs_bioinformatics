import scipy.stats

import narwhals

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_P_COL, GWASLAB_BETA_COL, GWASLAB_SE_COL


class ComputePFromBetaSEPipeIfNeeded(DataProcessingPipe):
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:

        schema = x.collect_schema()
        if GWASLAB_P_COL in schema:
            return x
        assert GWASLAB_SE_COL in schema
        assert GWASLAB_BETA_COL in schema
        collected = x.collect().to_pandas()
        z = collected[GWASLAB_BETA_COL]/collected[GWASLAB_SE_COL]
        collected[GWASLAB_P_COL]= 2* scipy.stats.norm.sf(abs(z))
        return narwhals.from_native(collected).lazy()
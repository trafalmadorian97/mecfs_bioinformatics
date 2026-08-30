"""Force-rebuild the l1 demonstrator explainability plot to visually check the
matplotlib rework (code changes don't invalidate the trace, so pass the plot as
must_rebuild_transitive). Requires the /mnt/d asset-store remap to be mounted.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr1_174_128_548 import (
    POLYFUN_EXPLAIN_CHR1_174,
)

plot = POLYFUN_EXPLAIN_CHR1_174.groups[0].plot  # l1
print("rebuilding target:", plot.asset_id)
result = DEFAULT_RUNNER.run(targets=[plot], must_rebuild_transitive=[plot])
print("done:", result[plot.asset_id])

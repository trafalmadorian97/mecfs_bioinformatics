"""Force-rebuild the l1 demonstrator explainability plot to visually check the
matplotlib rework (code changes don't invalidate the trace, so pass the plot as
must_rebuild_transitive). Requires the /mnt/d asset-store remap to be mounted.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr1_174_128_548 import (
    POLYFUN_EXPLAIN_CHR1_174,
)

group = POLYFUN_EXPLAIN_CHR1_174.groups[0]  # l1
plot = group.plot
contrast = group.contrast
print("rebuilding target:", plot.asset_id)
# Force both the contrast task (which now writes callouts.parquet) and the plot:
# code changes don't invalidate the trace, so the on-disk contrast asset would
# otherwise be reused without the new output.
result = DEFAULT_RUNNER.run(
    targets=[plot], must_rebuild_transitive=[plot, contrast]
)
print("done:", result[plot.asset_id])

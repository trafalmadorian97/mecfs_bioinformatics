"""
Timing probe for the R-backed SuSiE fine-mapping task, to size the cost of calling SuSiE
repeatedly in a SuSiE-coloc pipeline.

Forces a rebuild of the base (L=10) SuSiE fine-mapping task at the DecodeME chr1:174,128,548
locus (Broad UKBB LD matrix chr1_173000001_176000001, already cached locally), and nothing
downstream of it. The debug log timestamps then give a per-stage breakdown: loading and
aligning the LD matrix, the PSD repair, estimate_s_rss and kriging_rss diagnostics, and
susie_rss itself.
"""

import time

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr1_174_128_548_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES,
)

if __name__ == "__main__":
    task = DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES.susie_finemap_task
    start = time.perf_counter()
    DEFAULT_RUNNER.run(
        [task],
        incremental_save=True,
        must_rebuild_transitive=[task],
    )
    print(f"TOTAL wall-clock for runner.run: {time.perf_counter() - start:.1f} s")

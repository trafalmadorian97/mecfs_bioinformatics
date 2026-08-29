from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr1_174_128_548 import (
    POLYFUN_EXPLAIN_CHR1_174,
)


def test_demonstrator_wires_eight_susie_runs():
    # 4 run configs * (uniform, polyfun, contrast, plot) = 16 terminal tasks.
    assert len(POLYFUN_EXPLAIN_CHR1_174.terminal_tasks()) == 16
    susie = [
        t
        for g in POLYFUN_EXPLAIN_CHR1_174.groups
        for t in (g.susie_uniform, g.susie_polyfun)
    ]
    assert len({t.asset_id for t in susie}) == 8
    # The polyfun member of every pair carries the precomputed prior; the uniform
    # member does not.
    for g in POLYFUN_EXPLAIN_CHR1_174.groups:
        assert g.susie_polyfun.prior_info is not None
        assert g.susie_uniform.prior_info is None

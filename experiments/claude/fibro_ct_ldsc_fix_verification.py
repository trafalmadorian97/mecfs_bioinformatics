"""
Verify the degenerate-Z fix on the real assets, driving the actual task code path
(load_and_preprocess_sumstats + get_compatible_snps_polars + estimate_rg_by_ldsc)
rather than a hand-rolled imitation.

Reproduces the Fibromyalgia x Schizophrenia pair from CT_LDSC_INITIAL_ASSET_GENERATOR,
including the schizophrenia total-N pipe that the production config applies.

Run:
    pixi r python experiments/claude/fibro_ct_ldsc_fix_verification.py \
        2>&1 | tee experiments/claude/logs/fibro_ct_ldsc_fix_verification.log
"""

from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
    FilterSettings,
    SumstatsSource,
    get_compatible_snps_polars,
    get_prev_options,
    load_and_preprocess_sumstats,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
)

ASSET_STORE = Path("assets/base_asset_store")
FIBRO_PICKLE = (
    ASSET_STORE
    / "gwas/fibromyalgia/kerrebijn_et_al/gwaslab_sumstats/kerrebijin_fibro_sumstats_37.pickle"
)
SCZ_PICKLE = (
    ASSET_STORE
    / "gwas/schizophrenia/pgc_2022/gwaslab_sumstats/pgc_2022_sch_sumstats_37.pickle"
)
LD_DIR = (
    ASSET_STORE
    / "reference_data/linkage_disequilibrium_scores/thousand_genomes_phase_3_v1/extracted"
    / "thousand_genomes_phase_3_v1_eur_ld_scores_extracted"
)


def make_source(alias: str, pickle_path: Path, sample_info, pipe) -> tuple:
    asset_id = AssetId(f"{alias.lower()}_sumstats")
    task = FakeTask(
        meta=GWASLabSumStatsMeta(id=asset_id, trait=alias.lower(), project="verify")
    )
    source = SumstatsSource(
        task=task, alias=alias, sample_info=sample_info, pipe=pipe
    )

    def fetch(requested: AssetId) -> Asset:
        assert requested == asset_id
        return FileAsset(pickle_path)

    return source, fetch


def main() -> None:
    fibro_source, fibro_fetch = make_source(
        "Fibromyalgia",
        FIBRO_PICKLE,
        BinaryPhenotypeSampleInfo(
            sample_prevalence=0.5, estimated_population_prevalence=0.027
        ),
        IdentityPipe(),
    )
    scz_source, scz_fetch = make_source(
        "Schizophrenia",
        SCZ_PICKLE,
        BinaryPhenotypeSampleInfo(
            sample_prevalence=0.408, estimated_population_prevalence=0.01
        ),
        SetColToConstantPipe(GWASLAB_SAMPLE_SIZE_COLUMN, constant=130644),
    )

    settings = FilterSettings()
    fibro, fibro_name, fibro_info = load_and_preprocess_sumstats(
        source=fibro_source, fetch=fibro_fetch, settings=settings, build="19"
    )
    scz, scz_name, scz_info = load_and_preprocess_sumstats(
        source=scz_source, fetch=scz_fetch, settings=settings, build="19"
    )

    print(f"\nfibro rows reaching LDSC: {len(fibro.data)}")
    print(f"fibro zero-SE remaining:  {int((fibro.data['SE'] == 0).sum())}")

    compatible = get_compatible_snps_polars(fibro.data, scz.data)
    scz.data = scz.data.loc[
        scz.data[GWASLAB_RSID_COL].isin(compatible[GWASLAB_RSID_COL])
    ]

    options = get_prev_options(trait_1_prev=fibro_info, trait_2_prev=scz_info)
    fibro.estimate_rg_by_ldsc(
        other_traits=[scz],
        rg=f"{fibro_name},{scz_name}",
        ref_ld_chr=str(LD_DIR) + "/LDscore.@",
        w_ld_chr=str(LD_DIR) + "/LDscore.@",
        build="19",
        **options,
    )

    print("\n===== RESULT =====")
    print(fibro.ldsc_rg.to_string())


if __name__ == "__main__":
    main()

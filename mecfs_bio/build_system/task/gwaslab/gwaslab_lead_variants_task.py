from pathlib import Path

import gwaslab as gl
from attrs import field, frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_lead_variants_meta import (
    GWASLabLeadVariantsMeta,
)
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.read_spec.read_sumstats import read_sumstats
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class GwasLabLeadVariantsTask(Task):
    """
    A task to generate a list of lead variants from summary statistics.
    Uses Gwaslab.
    see: https://cloufield.github.io/gwaslab/utility_get_lead_novel/
    """

    _sumstats_task: GWASLabCreateSumstatsTask
    short_id: AssetId = field(converter=AssetId)
    sig_level: float = 5e-8

    @property
    def meta(self) -> GWASLabLeadVariantsMeta:
        return GWASLabLeadVariantsMeta(
            trait=self._input_meta.trait,
            project=self._input_meta.project,
            id=self.short_id,
        )

    @property
    def _input_meta(self) -> GWASLabSumStatsMeta:
        input_meta = self._sumstats_task.meta
        assert isinstance(input_meta, GWASLabSumStatsMeta)
        return input_meta

    @property
    def _input_asset_id(self) -> AssetId:
        return self._input_meta.asset_id

    @property
    def deps(self) -> list["Task"]:
        return [self._sumstats_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        sumstats_asset = fetch(self._input_asset_id)
        sumstats: gl.Sumstats = read_sumstats(sumstats_asset)
        variant_df = sumstats.get_lead(anno=True, sig_level=self.sig_level)

        """
        GWASLab (sumstats.get_lead()) does not add the GENE column if no significant variants are found.
        See: https://github.com/Cloufield/gwaslab/blob/eb4daf636f6e171d74a80c16723293734e851ee2/src/gwaslab/util/util_in_get_sig.py#L164-L177. 
        However, the GENE column is required by downstream tasks.
        """
        if not "GENE" in variant_df.columns:
            variant_df["GENE"] = None

        out_path = scratch_dir / "lead_variants.csv"
        variant_df.to_csv(out_path, index=False)
        return FileAsset(out_path)

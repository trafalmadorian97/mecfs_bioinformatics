import structlog

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.task.gwaslab.gwaslab_region_plots_task import (
    plot_region_around_variant,
)

logger = structlog.get_logger()

from pathlib import Path

import gwaslab as gl
from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_region_plots_meta import (
    GWASLabRegionPlotsMeta,
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
from mecfs_bio.build_system.task.gwaslab.gwaslab_util import (
    gwaslab_download_ref_if_missing,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLabVCFRefFile,
)


@frozen
class GwasLabRegionPlotsFromLeadVariantsTask(Task):
    """
    A task to generate region plots around arbitrary locus in genome.

    Useful for visualizing the local significance structure
     see https://cloufield.github.io/gwaslab/tutorial_3.4/#quick-regional-plot-without-ld-information
    Gwaslab can also use a vcf reference file to plot the linkage disequilibrium structure around the lead variants
    (vcf_name_for_lead_variants).  Doing this at a reasonable speed requires the installation of the "tabix" binary.
    """

    _meta: Meta
    _sumstats_task: GWASLabCreateSumstatsTask
    vcf_name_for_ld: GWASLabVCFRefFile | None
    chrom: int
    pos: int
    buffer: int = 500_000

    @property
    def meta(self) -> Meta:
        return self._meta
        # return GWASLabRegionPlotsMeta(
        #     trait=self._lead_variants_task_meta.trait,
        #     project=self._lead_variants_task_meta.project,
        #     short_id=self.short_id,
        # )

    @property
    def _sumstats_meta(self) -> GWASLabSumStatsMeta:
        meta = self._sumstats_task.meta
        assert isinstance(meta, GWASLabSumStatsMeta)
        return meta

    @property
    def _sumstats_id(self) -> AssetId:
        return self._sumstats_meta.asset_id

    @property
    def deps(self) -> list["Task"]:
        return [self._sumstats_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> DirectoryAsset:
        target_path = scratch_dir / "plot_dir"
        target_path.mkdir(parents=True, exist_ok=True)
        sumstats = read_sumstats(fetch(self._sumstats_id))

        plot_region_around_variant(
            sumstats=sumstats,
            vcf_name_for_ld=self.vcf_name_for_ld,
            chrom=self.chrom,
            pos=self.pos,
            buffer=self.buffer,
            output_path=target_path / "region_plot.png",
        )
        return DirectoryAsset(
            path=target_path,
        )

    @classmethod
    def create(
        cls,
        asset_id: str,
        sumstats_task: GWASLabCreateSumstatsTask,
        vcf_name_for_ld: GWASLabVCFRefFile | None,
        chrom: int,
        pos: int,
        buffer: int = 500_000,
    ):
        reference_meta = sumstats_task.meta
        assert isinstance(reference_meta, GWASLabSumStatsMeta)
        meta = GWASLabRegionPlotsMeta(
            trait=reference_meta.trait,
            project=reference_meta.project,
            short_id=AssetId(asset_id),
        )
        return cls(
            meta=meta,
            sumstats_task=sumstats_task,
            vcf_name_for_ld=vcf_name_for_ld,
            chrom=chrom,
            pos=pos,
            buffer=buffer,
        )


def _plot_region_around_variant(
    sumstats: gl.Sumstats,
    chrom: int,
    pos: int,
    buffer: int,
    output_path: Path,
    vcf_name_for_ld: GWASLabVCFRefFile | None,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if vcf_name_for_ld is not None:
        gwaslab_download_ref_if_missing(vcf_name_for_ld)
        vcf_path = gl.get_path(vcf_name_for_ld)
    else:
        vcf_path = None
    scaled = "MLOG10P" in sumstats.data.columns
    sumstats.plot_mqq(
        mode="r",
        skip=2,
        cut=20,
        scaled=scaled,
        region_grid=True,
        region=(chrom, max(pos - buffer, 0), pos + buffer),
        save=str(
            output_path,
        ),
        save_args={"dpi": 400, "facecolor": "white"},
        vcf_path=vcf_path,
    )

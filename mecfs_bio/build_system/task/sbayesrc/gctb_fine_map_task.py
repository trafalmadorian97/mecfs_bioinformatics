"""
Task assembling a RemoteJob that runs GCTB genome-wide fine-mapping (gctb --gwfm RC)
and dispatching it via the WF's remote executor.

"""

import json
from pathlib import Path, PurePath

from attrs import frozen

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.sbayesrc.gctb_gwfm_constants import (
    DEFAULT_DISK_GB,
    DEFAULT_MEMORY_GB,
    DEFAULT_PEP,
    DEFAULT_PIP,
    DEFAULT_VCPUS,
    GCTB_CS_TEMPLATE,
    GCTB_GWFM_TEMPLATE,
    GCTB_MAKE_LDM_EIGEN_TEMPLATE,
    GWFM_ANNOT_FILE_NAME,
    GWFM_ANNOT_ZIP_NAME,
    GWFM_GENE_MAP_FILE_NAME,
    GWFM_LDM_DIR_NAME,
    GWFM_LDM_ZIP_NAME,
    GWFM_OUT_PREFIX,
    GWFM_PWLD_RELPATH,
    MARKER_PREFIX_KEY,
    MATCHED_LDM_OUT,
    REMOTE_MA_PATH,
    REMOTE_OUT_DIR,
    REMOTE_REF_DIR,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.build_system.wf.remote_executor.remote_job import (
    RemoteJob,
    RemoteResources,
)


def _build_commands(threads: int, pip: float, pep: float) -> list[str]:
    """
    Build the ordered gctb container command sequence for genome-wide fine-mapping.

    """
    ref = PurePath(REMOTE_REF_DIR)
    raw_ldm = ref / GWFM_LDM_DIR_NAME
    annot = ref / GWFM_ANNOT_FILE_NAME
    gene_map = ref / GWFM_GENE_MAP_FILE_NAME
    pwld = ref / GWFM_PWLD_RELPATH
    return [
        f"unzip -o {ref / GWFM_LDM_ZIP_NAME} -d {ref}",
        f"unzip -o {ref / GWFM_ANNOT_ZIP_NAME} -d {ref}",
        f"mkdir -p {REMOTE_OUT_DIR}",
        GCTB_MAKE_LDM_EIGEN_TEMPLATE.format(
            ldm=raw_ldm,
            ma=REMOTE_MA_PATH,
            threads=threads,
            out=MATCHED_LDM_OUT,
        ),
        GCTB_GWFM_TEMPLATE.format(
            ldm=MATCHED_LDM_OUT,
            ma=REMOTE_MA_PATH,
            annot=annot,
            gene_map=gene_map,
            threads=threads,
            out=GWFM_OUT_PREFIX,
        ),
        GCTB_CS_TEMPLATE.format(
            pwld=pwld,
            pip=pip,
            pep=pep,
            gene_map=gene_map,
            out=GWFM_OUT_PREFIX,
        ),
    ]


@frozen
class GctbFineMapTask(Task):
    """
    Assembles a RemoteJob running GCTB genome-wide fine-mapping and dispatches it via the
    WF's remote executor, returning the retrieved output directory as a DirectoryAsset.
    """

    meta: Meta
    ma_task: Task
    reference_task: Task
    image: str
    threads: int = DEFAULT_VCPUS
    region: str | None = None
    pip: float = DEFAULT_PIP
    pep: float = DEFAULT_PEP

    @property
    def deps(self) -> list["Task"]:
        return [self.ma_task, self.reference_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> DirectoryAsset:
        ma_asset = fetch(self.ma_task.asset_id)
        assert isinstance(ma_asset, FileAsset)
        marker_asset = fetch(self.reference_task.asset_id)
        assert isinstance(marker_asset, FileAsset)
        marker = json.loads(Path(marker_asset.path).read_text())
        s3_prefix = marker[MARKER_PREFIX_KEY]

        job = RemoteJob(
            image=self.image,
            commands=_build_commands(self.threads, self.pip, self.pep),
            input_files={ma_asset.path: PurePath(REMOTE_MA_PATH)},
            s3_inputs={s3_prefix: PurePath(REMOTE_REF_DIR)},
            output_files=[PurePath(REMOTE_OUT_DIR)],
            resources=RemoteResources(
                DEFAULT_MEMORY_GB, self.threads, DEFAULT_DISK_GB, self.region
            ),
        )
        wf.remote_executor.run(job, scratch_dir)
        return DirectoryAsset(scratch_dir / REMOTE_OUT_DIR)

    @classmethod
    def create(
        cls,
        id: str,
        ma_task: Task,
        reference_task: Task,
        image: str,
        threads: int = DEFAULT_VCPUS,
        region: str | None = None,
        pip: float = DEFAULT_PIP,
        pep: float = DEFAULT_PEP,
    ) -> "GctbFineMapTask":
        ma_meta = ma_task.meta
        assert isinstance(ma_meta, FilteredGWASDataMeta)
        meta = ResultDirectoryMeta(
            id=id,
            trait=ma_meta.trait,
            project=ma_meta.project,
            sub_dir=PurePath("gwfm"),
        )
        return cls(
            meta=meta,
            ma_task=ma_task,
            reference_task=reference_task,
            image=image,
            threads=threads,
            region=region,
            pip=pip,
            pep=pep,
        )

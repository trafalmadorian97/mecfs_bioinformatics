import tempfile
from pathlib import Path

import gwaslab
import structlog

from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.filtered_snps_gwas_1 import (
    DECODE_ME_FILTER_SNPS_GWAS_1_TASK,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
    GwasLabTransformSpec,
    GWASLabVCFRef,
    HarmonizationOptions,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_transform_sumstats import (
    GWASLabTransformSumstatsTask,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.head_pipe import HeadPipe

logger = structlog.get_logger()


def count_ambiguous_or_missing_cases(sumstats: gwaslab.Sumstats) -> int:
    """
    See status codes here: https://cloufield.github.io/gwaslab/StatusCode/
    status code 7 indicates and indistinguishable variant, while status code 8 indicates a variant not in the reference
    """
    ambiguous_or_missing = (
        (sumstats.data["STATUS"].str.slice(-1) == "7")
        | (sumstats.data["STATUS"].str.slice(-1) == "8")
    ).sum()
    return ambiguous_or_missing


_decode_me_select_first_rows = PipeDataFrameTask.create(
    source_task=DECODE_ME_FILTER_SNPS_GWAS_1_TASK,
    asset_id="testing_decode_me_slice_task",
    out_format=ParquetOutFormat(),
    pipes=[HeadPipe(num_rows=10_000)],
)


_decode_me_liftover_to_37_first_rows = GWASLabCreateSumstatsTask(
    df_source_task=_decode_me_select_first_rows,
    asset_id=AssetId("testing_first_rows_decode_me_gwas_1_sumstats_liftover_to_37"),
    basic_check=True,
    genome_build="infer",
    liftover_to="19",
)

_harmonized_task = GWASLabTransformSumstatsTask.create_from_source_task(
    _decode_me_liftover_to_37_first_rows,
    asset_id="testing_first_rows__harmonized",
    spec=GwasLabTransformSpec(
        harmonize_options=HarmonizationOptions(
            ref_infer=GWASLabVCFRef(name="1kg_eur_hg19", ref_alt_freq="AF"),
            ref_seq="ucsc_genome_hg19",
            check_ref_files=True,
            drop_missing_from_ref_seq=True,
            cores=4,
        )
    ),
)


def test_harmonize_drop_ambiguous():
    """
    The purpose of this test is to verify that our harmonization task correctly drops
    palindromic variants whose strand cannot be inferred, and variants not present in the VCF reference
    """
    with tempfile.TemporaryDirectory() as tempdirname:
        tempdir = Path(tempdirname)
        info_store = tempdir / "info_store.yaml"
        asset_root = tempdir / "asset_store"
        asset_root.mkdir(parents=True, exist_ok=True)
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run([_harmonized_task])
        result = result[_harmonized_task.asset_id]
        assert isinstance(result, FileAsset)
        sumstats: gwaslab.Sumstats = gwaslab.load_pickle(str(result.path))
        count = count_ambiguous_or_missing_cases(sumstats)
        logger.debug(f"Result dataframe has shape {sumstats.data.shape}")
        assert count == 0

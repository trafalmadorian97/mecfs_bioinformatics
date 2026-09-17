"""
Run one harmonization scenario in this process and print peak anonymous RSS as JSON.

RssAnon excludes file-backed pages, so the memory-mapped FASTA does not count against
the task. Scenarios:
- genome_reference_all: GenomeReferenceHarmonizationTask on DecodeME build 37
- genome_reference_chr1_2: the same, restricted to chromosomes 1 and 2 by a pipe
- gwaslab: the existing gwaslab harmonization task on the same input

Run (normally via memory_benchmark.py):
  pixi r python -m experiments.claude.genome_reference_harmonization.measure_harmonization_memory genome_reference_all


Historical record (V4): this ran against the pre-switch code, before Task 11 replaced
the gwaslab harmonize_task it imports. It no longer imports; the committed
memory_benchmark.log at that commit is the evidence.
"""

import json
import sys
import tempfile
import threading
import time
from pathlib import Path

from experiments.claude.genome_reference_harmonization.compare_with_gwaslab import CASES
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import DirMeta
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL

SAMPLE_SECONDS = 0.2
_DECODE_ME_PRE_TABLE = CASES[0].pre_harmonization


def _rss_anon_kib() -> int:
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("RssAnon:"):
            return int(line.split()[1])
    raise RuntimeError("RssAnon not reported by /proc/self/status")


class _PeakSampler(threading.Thread):
    def __init__(self) -> None:
        super().__init__(daemon=True)
        self.peak_kib = 0
        self.stopped = threading.Event()

    def run(self) -> None:
        while not self.stopped.is_set():
            self.peak_kib = max(self.peak_kib, _rss_anon_kib())
            time.sleep(SAMPLE_SECONDS)


def _fetch_from_store(task: Task):
    by_id = {dep.asset_id: dep for dep in task.deps}

    def fetch(asset_id: AssetId) -> Asset:
        meta = by_id[asset_id].meta
        path = DEFAULT_RUNNER.meta_to_path(meta)
        return DirectoryAsset(path) if isinstance(meta, DirMeta) else FileAsset(path)

    return fetch


def _scenario_task(name: str) -> Task:
    if name == "gwaslab":
        return DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.harmonize_task
    pipe = (
        FilterRowsByValue(target_column=GWASLAB_CHROM_COL, valid_values=[1, 2])
        if name == "genome_reference_chr1_2"
        else IdentityPipe()
    )
    assert name in {"genome_reference_all", "genome_reference_chr1_2"}, (
        f"unknown scenario {name}"
    )
    return GenomeReferenceHarmonizationTask.create(
        asset_id=f"experiment_memory_{name}",
        sumstats_task=_DECODE_ME_PRE_TABLE,
        fasta_task=UCSC_HG19_INDEXED_FASTA,
        panel_task=THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
        pipe=pipe,
    )


def main() -> None:
    name = sys.argv[1]
    task = _scenario_task(name)
    DEFAULT_RUNNER.run(list(task.deps))  # materialize inputs before sampling
    sampler = _PeakSampler()
    baseline_kib = _rss_anon_kib()
    sampler.start()
    started = time.monotonic()
    with tempfile.TemporaryDirectory() as scratch:
        task.execute(
            scratch_dir=Path(scratch), fetch=_fetch_from_store(task), wf=make_wf()
        )
    elapsed = time.monotonic() - started
    sampler.stopped.set()
    sampler.join()
    print(
        json.dumps(
            {
                "scenario": name,
                "baseline_rss_anon_gib": baseline_kib / 2**20,
                "peak_rss_anon_gib": sampler.peak_kib / 2**20,
                "seconds": elapsed,
            }
        )
    )


if __name__ == "__main__":
    main()

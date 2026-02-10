from pathlib import Path

from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import GENE_THESAURUS
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import ImoHasher
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from test_mecfs_bio.system.util import log_on_error


def go():
    tmp_path = Path("data")/"tmp"

    info_store = tmp_path / "info_store.yaml"
    asset_root = tmp_path / "asset_store"

    with log_on_error(info_store):
        asset_root.mkdir(parents=True, exist_ok=True)
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run(
            [
                GENE_THESAURUS
            ],
            incremental_save=True,
            must_rebuild_transitive=[GENE_THESAURUS]
        )
if __name__ == "__main__":
    go()
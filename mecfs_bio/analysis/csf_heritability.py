from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.csf_pqtl.csf_database.hapmap3.hapmap3_csf_heritability import (
    HAPMAP_3_CSF_HERITABILITY,
)


def go() -> None:
    DEFAULT_RUNNER.run([HAPMAP_3_CSF_HERITABILITY])


if __name__ == "__main__":
    go()

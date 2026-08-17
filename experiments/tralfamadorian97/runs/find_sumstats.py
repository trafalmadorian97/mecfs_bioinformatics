from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import \
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS


def go():
    DEFAULT_RUNNER.run(
        [DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS.join_task]
    )

if __name__ == "__main__":
    go()
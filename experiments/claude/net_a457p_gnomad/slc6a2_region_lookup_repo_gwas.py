"""Look up the SLC6A2 (NET) region in GWAS already materialized in this repo's asset store.

For each GWAS: SNP-level summary (min P and top SNP within SLC6A2 +/- 50 kb and +/- 500 kb, GRCh37), P values at
named SLC6A2 variants (eQTL leads, coding variants, candidate-gene SNPs from the literature), and the MAGMA
gene-level P for SLC6A2 when the gene analysis is materialized.

Assets are read from the main checkout's asset store (configured by its default_runner_config.yaml), so this runs
from a worktree too. Missing assets are reported, not built.

Run: pixi r python experiments/claude/net_a457p_gnomad/slc6a2_region_lookup_repo_gwas.py
"""

from pathlib import Path

import attrs
import polars as pl
import yaml
from scipy.stats import norm

from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_standard_analysis import (
    KEATON_DBP_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.brainstem.whole_brainstem.xue_et_al.analysis.xue_whole_brainstem_standard_analysis import (
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.analysis.standard_analysis_kerrebijin_fibro import (
    KERREBIJN_ET_AL_FIBRO_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.heart_rate_recovery.verweij_et_al_2018.analysis.verweiji_standard_analysis import (
    VERWEIJI_ET_AL_HRR_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_37_parquet_file import (
    DECODE_ME_GWAS_1_LIFTOVER_TO_37_PARQUET,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_ensembl_gene_analysis import (
    DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_female.analysis.decode_me_female_standard_analysis import (
    DECODE_ME_FEMALE_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_gwas_2.analysis.decode_me_gwas_2_standard_analysis import (
    DECODE_ME_GWAS_2_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_infectious_onset.analysis.decode_me_infectious_onset_standard_analysis import (
    DECODE_ME_INFECTIOUS_ONSET_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_male.analysis.decode_me_male_standard_analysis import (
    DECODE_ME_MALE_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.million_veterans.analysis.million_veterans_cfs_standard_analysis import (
    MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP,
)
from mecfs_bio.assets.gwas.me_cfs.neale_lab.analysis.neale_lab_cfs_standard_analysis import (
    NEALE_LAB_CFS_STANDARD_ANALYSIS_TASK_GROUP,
)
from mecfs_bio.assets.gwas.migraine.uk_biobank_2025.analysis.uk_biobank_2025_migraine_standard_analysis import (
    UK_BIOBANK_2025_EUR_MIGRAINE_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import (
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.standard_analysis_sc_pgc_2022 import (
    SCH_PGC_2022_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.analysis.syncope_magma_task_generator import (
    AEGISDOTTIR_COMBINED_MAGMA_TASKS,
)
from mecfs_bio.build_system.rebuilder.metadata_to_path.remapping_meta_to_path import (
    PathRemapRule,
    RemappingMetaToPath,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import (
    read_magma_gene_analysis_result,
)

MAIN_CHECKOUT = Path("/home/paiforsyth/src/traf/biostatistics")
OUT_PATH = Path(__file__).parent / "slc6a2_region_lookup_repo_gwas.tsv"

SLC6A2_CHROM = 16
SLC6A2_START_37 = 55_689_516
SLC6A2_END_37 = 55_740_102
"""GRCh37 gene body (GRCh38 55,655,604-55,706,192 from GTEx/GENCODE, minus the 33,912 bp b37->b38 offset here)."""
NEAR_FLANK = 50_000
WIDE_FLANK = 500_000
SLC6A2_GENE_IDS = ("6530", "ENSG00000103546")
"""MAGMA gene ids: Entrez and Ensembl, depending on the gene-location file used."""


@attrs.frozen(slots=True)
class NamedVariant:
    rsid: str
    pos_37: int
    """GRCh37 position; matched by position because some tables (DecodeME) carry no rsID column."""
    note: str


NAMED_VARIANTS = (
    NamedVariant(
        "rs55830311", 55_682_671, "eQTL lead esophagus; height/lean-mass coloc"
    ),
    NamedVariant("rs28386840", 55_686_818, "-3081A/T promoter, ADHD/MPH candidate"),
    NamedVariant("rs2242446", 55_690_425, "-182T/C promoter, ADHD candidate"),
    NamedVariant(
        "rs9937434", 55_693_855, "eQTL lead adrenal/skin; alt C lowers expression"
    ),
    NamedVariant("rs17841327", 55_694_253, "eQTL lead skin; alt C lowers expression"),
    NamedVariant("rs1805065", 55_703_498, "T99I"),
    NamedVariant("rs552629380", 55_705_760, "start lost"),
    NamedVariant("rs11568343", 55_705_838, "C27Y"),
    NamedVariant("rs187714", 55_706_499, "eQTL lead testis; birth-weight coloc"),
    NamedVariant("rs3785151", 55_712_519, "MVP sensorineural hearing loss lead"),
    NamedVariant("rs1805066", 55_719_143, "V245I"),
    NamedVariant("rs12920735", 55_725_394, "Keaton DBP lead in SLC6A2"),
    NamedVariant("rs5564", 55_725_975, "Bayles 2012 POTS-associated synonymous SNP"),
    NamedVariant("rs5569", 55_731_835, "G1287A, ADHD/MPH candidate; GLGC HDL lead"),
    NamedVariant("rs121918126", 55_731_917, "A457P site (also A457T/S)"),
    NamedVariant("rs42460", 55_737_656, "tuQTL testis; placental-weight coloc"),
)


@attrs.frozen(slots=True)
class GwasSource:
    label: str
    snp_table_task: Task
    magma_gene_analysis_task: Task | None


def from_standard_analysis(label: str, group: object) -> GwasSource:
    """Accepts both StandardAnalysisTaskGroup and StandardAnalysisTaskGroupAddRSIDS."""
    inner = getattr(group, "tasks", group)
    magma_tasks = inner.magma_tasks
    return GwasSource(
        label=label,
        snp_table_task=magma_tasks.parquet_file_task,
        magma_gene_analysis_task=magma_tasks.inner.gene_analysis_task,
    )


SOURCES = (
    GwasSource(
        "DecodeME GWAS-1 (ME/CFS)",
        DECODE_ME_GWAS_1_LIFTOVER_TO_37_PARQUET,
        DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
    ),
    from_standard_analysis(
        "DecodeME GWAS-2 (ME/CFS)", DECODE_ME_GWAS_2_STANDARD_ANALYSIS
    ),
    from_standard_analysis(
        "DecodeME female (ME/CFS)", DECODE_ME_FEMALE_STANDARD_ANALYSIS
    ),
    from_standard_analysis("DecodeME male (ME/CFS)", DECODE_ME_MALE_STANDARD_ANALYSIS),
    from_standard_analysis(
        "DecodeME infectious onset (ME/CFS)",
        DECODE_ME_INFECTIOUS_ONSET_STANDARD_ANALYSIS,
    ),
    from_standard_analysis(
        "Million Veterans CFS", MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP
    ),
    from_standard_analysis(
        "Neale lab UKB CFS", NEALE_LAB_CFS_STANDARD_ANALYSIS_TASK_GROUP
    ),
    from_standard_analysis(
        "Heart rate recovery (Verweij 2018)", VERWEIJI_ET_AL_HRR_STANDARD_ANALYSIS
    ),
    GwasSource(
        "Syncope (Aegisdottir)",
        AEGISDOTTIR_COMBINED_MAGMA_TASKS.parquet_file_task,
        AEGISDOTTIR_COMBINED_MAGMA_TASKS.inner.gene_analysis_task,
    ),
    from_standard_analysis(
        "Diastolic blood pressure (Keaton)", KEATON_DBP_STANDARD_ANALYSIS
    ),
    from_standard_analysis(
        "Multisite chronic pain (Johnston)", JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS
    ),
    from_standard_analysis(
        "Fibromyalgia (Kerrebijn)", KERREBIJN_ET_AL_FIBRO_STANDARD_ANALYSIS
    ),
    from_standard_analysis(
        "Migraine (UKB 2025 EUR)", UK_BIOBANK_2025_EUR_MIGRAINE_STANDARD_ANALYSIS
    ),
    from_standard_analysis("Schizophrenia (PGC 2022)", SCH_PGC_2022_STANDARD_ANALYSIS),
    from_standard_analysis(
        "Whole brainstem volume (Xue)", XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS
    ),
)


def main_checkout_meta_to_path() -> RemappingMetaToPath:
    config = yaml.safe_load((MAIN_CHECKOUT / "default_runner_config.yaml").read_text())
    return RemappingMetaToPath(
        default_root=MAIN_CHECKOUT / config["asset_root"],
        rules=PathRemapRule.tuple_from_config(config.get("path_remap", [])),
    )


def parquet_file_in(path: Path) -> Path:
    if path.is_file():
        return path
    files = sorted(path.glob("*.parquet"))
    assert len(files) == 1, f"expected one parquet file in {path}, found {files}"
    return files[0]


def region_snps(table_path: Path) -> pl.DataFrame:
    frame = pl.scan_parquet(parquet_file_in(table_path))
    columns = frame.collect_schema().names()
    keep = [
        c
        for c in (
            "SNPID",
            "rsID",
            "CHR",
            "POS",
            "EA",
            "NEA",
            "EAF",
            "BETA",
            "SE",
            "P",
            "MLOG10P",
            "N",
        )
        if c in columns
    ]
    snps = (
        frame.select(keep)
        .filter(
            (pl.col("CHR").cast(pl.Int64) == SLC6A2_CHROM)
            & pl.col("POS").is_between(
                SLC6A2_START_37 - WIDE_FLANK, SLC6A2_END_37 + WIDE_FLANK
            )
        )
        .collect()
    )
    return with_p_value(snps)


def with_p_value(snps: pl.DataFrame) -> pl.DataFrame:
    """Some tables carry only MLOG10P or only BETA/SE; derive a two-sided P in that case."""
    if "P" in snps.columns:
        return snps
    if "MLOG10P" in snps.columns:
        return snps.with_columns((10.0 ** -pl.col("MLOG10P")).alias("P"))
    assert {"BETA", "SE"} <= set(snps.columns), f"cannot derive P from {snps.columns}"
    z = (snps["BETA"] / snps["SE"]).abs().to_numpy()
    return snps.with_columns(pl.Series("P", 2 * norm.sf(z)))


def window_summary(snps: pl.DataFrame, flank: int) -> dict:
    window = snps.filter(
        pl.col("POS").is_between(SLC6A2_START_37 - flank, SLC6A2_END_37 + flank)
    )
    if window.is_empty():
        return {"n_snps": 0}
    top = window.sort("P").row(0, named=True)
    return {
        "n_snps": window.height,
        "min_p": top["P"],
        "top_snp": top.get("rsID") or top.get("SNPID"),
        "top_pos": top["POS"],
        "top_beta": top.get("BETA"),
    }


def magma_slc6a2(result_dir: Path) -> dict:
    genes = pl.from_pandas(read_magma_gene_analysis_result(result_dir))
    row = genes.filter(pl.col("GENE").cast(pl.Utf8).is_in(SLC6A2_GENE_IDS))
    if row.is_empty():
        return {"magma_p": None}
    hit = row.row(0, named=True)
    return {"magma_p": hit["P"], "magma_z": hit["ZSTAT"], "magma_n_genes": genes.height}


def check_build_37(snps: pl.DataFrame) -> bool | None:
    """rs9937434 is at 55,693,855 on GRCh37 and 55,659,943 on GRCh38."""
    if "rsID" not in snps.columns:
        return None
    marker = snps.filter(pl.col("rsID") == "rs9937434")
    if marker.is_empty():
        return None
    return marker["POS"][0] == 55_693_855


def summarize_source(
    source: GwasSource, meta_to_path: RemappingMetaToPath
) -> list[dict]:
    table_path = meta_to_path(source.snp_table_task.meta)
    if not table_path.exists():
        print(f"MISSING {source.label}: {table_path}")
        return [{"gwas": source.label, "status": "not materialized"}]
    snps = region_snps(table_path)
    base = {"gwas": source.label, "status": "ok", "is_build_37": check_build_37(snps)}
    base |= {f"near_{k}": v for k, v in window_summary(snps, NEAR_FLANK).items()}
    base |= {f"wide_{k}": v for k, v in window_summary(snps, WIDE_FLANK).items()}
    if source.magma_gene_analysis_task is not None:
        magma_path = meta_to_path(source.magma_gene_analysis_task.meta)
        base |= (
            magma_slc6a2(magma_path)
            if magma_path.exists()
            else {"magma_p": "not materialized"}
        )
    rows = [base]
    by_position = {v.pos_37: v for v in NAMED_VARIANTS}
    for row in snps.filter(pl.col("POS").is_in(list(by_position))).iter_rows(
        named=True
    ):
        named = by_position[row["POS"]]
        rows.append(
            {
                "gwas": source.label,
                "status": "named_variant",
                "variant": named.rsid,
                "variant_note": named.note,
                "EA": row.get("EA"),
                "NEA": row.get("NEA"),
                "EAF": row.get("EAF"),
                "BETA": row.get("BETA"),
                "SE": row.get("SE"),
                "P": row["P"],
            }
        )
    return rows


def main() -> None:
    meta_to_path = main_checkout_meta_to_path()
    rows = [row for source in SOURCES for row in summarize_source(source, meta_to_path)]
    result = pl.DataFrame(rows, infer_schema_length=None)
    result.write_csv(OUT_PATH, separator="\t")
    with pl.Config(tbl_rows=400, tbl_cols=30, fmt_str_lengths=40, tbl_width_chars=250):
        print(
            result.filter(pl.col("status") != "named_variant").drop(
                "variant",
                "variant_note",
                "EA",
                "NEA",
                "EAF",
                "BETA",
                "SE",
                "P",
                strict=False,
            )
        )
        print(
            result.filter(pl.col("status") == "named_variant")
            .select("gwas", "variant", "variant_note", "EA", "EAF", "BETA", "P")
            .sort("P")
        )


if __name__ == "__main__":
    main()

"""
Compute local genetic correlation using LAVA (Local Analysis of [co]Variant Association).


Implemented by Github copilot

See: Werme et al. "An integrated framework for local genetic correlation analysis"
Nature Genetics 54 (2022): 274-282.

LAVA estimates local heritability (h2) for each phenotype at each genomic locus,
and local genetic correlation (rg) between pairs of phenotypes at loci where both
show significant heritability.

This task wraps the R LAVA package via rpy2.
"""

import gc
import tempfile
from pathlib import Path, PurePath
from typing import Sequence

import numpy as np
import pandas as pd
import polars as pl
import rpy2.robjects as ro
import structlog
from attrs import frozen
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from tqdm import tqdm

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

logger = structlog.get_logger()

# LAVA expected column names
LAVA_SNP_COL = "SNP"
LAVA_A1_COL = "A1"
LAVA_A2_COL = "A2"
LAVA_N_COL = "N"
LAVA_Z_COL = "Z"

# Output filenames
UNIV_RESULTS_FILENAME = "lava_univariate.csv"
BIVAR_RESULTS_FILENAME = "lava_bivariate.csv"


@frozen
class LavaBinarySampleInfo:
    """
    Sample information for a binary (case/control) phenotype.

    LAVA uses cases and controls to determine if a phenotype is binary,
    and optionally uses prevalence to convert observed h2 to liability scale.
    """

    cases: int
    controls: int
    prevalence: float | None = None


@frozen
class LavaContinuousSampleInfo:
    """Marker indicating a continuous (quantitative) phenotype."""

    pass


LavaSampleInfo = LavaBinarySampleInfo | LavaContinuousSampleInfo


@frozen
class LavaPhenotypeDataSource:
    """
    A source Task providing tabular GWAS summary statistics for a phenotype.

    Column names are assumed to be in GWASLAB format
    (see: mecfs_bio/constants/gwaslab_constants.py).
    Data will be copied and converted to the format expected by LAVA.
    """

    task: Task
    alias: str
    pipe: DataProcessingPipe = IdentityPipe()
    sample_info: LavaSampleInfo = LavaContinuousSampleInfo()

    @property
    def asset_id(self) -> AssetId:
        return self.task.asset_id


@frozen
class LDReferenceInfo:
    """
    Reference LD data for LAVA.

    filename_prefix is the prefix of the PLINK-format LD reference files
    (i.e. the path without .bed/.bim/.fam extensions).
    """

    ld_ref_task: Task
    filename_prefix: str


@frozen
class LavaTask(Task):
    """
    Given a locus definition file, does the following for each locus:
    - Estimates heritability for all phenotypes using LAVA
    - For each pair of phenotypes whose heritabilities are both significant
      at the locus, calculates local genetic correlation using LAVA

    Outputs two CSV files to scratch_dir:
    - lava_univariate.csv: heritability estimates per phenotype per locus
    - lava_bivariate.csv: genetic correlation estimates per phenotype pair per locus
    """

    _meta: Meta
    sources: Sequence[LavaPhenotypeDataSource]
    ld_reference_info: LDReferenceInfo
    lava_locus_definitions_task: Task
    ct_ldsc_task_for_overlap: Task | None = None
    univ_p_threshold: float = 0.05
    max_loci: int | None = None

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        result: list[Task] = [item.task for item in self.sources]
        result.append(self.ld_reference_info.ld_ref_task)
        result.append(self.lava_locus_definitions_task)
        if self.ct_ldsc_task_for_overlap is not None:
            result.append(self.ct_ldsc_task_for_overlap)
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        lava = importr("LAVA")
        conv = ro.default_converter + pandas2ri.converter
        r_is_null = ro.r("is.null")

        # Fetch LD reference
        ld_ref_asset = fetch(self.ld_reference_info.ld_ref_task.asset_id)
        assert isinstance(ld_ref_asset, DirectoryAsset)
        ref_prefix = str(ld_ref_asset.path / self.ld_reference_info.filename_prefix)

        # Fetch locus definitions
        locus_asset = fetch(self.lava_locus_definitions_task.asset_id)
        assert isinstance(locus_asset, FileAsset)

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)

            # Write LAVA-compatible summary statistics and build input info
            info_rows: list[dict] = []
            pheno_names: list[str] = []
            for source in self.sources:
                sumstats_path = _write_lava_sumstats(
                    source=source, fetch=fetch, tmp_dir=tmp_path
                )
                info_rows.append(_make_info_row(source, sumstats_path))
                pheno_names.append(source.alias)

            # Write the input info file
            info_file_path = tmp_path / "input.info.txt"
            info_df = pd.DataFrame(info_rows)
            info_df.to_csv(info_file_path, sep="\t", index=False)
            logger.debug(f"LAVA input info file:\n{info_df}")

            # Build sample overlap matrix from CT-LDSC results if available
            sample_overlap_path = _get_sample_overlap_path(
                self.ct_ldsc_task_for_overlap, fetch, pheno_names, tmp_path
            )

            # Process input through LAVA
            sample_overlap_r = (
                ro.NULL if sample_overlap_path is None else sample_overlap_path
            )
            lava_input = lava.process_input(
                input_info_file=str(info_file_path),
                sample_overlap_file=sample_overlap_r,
                ref_prefix=ref_prefix,
                phenos=ro.StrVector(pheno_names),
            )

            # Read loci
            loci_r = lava.read_loci(str(locus_asset.path))
            ro.globalenv["lava_loci"] = loci_r
            n_loci = int(ro.r("nrow(lava_loci)")[0])  # type: ignore
            if self.max_loci is not None:
                n_loci = min(n_loci, self.max_loci)
            logger.debug(f"Processing {n_loci} loci")

            # Process each locus
            univ_results: list[pd.DataFrame] = []
            bivar_results: list[pd.DataFrame] = []

            ro.globalenv["lava_input"] = lava_input
            for i in tqdm(range(1, n_loci + 1)):
                locus_row = ro.r(f"lava_loci[{i},]")
                locus = lava.process_locus(locus_row, lava_input)

                if bool(r_is_null(locus)[0]):  # type: ignore
                    continue

                loc_info = _extract_locus_info(locus)

                # Run univariate + bivariate with automatic filtering
                result = lava.run_univ_bivar(locus, univ_thresh=self.univ_p_threshold)

                # Extract univariate results
                univ_r = result.rx2("univ")
                with localconverter(conv):
                    univ_df: pd.DataFrame = ro.conversion.get_conversion().rpy2py(
                        univ_r
                    )
                for key, val in loc_info.items():
                    univ_df[key] = val
                univ_results.append(univ_df)

                # Extract bivariate results if present
                bivar_r = result.rx2("bivar")
                if not bool(r_is_null(bivar_r)[0]):  # type: ignore
                    with localconverter(conv):
                        bivar_df: pd.DataFrame = ro.conversion.get_conversion().rpy2py(
                            bivar_r
                        )
                    for key, val in loc_info.items():
                        bivar_df[key] = val
                    bivar_results.append(bivar_df)

        _write_output(scratch_dir, univ_results, bivar_results)

        # Clean up R memory
        gc.collect()
        ro.r("gc()")

        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        sources: Sequence[LavaPhenotypeDataSource],
        ld_reference_info: LDReferenceInfo,
        lava_locus_definitions_task: Task,
        ct_ldsc_task_for_overlap: Task | None = None,
        univ_p_threshold: float = 0.05,
    ) -> "LavaTask":
        meta = ResultDirectoryMeta(
            id=asset_id,
            trait="multi_trait",
            project="lava_local_rg",
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            sources=sources,
            ld_reference_info=ld_reference_info,
            lava_locus_definitions_task=lava_locus_definitions_task,
            ct_ldsc_task_for_overlap=ct_ldsc_task_for_overlap,
            univ_p_threshold=univ_p_threshold,
        )


def _extract_locus_info(locus: ro.RObject) -> dict:
    """Extract locus metadata from a processed LAVA locus object.

    LAVA's process.locus() returns an R environment (not a list),
    so we use ro.r("$") to access fields reliably.
    """

    return {
        "locus": int(locus["id"][0]),  # type: ignore
        "chr": int(locus["chr"][0]),  # type: ignore
        "start": int(locus["start"][0]),  # type: ignore
        "stop": int(locus["stop"][0]),  # type: ignore
        "n.snps": int(locus["n.snps"][0]),  # type: ignore
        "n.pcs": int(locus["K"][0]),  # type: ignore
    }


def _write_lava_sumstats(
    source: LavaPhenotypeDataSource, fetch: Fetch, tmp_dir: Path
) -> Path:
    """
    Read gwaslab-format summary statistics and write a LAVA-compatible copy.

    LAVA expects columns: SNP, A1, A2, N, Z
    gwaslab provides: rsID, EA, NEA, N, BETA, SE
    Z-scores are computed as BETA / SE.
    """
    asset = fetch(source.asset_id)
    df = (
        source.pipe.process(scan_dataframe_asset(asset, source.task.meta))
        .collect()
        .to_polars()
    )

    required_cols = [
        GWASLAB_RSID_COL,
        GWASLAB_EFFECT_ALLELE_COL,
        GWASLAB_NON_EFFECT_ALLELE_COL,
        GWASLAB_SAMPLE_SIZE_COLUMN,
        GWASLAB_BETA_COL,
        GWASLAB_SE_COL,
    ]
    for col in required_cols:
        assert col in df.columns, (
            f"Missing required column '{col}' in source '{source.alias}'"
        )

    lava_df = df.select(
        [
            pl.col(GWASLAB_RSID_COL).alias(LAVA_SNP_COL),
            pl.col(GWASLAB_EFFECT_ALLELE_COL).alias(LAVA_A1_COL),
            pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).alias(LAVA_A2_COL),
            pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).alias(LAVA_N_COL),
            (pl.col(GWASLAB_BETA_COL) / pl.col(GWASLAB_SE_COL)).alias(LAVA_Z_COL),
        ]
    )

    output_path = tmp_dir / f"{source.alias}.sumstats.txt"
    lava_df.to_pandas().to_csv(output_path, sep="\t", index=False)
    logger.debug(
        f"Wrote LAVA sumstats for '{source.alias}' to {output_path} ({len(lava_df)} variants)"
    )
    return output_path


def _make_info_row(source: LavaPhenotypeDataSource, sumstats_path: Path) -> dict:
    """Create a row for the LAVA input info file."""
    info = source.sample_info
    if isinstance(info, LavaBinarySampleInfo):
        return {
            "phenotype": source.alias,
            "cases": info.cases,
            "controls": info.controls,
            "prevalence": info.prevalence if info.prevalence is not None else "NA",
            "filename": str(sumstats_path),
        }
    return {
        "phenotype": source.alias,
        "cases": "NA",
        "controls": "NA",
        "prevalence": "NA",
        "filename": str(sumstats_path),
    }


def _get_sample_overlap_path(
    ct_ldsc_task: Task | None,
    fetch: Fetch,
    pheno_names: list[str],
    tmp_dir: Path,
) -> str | None:
    """Build a sample overlap file from CT-LDSC results for LAVA.

    Reads the cross-trait LDSC output (a CSV with columns p1, p2, gcov_int,
    h2_int, …), builds a covariance matrix from the genetic covariance
    intercepts, standardises it via cov2cor, and writes it in the format
    expected by LAVA's ``process.input(sample_overlap_file=…)``.

    The procedure follows the LAVA sample-overlap tutorial:
    https://github.com/josefin-werme/LAVA/blob/main/vignettes/sample_overlap.md

    When self-LDSC entries (p1 == p2) are absent, the diagonal of the
    covariance matrix is filled from the ``h2_int`` column of rows where the
    phenotype appears as ``p2``.
    """
    if ct_ldsc_task is None:
        return None

    ct_ldsc_asset = fetch(ct_ldsc_task.asset_id)
    ct_ldsc_df = (
        scan_dataframe_asset(ct_ldsc_asset, ct_ldsc_task.meta).collect().to_pandas()
    )

    n = len(pheno_names)
    mat = np.full((n, n), np.nan)
    name_to_idx = {name: i for i, name in enumerate(pheno_names)}

    # Fill the matrix from gcov_int (including diagonal for self-LDSC entries)
    for _, row in ct_ldsc_df.iterrows():
        p1, p2 = str(row["p1"]), str(row["p2"])
        if p1 in name_to_idx and p2 in name_to_idx:
            i, j = name_to_idx[p1], name_to_idx[p2]
            mat[i, j] = float(row["gcov_int"])
            if np.isnan(mat[j, i]):
                mat[j, i] = float(row["gcov_int"])

    # Fill diagonal from h2_int when self-LDSC entries are missing
    for _, row in ct_ldsc_df.iterrows():
        p2 = str(row["p2"])
        if p2 in name_to_idx:
            idx = name_to_idx[p2]
            if np.isnan(mat[idx, idx]):
                mat[idx, idx] = float(row["h2_int"])

    # Symmetrise (prefer upper-triangle values)
    for i in range(n):
        for j in range(i + 1, n):
            if np.isnan(mat[i, j]) and not np.isnan(mat[j, i]):
                mat[i, j] = mat[j, i]
            elif np.isnan(mat[j, i]) and not np.isnan(mat[i, j]):
                mat[j, i] = mat[i, j]

    # Standardise using cov2cor: D^{-1/2} C D^{-1/2}
    d = np.sqrt(np.diag(mat))
    mat_cor = mat / np.outer(d, d)
    mat_cor = np.round(mat_cor, 5)

    overlap_path = tmp_dir / "sample_overlap.txt"
    overlap_df = pd.DataFrame(mat_cor, index=pheno_names, columns=pheno_names)
    overlap_df.to_csv(overlap_path, sep=" ")
    logger.debug(f"Sample overlap matrix:\n{overlap_df}")
    return str(overlap_path)


def _write_output(
    scratch_dir: Path,
    univ_results: list[pd.DataFrame],
    bivar_results: list[pd.DataFrame],
) -> None:
    """Write univariate and bivariate results to CSV files."""
    if univ_results:
        all_univ = pd.concat(univ_results, ignore_index=True)
        logger.debug(f"Univariate results ({len(all_univ)} rows):\n{all_univ}")
        all_univ.to_csv(scratch_dir / UNIV_RESULTS_FILENAME, index=False)
    else:
        logger.warning("No univariate results produced")
        pd.DataFrame().to_csv(scratch_dir / UNIV_RESULTS_FILENAME, index=False)

    if bivar_results:
        all_bivar = pd.concat(bivar_results, ignore_index=True)
        logger.debug(f"Bivariate results ({len(all_bivar)} rows):\n{all_bivar}")
        all_bivar.to_csv(scratch_dir / BIVAR_RESULTS_FILENAME, index=False)
    else:
        logger.debug("No bivariate results produced")
        pd.DataFrame().to_csv(scratch_dir / BIVAR_RESULTS_FILENAME, index=False)

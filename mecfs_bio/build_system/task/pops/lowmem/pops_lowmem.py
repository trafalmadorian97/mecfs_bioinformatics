"""
Low-peak-memory reimplementation of POPs' ridge-regression fit.

Stock POPs (pops.py) materializes the full genes-by-features matrix and runs an
SVD of it inside RidgeCV, so peak memory grows with the number of selected
features and out-of-memory-kills the box for high-feature traits (for example
uncapped IBD selects ~14,564 features and needs ~14.6 GB). This script keeps the
entire stock front-end (reading MAGMA, regularizing the error covariance,
projecting out covariates, computing marginal associations, and selecting
features) but replaces only the training block with a streaming kernel-ridge fit
whose peak memory is O(n_genes squared), independent of the feature count.

The exact ridge identity used is the push-through identity, valid for all n, p:

    coef = (X'X + a I_p)^-1 X'y = X' (X X' + a I_n)^-1 y = X' (K + a I_n)^-1 y

with K = X X' the gene-by-gene Gram matrix. Because the munged feature matrix is
stored as column chunks, K accumulates additively across chunks in a single
streaming pass, so we never hold the full X or its SVD workspace. Alpha selection
reproduces scikit-learn's RidgeCV generalized-cross-validation objective in its
eigen mode, which itself operates on K = X X', so the selected alpha and the
resulting coefficients match stock RidgeCV to numerical precision. Only the ridge
method admits this exact dual, so this script supports method=ridge only.

This script is not run directly against the POPs source directory; it imports the
pinned stock pops.py as a module (via --pops_source_dir) and reuses its functions,
so the front-end is literally the same validated code. It is invoked by
PopsLowMemRunTask and accepts the same arguments as pops.py plus --pops_source_dir.
"""

import argparse
import importlib
import logging
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.linalg

# The 25 ridge penalties stock POPs uses (initialize_regressor, method="ridge").
RIDGE_ALPHAS = np.logspace(-2, 10, num=25)


def extract_source_dir(argv: list[str]) -> tuple[Path, list[str]]:
    """Pull --pops_source_dir off the argument vector, returning it and the
    remaining arguments (which are forwarded to stock pops.py's parser)."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--pops_source_dir", required=True)
    namespace, rest = parser.parse_known_args(argv)
    return Path(namespace.pops_source_dir), rest


def build_training_transforms(
    y_proj, y_ids, error_cov, gene_annot_df, training_y_gene_inds, stock
):
    """Reproduce the target-side setup of stock build_training without touching X.

    Returns the training gene ids, the whitened intercept-projected training Y, the
    whitened intercept, and the block-Cholesky-inverse whitening operator plus its
    chromosome block labels (both None when there is no error covariance). These
    are exactly the transforms stock build_training applies to X, factored out so
    they can be applied per feature chunk.
    """
    training_genes = y_ids[training_y_gene_inds]
    sub_y = y_proj[training_y_gene_inds].astype(np.float64)
    intercept = np.ones((sub_y.shape[0], 1))
    linv = None
    labels = None
    if error_cov is not None:
        sub_error_cov = error_cov[np.ix_(training_y_gene_inds, training_y_gene_inds)]
        labels = gene_annot_df.loc[training_genes].CHR.values
        linv = stock.block_Linv(sub_error_cov, labels)
        sub_y = stock.block_AB(linv, labels, sub_y)
        intercept = stock.block_AB(linv, labels, intercept)
    # project_out_intercept=True in stock build_training
    sub_y = stock.project_out_V(sub_y.reshape(-1, 1), intercept).flatten()
    return training_genes, sub_y, intercept, linv, labels


def transform_train_chunk(x_chunk, intercept, linv, labels, stock):
    """Apply stock build_training's X transforms (error-covariance whitening then
    intercept projection) to one feature chunk's training rows. These are
    left-multiplications on the row space, so applying them per column chunk and
    summing the chunk Gram contributions yields the same K as transforming the
    full X at once."""
    if linv is not None:
        x_chunk = stock.block_AB(linv, labels, x_chunk)
    x_chunk = stock.project_out_V(x_chunk, intercept)
    return x_chunk


def stream_build_gram(
    feature_mat_prefix,
    num_feature_chunks,
    selected_features,
    training_genes,
    intercept,
    linv,
    labels,
    stock,
):
    """Accumulate the training gene-by-gene Gram matrix K = X X' in a single pass
    over the munged feature chunks, holding only one chunk plus K at a time.

    X here is the whitened, intercept-projected training feature matrix restricted
    to the selected features, exactly as stock build_training would build it.
    """
    selected_set = set(selected_features)
    rows = np.loadtxt(f"{feature_mat_prefix}.rows.txt", dtype=str).flatten()
    x_train_inds = stock.get_indices_in_target_order(rows, training_genes)
    n_train = len(training_genes)
    gram = np.zeros((n_train, n_train), dtype=np.float64)
    for i in range(num_feature_chunks):
        cols = np.loadtxt(f"{feature_mat_prefix}.cols.{i}.txt", dtype=str).flatten()
        keep = np.array([c in selected_set for c in cols])
        if not keep.any():
            continue
        mat = np.load(f"{feature_mat_prefix}.mat.{i}.npy")
        x_chunk = mat[x_train_inds][:, keep]
        x_chunk = transform_train_chunk(x_chunk, intercept, linv, labels, stock)
        gram += x_chunk @ x_chunk.T
    return gram, rows, x_train_inds


def ridge_gcv_from_gram(gram, y, alphas):
    """Select the ridge penalty by scikit-learn's generalized-cross-validation and
    return the chosen alpha, its GCV score, and the dual coefficients (K + a I)^-1 y.

    This reproduces sklearn RidgeCV(fit_intercept=False, scoring=None) in its eigen
    mode: for each alpha the leave-one-out squared errors are (c / diag(G^-1))^2
    where c = (K + a I)^-1 y and diag(G^-1)_i = sum_k Q[i,k]^2 / (eigval_k + a); the
    score is the negative mean of those errors, and the best alpha maximizes it
    (strict comparison, so the first/smallest alpha wins ties, matching sklearn's
    iteration order).
    """
    eigvals, q_mat = scipy.linalg.eigh(gram)
    qt_y = q_mat.T @ y
    q_sq = q_mat**2
    best_alpha = None
    best_score = None
    best_c = None
    for alpha in alphas:
        w = 1.0 / (eigvals + alpha)
        c = q_mat @ (w * qt_y)
        g_inverse_diag = (w * q_sq).sum(axis=1)
        score = -np.mean((c / g_inverse_diag) ** 2)
        if best_score is None or score > best_score:
            best_alpha, best_score, best_c = float(alpha), float(score), c
    return best_alpha, best_score, best_c


def stream_coefs_and_preds(
    feature_mat_prefix,
    num_feature_chunks,
    selected_features,
    rows,
    x_train_inds,
    intercept,
    linv,
    labels,
    dual_c,
    stock,
):
    """Second streaming pass: build the per-feature coefficients and all-gene PoPS
    scores without materializing the full feature matrix.

    For each chunk, the primal coefficients for its selected columns are
    X_chunk_whitened' dual_c (equivalent to stock RidgeCV.coef_ for those columns),
    and the raw (un-whitened) all-gene columns times those coefficients accumulate
    into the predictions, matching stock pops_predict (mat.dot(coef)). Columns are
    emitted in chunk-then-file order, identical to stock load_feature_matrix.
    """
    selected_set = set(selected_features)
    all_cols = []
    all_betas = []
    preds = np.zeros(len(rows), dtype=np.float64)
    for i in range(num_feature_chunks):
        cols = np.loadtxt(f"{feature_mat_prefix}.cols.{i}.txt", dtype=str).flatten()
        keep = np.array([c in selected_set for c in cols])
        if not keep.any():
            continue
        mat = np.load(f"{feature_mat_prefix}.mat.{i}.npy")
        mat_sel = mat[:, keep]
        cols_sel = cols[keep]
        x_chunk = transform_train_chunk(
            mat_sel[x_train_inds], intercept, linv, labels, stock
        )
        beta_chunk = x_chunk.T @ dual_c
        preds += mat_sel @ beta_chunk
        all_cols.append(cols_sel)
        all_betas.append(beta_chunk)
    cols = np.hstack(all_cols)
    betas = np.hstack(all_betas)
    return cols, betas, preds


def build_coefs_df(cols, betas, alpha, best_score):
    """Assemble the .coefs data frame in stock compute_coefficients' ridge format:
    METHOD/SELECTED_CV_ALPHA/BEST_CV_SCORE metadata rows followed by per-feature
    betas, indexed by parameter name."""
    coefs_df = pd.DataFrame(
        [
            ["METHOD", "RidgeCV"],
            ["SELECTED_CV_ALPHA", alpha],
            ["BEST_CV_SCORE", best_score],
        ]
    )
    coefs_df = pd.concat([coefs_df, pd.DataFrame([cols, betas]).T])
    coefs_df.columns = ["parameter", "beta"]
    coefs_df = coefs_df.set_index("parameter")
    return coefs_df


def select_features(config_dict, marginal_assoc_df, stock):
    """Reproduce stock main's feature-selection annotation of marginal_assoc_df and
    return the selected feature names. Forward stepwise selection is intentionally
    unsupported here (it is never used by our pipeline)."""
    assert config_dict["feature_selection_fss_num_features"] is None, (
        "Forward stepwise selection is not supported by the low-memory POPs "
        "implementation."
    )
    selected_features = stock.select_features_from_marginal_assoc_df(
        marginal_assoc_df,
        config_dict["subset_features_path"],
        config_dict["control_features_path"],
        config_dict["feature_selection_p_cutoff"],
        config_dict["feature_selection_max_num"],
    )
    marginal_assoc_df["selected"] = marginal_assoc_df.index.isin(selected_features)
    marginal_assoc_df["selected"] = marginal_assoc_df["selected"] & ~pd.isnull(
        marginal_assoc_df.pval
    )
    return marginal_assoc_df[marginal_assoc_df.selected].index.values


def main_lowmem(config_dict, stock):
    """Run POPs with the streaming kernel-ridge fit. Mirrors stock pops.py main but
    replaces the training block (load_feature_matrix, build_training,
    compute_coefficients, pops_predict) with the streaming routines above."""
    assert config_dict["method"] == "ridge", (
        "The low-memory POPs implementation only supports method=ridge (the exact "
        "kernel dual does not apply to lasso or linreg)."
    )

    if config_dict["verbose"]:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.info("Verbose output enabled.")
    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")
    np.random.seed(config_dict["random_seed"])
    random.seed(config_dict["random_seed"])
    logging.info(f"Low-memory POPs. Config dict = {str(config_dict)}")

    # --------------------------- Read/process data (stock front-end) --------------------------- #
    gene_annot_df = stock.read_gene_annot_df(config_dict["gene_annot_path"])
    all_chromosomes = sorted(gene_annot_df.CHR.unique(), key=stock.natural_key)
    for key in [
        "project_out_covariates_chromosomes",
        "feature_selection_chromosomes",
        "training_chromosomes",
    ]:
        if config_dict[key] is None:
            config_dict[key] = all_chromosomes
    assert config_dict["magma_prefix"] is not None, "--magma_prefix is required."
    logging.info("MAGMA scores provided, loading MAGMA.")
    y_vec, covariates, error_cov, y_ids = stock.read_magma(
        config_dict["magma_prefix"],
        config_dict["use_magma_covariates"],
        config_dict["use_magma_error_cov"],
    )
    if error_cov is not None:
        logging.info("Regularizing MAGMA error covariance.")
        error_cov = stock.regularize_error_cov(error_cov, y_vec, y_ids, gene_annot_df)

    project_out_covariates_y_gene_inds = stock.get_gene_indices_to_use(
        y_ids,
        gene_annot_df,
        config_dict["project_out_covariates_chromosomes"],
        config_dict["project_out_covariates_remove_hla"],
    )
    feature_selection_y_gene_inds = stock.get_gene_indices_to_use(
        y_ids,
        gene_annot_df,
        config_dict["feature_selection_chromosomes"],
        config_dict["feature_selection_remove_hla"],
    )
    training_y_gene_inds = stock.get_gene_indices_to_use(
        y_ids,
        gene_annot_df,
        config_dict["training_chromosomes"],
        config_dict["training_remove_hla"],
    )

    if covariates is not None:
        logging.info(
            f"Projecting {covariates.shape[1]} covariates out of target scores."
        )
        y_proj = stock.project_out_covariates(
            y_vec,
            covariates,
            error_cov,
            y_ids,
            gene_annot_df,
            project_out_covariates_y_gene_inds,
        )
    else:
        y_proj = y_vec

    # --------------------------- Feature selection (stock front-end) --------------------------- #
    logging.info("Computing marginal association table.")
    marginal_assoc_df = stock.compute_marginal_assoc(
        config_dict["feature_mat_prefix"],
        config_dict["num_feature_chunks"],
        y_proj,
        y_ids,
        None,
        error_cov,
        gene_annot_df,
        feature_selection_y_gene_inds,
    )
    selected_features = select_features(config_dict, marginal_assoc_df, stock)
    logging.info(f"{len(selected_features)} features selected.")

    # --------------------------- Training (streaming kernel ridge) --------------------------- #
    training_genes, y_train, intercept, linv, labels = build_training_transforms(
        y_proj, y_ids, error_cov, gene_annot_df, training_y_gene_inds, stock
    )
    logging.info(
        "Accumulating {0}x{0} training Gram matrix over feature chunks.".format(
            len(training_genes)
        )
    )
    gram, rows, x_train_inds = stream_build_gram(
        config_dict["feature_mat_prefix"],
        config_dict["num_feature_chunks"],
        selected_features,
        training_genes,
        intercept,
        linv,
        labels,
        stock,
    )
    # The dense error covariance and Cholesky whitener are no longer needed.
    del error_cov
    logging.info(
        "Selecting ridge penalty via generalized cross-validation on the Gram matrix."
    )
    alpha, best_score, dual_c = ridge_gcv_from_gram(gram, y_train, RIDGE_ALPHAS)
    del gram
    logging.info(f"Selected alpha = {alpha}, GCV score = {best_score}.")
    logging.info("Computing coefficients and PoPS scores in a second streaming pass.")
    cols, betas, preds = stream_coefs_and_preds(
        config_dict["feature_mat_prefix"],
        config_dict["num_feature_chunks"],
        selected_features,
        rows,
        x_train_inds,
        intercept,
        linv,
        labels,
        dual_c,
        stock,
    )

    coefs_df = build_coefs_df(cols, betas, alpha, best_score)
    preds_df = pd.DataFrame([rows, preds]).T
    preds_df.columns = ["ENSGID", "PoPS_Score"]

    # --------------------------- Annotate predictions (stock main) --------------------------- #
    preds_df = preds_df.merge(
        pd.DataFrame(np.array([y_ids, y_vec]).T, columns=["ENSGID", "Y"]),
        how="left",
        on="ENSGID",
    )
    if covariates is not None:
        preds_df = preds_df.merge(
            pd.DataFrame(np.array([y_ids, y_proj]).T, columns=["ENSGID", "Y_proj"]),
            how="left",
            on="ENSGID",
        )
        preds_df["project_out_covariates_gene"] = preds_df.ENSGID.isin(
            y_ids[project_out_covariates_y_gene_inds]
        )
    preds_df["feature_selection_gene"] = preds_df.ENSGID.isin(
        y_ids[feature_selection_y_gene_inds]
    )
    preds_df["training_gene"] = preds_df.ENSGID.isin(y_ids[training_y_gene_inds])

    # --------------------------- Save --------------------------- #
    logging.info("Writing output files.")
    preds_df.to_csv(config_dict["out_prefix"] + ".preds", sep="\t", index=False)
    coefs_df.to_csv(config_dict["out_prefix"] + ".coefs", sep="\t")
    marginal_assoc_df.to_csv(config_dict["out_prefix"] + ".marginals", sep="\t")


def import_stock_pops(source_dir: Path):
    """Import the pinned stock pops.py from the extracted POPs source directory.

    Loaded dynamically (not a static import) because the source directory is only
    known at runtime; this reuses the validated stock front-end functions verbatim.
    """
    sys.path.insert(0, str(source_dir))
    return importlib.import_module("pops")


if __name__ == "__main__":
    source_dir, rest_argv = extract_source_dir(sys.argv[1:])
    stock_pops = import_stock_pops(source_dir)
    config = vars(stock_pops.get_pops_args(rest_argv))
    main_lowmem(config, stock_pops)

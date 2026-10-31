"""
Numerical verification of the claims in coloc_abf_derivation.typ.

Every check is against something independent of the closed form it tests: brute-force
enumeration of all (Q + 1)^2 configurations with full multivariate normal densities, quadrature
of the exact logistic likelihood, or the implementation on branch large-scale-protein-coloc.

Run from the repo root (the implementation is read from the branch with git show):

    pixi r python experiments/claude/coloc_abf_derivation/verify_coloc_abf_derivation.py \
        2>&1 | tee experiments/claude/logs/coloc_abf_derivation_verify.log
"""

import importlib.util
import itertools
import math
import subprocess
import sys
import tempfile
from pathlib import Path
from types import ModuleType

import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal, norm

BRANCH = "large-scale-protein-coloc"
BRANCH_FILE = "mecfs_bio/build_system/task/protein_coloc/coloc_abf.py"
P1, P2, P12 = 1e-4, 1e-4, 1e-5
OMEGA_QUANT = 0.15
OMEGA_CC = 0.2


def load_branch_module() -> ModuleType:
    """Load coloc_abf.py from the branch, which is not checked out on main."""
    source = subprocess.run(
        ["git", "show", f"{BRANCH}:{BRANCH_FILE}"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    path = Path(tempfile.mkdtemp()) / "coloc_abf_branch.py"
    path.write_text(source)
    spec = importlib.util.spec_from_file_location("coloc_abf_branch", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["coloc_abf_branch"] = module
    spec.loader.exec_module(module)
    return module


def log_abf(beta_hat: np.ndarray, v: np.ndarray, omega: float) -> np.ndarray:
    r = omega**2 / (omega**2 + v)
    return 0.5 * (np.log1p(-r) + r * beta_hat**2 / v)


def correlated_genotypes(
    n: int, q: int, rho: float, rng: np.random.Generator
) -> np.ndarray:
    """Additive dosages with AR(1) LD, from two thresholded correlated haplotypes."""

    def haplotype() -> np.ndarray:
        latent = np.zeros((n, q))
        latent[:, 0] = rng.standard_normal(n)
        for j in range(1, q):
            latent[:, j] = rho * latent[:, j - 1] + math.sqrt(
                1 - rho**2
            ) * rng.standard_normal(n)
        thresholds = norm.ppf(rng.uniform(0.1, 0.5, q))
        return (latent < thresholds).astype(float)

    return haplotype() + haplotype()


def null_space_basis(c: np.ndarray) -> np.ndarray:
    """Rows are an orthonormal basis of the orthogonal complement of the columns of c."""
    full, _ = np.linalg.qr(c, mode="complete")
    return full[:, c.shape[1] :].T


def check_theorem_exact_linear_bf() -> None:
    """Exact BF equals the ABF for the Gaussian model, under LD and covariates."""
    rng = np.random.default_rng(1)
    n, q, sigma = 80, 15, 1.3
    g = correlated_genotypes(n, q, rho=0.9, rng=rng)
    c = np.column_stack([np.ones(n), rng.standard_normal((n, 2))])
    y = c @ np.array([0.5, -1.0, 0.3]) + 0.7 * g[:, 4] + sigma * rng.standard_normal(n)
    u_basis = null_space_basis(c)
    u = u_basis @ y
    h = u_basis @ g
    keep = np.sum(h**2, axis=0) > 1e-8
    h = h[:, keep]
    worst = 0.0
    null_log_density = multivariate_normal(
        mean=np.zeros(len(u)), cov=sigma**2 * np.eye(len(u))
    ).logpdf(u)
    for i in range(h.shape[1]):
        hi = h[:, i]
        cov = sigma**2 * np.eye(len(u)) + OMEGA_QUANT**2 * np.outer(hi, hi)
        exact = (
            multivariate_normal(mean=np.zeros(len(u)), cov=cov).logpdf(u)
            - null_log_density
        )
        ols_beta = float(hi @ u / (hi @ hi))
        v = sigma**2 / float(hi @ hi)
        formula = float(log_abf(np.array(ols_beta), np.array(v), OMEGA_QUANT))
        worst = max(worst, abs(exact - formula))
    assert worst < 1e-8, worst
    print(f"[Thm exact linear BF] max |exact - ABF| over variants = {worst:.2e}")


def check_flat_nuisance_equals_projection() -> None:
    """A vague prior on the covariate effects gives the projected-data BF in the limit."""
    rng = np.random.default_rng(2)
    n, sigma = 40, 1.0
    g = rng.binomial(2, 0.3, n).astype(float)
    c = np.column_stack([np.ones(n), rng.standard_normal(n)])
    y = c @ np.array([1.0, -0.5]) + 0.3 * g + sigma * rng.standard_normal(n)
    tau = 1e3

    def log_density(cov: np.ndarray) -> float:
        _, log_det = np.linalg.slogdet(cov)
        return float(-0.5 * (log_det + y @ np.linalg.solve(cov, y)))

    base = sigma**2 * np.eye(n) + tau**2 * c @ c.T
    alt = base + OMEGA_QUANT**2 * np.outer(g, g)
    vague = log_density(alt) - log_density(base)
    basis = null_space_basis(c)
    hi, u = basis @ g, basis @ y
    projected = float(
        log_abf(
            np.array(float(hi @ u / (hi @ hi))),
            np.array(sigma**2 / float(hi @ hi)),
            OMEGA_QUANT,
        )
    )
    assert abs(vague - projected) < 1e-4, (vague, projected)
    print(f"[Rem flat nuisance] |vague - projected| = {abs(vague - projected):.2e}")


def check_prior_lemma() -> None:
    """Enumerate categorical vectors; conditioning gives the claimed configuration odds."""
    q = 6
    p1, p2, p12 = 0.03, 0.05, 0.02
    p0 = 1 - p1 - p2 - p12
    prob = {"0": p0, "1": p1, "2": p2, "b": p12}
    totals = dict.fromkeys(["H0", "H1", "H2", "H3", "H4"], 0.0)
    for categories in itertools.product("012b", repeat=q):
        n_trait1 = sum(c in "1b" for c in categories)
        n_trait2 = sum(c in "2b" for c in categories)
        if n_trait1 > 1 or n_trait2 > 1:
            continue
        weight = math.prod(prob[c] for c in categories)
        n_both = sum(c == "b" for c in categories)
        if n_trait1 == 0 and n_trait2 == 0:
            totals["H0"] += weight
        elif n_trait2 == 0:
            totals["H1"] += weight
        elif n_trait1 == 0:
            totals["H2"] += weight
        elif n_both == 1:
            totals["H4"] += weight
        else:
            totals["H3"] += weight
    rho1, rho2, rho12 = p1 / p0, p2 / p0, p12 / p0
    claimed = {
        "H1": q * rho1,
        "H2": q * rho2,
        "H3": q * (q - 1) * rho1 * rho2,
        "H4": q * rho12,
    }
    for key, value in claimed.items():
        assert math.isclose(totals[key] / totals["H0"], value, rel_tol=1e-12)
    print("[Lem prior] enumerated hypothesis odds match the closed form")


def check_rescaled_prior_is_exact_for_implementation() -> None:
    """Weights (1, p1, p2, p1 p2, p12) are exact for p' = p / (1 + p1 + p2 + p12)."""
    denominator = 1 + P1 + P2 + P12
    p1_, p2_, p12_ = P1 / denominator, P2 / denominator, P12 / denominator
    p0_ = 1 - p1_ - p2_ - p12_
    assert math.isclose(p1_ / p0_, P1, rel_tol=1e-12)
    assert math.isclose(p12_ / p0_, P12, rel_tol=1e-12)
    assert math.isclose((p1_ / p0_) * (p2_ / p0_), P1 * P2, rel_tol=1e-12)
    p0 = 1 - P1 - P2 - P12
    print(
        f"[Prop implemented prior] default p0 = {p0:.5f}; rescale factor "
        f"1/(1+p1+p2+p12) = {1 / denominator:.6f}"
    )


def check_brute_force_posterior(branch: ModuleType) -> None:
    """All (Q + 1)^2 configurations with full densities versus the branch implementation."""
    rng = np.random.default_rng(3)
    q, n1, n2, sigma1, sigma2 = 12, 60, 50, 1.0, 1.0
    g1 = correlated_genotypes(n1, q, rho=0.85, rng=rng)
    g2 = correlated_genotypes(n2, q, rho=0.85, rng=rng)
    y1 = 0.9 * g1[:, 3] + sigma1 * rng.standard_normal(n1)
    y2 = 0.8 * g2[:, 3] + sigma2 * rng.standard_normal(n2)
    u1_basis = null_space_basis(np.ones((n1, 1)))
    u2_basis = null_space_basis(np.ones((n2, 1)))
    u1, h1 = u1_basis @ y1, u1_basis @ g1
    u2, h2 = u2_basis @ y2, u2_basis @ g2

    def log_density(
        u: np.ndarray, h: np.ndarray, index: int, sigma: float, omega: float
    ) -> float:
        cov = sigma**2 * np.eye(len(u))
        if index >= 0:
            cov = cov + omega**2 * np.outer(h[:, index], h[:, index])
        return float(multivariate_normal(mean=np.zeros(len(u)), cov=cov).logpdf(u))

    density1 = [log_density(u1, h1, i, sigma1, OMEGA_QUANT) for i in range(-1, q)]
    density2 = [log_density(u2, h2, i, sigma2, OMEGA_QUANT) for i in range(-1, q)]
    log_weights: dict[str, list[float]] = {f"H{h}": [] for h in range(5)}
    for s1, s2 in itertools.product(range(-1, q), repeat=2):
        log_likelihood = density1[s1 + 1] + density2[s2 + 1]
        if s1 < 0 and s2 < 0:
            log_weights["H0"].append(log_likelihood)
        elif s2 < 0:
            log_weights["H1"].append(log_likelihood + math.log(P1))
        elif s1 < 0:
            log_weights["H2"].append(log_likelihood + math.log(P2))
        elif s1 == s2:
            log_weights["H4"].append(log_likelihood + math.log(P12))
        else:
            log_weights["H3"].append(log_likelihood + math.log(P1) + math.log(P2))
    by_hypothesis = np.array([logsumexp(log_weights[f"H{h}"]) for h in range(5)])
    brute = np.exp(by_hypothesis - logsumexp(by_hypothesis))
    assert len(log_weights["H3"]) == q * (q - 1)

    norms1, norms2 = np.sum(h1**2, axis=0), np.sum(h2**2, axis=0)
    post = branch.ColocAbf().run(
        branch.ColocInput(
            beta1=(h1.T @ u1) / norms1,
            se1=sigma1 / np.sqrt(norms1),
            beta2=(h2.T @ u2) / norms2,
            se2=sigma2 / np.sqrt(norms2),
            maf1=None,
            n1=None,
            maf2=None,
            n2=None,
            model1=branch.QuantitativeTraitModel(sd_y=1.0),
            model2=branch.QuantitativeTraitModel(sd_y=1.0),
        )
    )
    implemented = np.array(
        [post.pp_h0, post.pp_h1, post.pp_h2, post.pp_h3, post.pp_h4]
    )
    assert np.max(np.abs(brute - implemented)) < 1e-10
    print(
        "[Thm posterior] brute force vs branch ColocAbf: max |diff| = "
        f"{np.max(np.abs(brute - implemented)):.2e}; PP = {np.round(brute, 4)}"
    )


def four_number_summary(
    labf1: np.ndarray, labf2: np.ndarray, p1: float, p2: float, p12: float
) -> dict[str, float]:
    """a1, a2, C, kappa and the five posteriors from the four-number form, in log space."""
    l1, l2, l12 = logsumexp(labf1), logsumexp(labf2), logsumexp(labf1 + labf2)
    a1 = float(np.exp(math.log(p1) + l1 - np.logaddexp(0.0, math.log(p1) + l1)))
    a2 = float(np.exp(math.log(p2) + l2 - np.logaddexp(0.0, math.log(p2) + l2)))
    overlap = float(np.exp(l12 - l1 - l2))
    kappa = p12 / (p1 * p2)
    denominator = 1 + (kappa - 1) * a1 * a2 * overlap
    return {
        "a1": a1,
        "a2": a2,
        "C": overlap,
        "kappa": kappa,
        "pp0": (1 - a1) * (1 - a2) / denominator,
        "pp1": a1 * (1 - a2) / denominator,
        "pp2": (1 - a1) * a2 / denominator,
        "pp3": a1 * a2 * (1 - overlap) / denominator,
        "pp4": kappa * a1 * a2 * overlap / denominator,
    }


def check_four_number_form(branch: ModuleType) -> None:
    """The (a1, a2, C, kappa) form equals combine_abf on random log Bayes factors."""
    rng = np.random.default_rng(4)
    worst = 0.0
    priors = branch.ColocPriors()
    for _ in range(2000):
        q = int(rng.integers(2, 300))
        scale = rng.choice([1.0, 5.0, 20.0])
        labf1 = rng.normal(-1.5, scale, q)
        labf2 = rng.normal(-1.5, scale, q)
        post = branch.combine_abf(labf1, labf2, priors)
        reduced = four_number_summary(labf1, labf2, P1, P2, P12)
        for attr, key in [
            ("pp_h0", "pp0"),
            ("pp_h1", "pp1"),
            ("pp_h2", "pp2"),
            ("pp_h3", "pp3"),
            ("pp_h4", "pp4"),
        ]:
            worst = max(worst, abs(getattr(post, attr) - reduced[key]))
        l1, l2, l12 = (
            np.exp(logsumexp(labf1)),
            np.exp(logsumexp(labf2)),
            np.exp(logsumexp(labf1 + labf2)),
        )
        total = (
            1 + P1 * l1 + P2 * l2 + P1 * P2 * (l1 * l2 - l12) + P12 * l12
        )
        factored = (1 + P1 * l1) * (1 + P2 * l2) + (P12 - P1 * P2) * l12
        assert math.isclose(total, factored, rel_tol=1e-9)
    assert worst < 1e-8, worst
    print(f"[Thm four-number form] max |diff| vs combine_abf = {worst:.2e}")


def check_overlap_bounds() -> None:
    """Cauchy-Schwarz bound, uniform baseline, and the H4:H3 thresholds quoted in the text."""
    rng = np.random.default_rng(5)
    for _ in range(1000):
        q = int(rng.integers(2, 200))
        pi1 = rng.dirichlet(np.full(q, 0.2))
        pi2 = rng.dirichlet(np.full(q, 0.2))
        overlap = float(pi1 @ pi2)
        bound = math.sqrt(float(pi1 @ pi1) * float(pi2 @ pi2))
        assert overlap <= bound + 1e-12
    q = 1000
    uniform = np.full(q, 1 / q)
    overlap = float(uniform @ uniform)
    assert math.isclose((q - 1) * overlap / (1 - overlap), 1.0, rel_tol=1e-12)

    def threshold(tau: float, kappa: float) -> float:
        return tau / (tau + kappa * (1 - tau))

    for p12 in (1e-5, 5e-6, 1e-6):
        kappa = p12 / (P1 * P2)
        print(
            f"[Overlap threshold] p12={p12:.0e} kappa={kappa:.0f}: "
            f"C* (tau=0.8) = {threshold(0.8, kappa):.5f}; "
            f"C* (tau=0.5) = {threshold(0.5, kappa):.5f}"
        )


def check_p12_rescaling(branch: ModuleType) -> None:
    """PP4 odds are linear in p12, and rho12* reaches the requested tau."""
    rng = np.random.default_rng(6)
    labf1, labf2 = rng.normal(0, 3, 500), rng.normal(0, 3, 500)
    labf1[10] = labf2[10] = 12.0
    base = branch.combine_abf(labf1, labf2, branch.ColocPriors(p12=1e-5))
    other = branch.combine_abf(labf1, labf2, branch.ColocPriors(p12=1e-6))
    odds_base = base.pp_h4 / (1 - base.pp_h4)
    odds_other = other.pp_h4 / (1 - other.pp_h4)
    assert math.isclose(odds_other, odds_base / 10, rel_tol=1e-9)
    tau = 0.9
    star = 1e-5 * (tau / (1 - tau)) / odds_base
    assert math.isclose(
        branch.combine_abf(labf1, labf2, branch.ColocPriors(p12=star)).pp_h4,
        tau,
        rel_tol=1e-9,
    )
    print("[Cor p12] PP4 odds linear in p12 and threshold formula verified")


def fit_logistic(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """(intercept, slope, Wald variance of the slope) by Newton."""
    design = np.column_stack([np.ones(len(x)), x])
    theta = np.array([math.log(y.mean() / (1 - y.mean())), 0.0])
    for _ in range(50):
        mu = 1 / (1 + np.exp(-design @ theta))
        weights = mu * (1 - mu)
        info = (design * weights[:, None]).T @ design
        step = np.linalg.solve(info, design.T @ (y - mu))
        theta = theta + step
        if np.max(np.abs(step)) < 1e-12:
            break
    mu = 1 / (1 + np.exp(-design @ theta))
    info = (design * (mu * (1 - mu))[:, None]).T @ design
    return float(theta[0]), float(theta[1]), float(np.linalg.inv(info)[1, 1])


def log_integrated_likelihood(
    x: np.ndarray, y: np.ndarray, beta: float, alpha_start: float
) -> float:
    """log of the integral over a flat intercept, by adaptive Gauss-Hermite quadrature."""
    alpha = alpha_start
    for _ in range(50):
        mu = 1 / (1 + np.exp(-(alpha + beta * x)))
        curvature = float(np.sum(mu * (1 - mu)))
        step = float(np.sum(y - mu)) / curvature
        alpha += step
        if abs(step) < 1e-12:
            break
    mu = 1 / (1 + np.exp(-(alpha + beta * x)))
    scale = 1 / math.sqrt(float(np.sum(mu * (1 - mu))))
    nodes, weights = np.polynomial.hermite.hermgauss(24)
    values = []
    for node, weight in zip(nodes, weights, strict=True):
        a = alpha + math.sqrt(2) * scale * node
        eta = a + beta * x
        values.append(
            math.log(weight)
            + node**2
            + float(np.sum(y * eta - np.logaddexp(0.0, eta)))
        )
    return float(logsumexp(values)) + math.log(math.sqrt(2) * scale)


def exact_logistic_log_bf(x: np.ndarray, y: np.ndarray, omega: float) -> float:
    """Exact log BF of the logistic model with a flat intercept, by nested quadrature."""
    intercept, slope, v = fit_logistic(x, y)
    r = omega**2 / (omega**2 + v)
    center, half_width = r * slope, 12 * math.sqrt(v)
    nodes = np.linspace(center - half_width, center + half_width, 241)
    log_values = [
        log_integrated_likelihood(x, y, float(b), intercept) + norm.logpdf(b, 0, omega)
        for b in nodes
    ]
    log_alt = float(logsumexp(log_values)) + math.log(nodes[1] - nodes[0])
    log_null = log_integrated_likelihood(x, y, 0.0, intercept)
    return log_alt - log_null


def check_logistic_approximation() -> None:
    """Exact logistic BF versus the ABF across effect sizes."""
    rng = np.random.default_rng(7)
    n, f = 20000, 0.3
    print("[Prop logistic] N=20000 f=0.3 case fraction~0.3 omega=0.2")
    print("  true_beta   z      logABF    logBF_exact  diff")
    for true_beta in (0.0, 0.05, 0.1, 0.2, 0.4):
        x = rng.binomial(2, f, n).astype(float)
        intercept = math.log(0.3 / 0.7) - true_beta * 2 * f
        y = (rng.uniform(size=n) < 1 / (1 + np.exp(-(intercept + true_beta * x)))).astype(
            float
        )
        _, slope, v = fit_logistic(x, y)
        approx = float(log_abf(np.array(slope), np.array(v), OMEGA_CC))
        exact = exact_logistic_log_bf(x, y, OMEGA_CC)
        print(
            f"  {true_beta:8.2f} {slope / math.sqrt(v):7.2f} {approx:9.3f} "
            f"{exact:11.3f} {approx - exact:8.3f}"
        )


def check_null_variance_case_control() -> None:
    """V = 1 / (2 N f (1 - f) s (1 - s)) matches the Wald variance near the null."""
    rng = np.random.default_rng(8)
    n, f, s = 200000, 0.25, 0.3
    x = rng.binomial(2, f, n).astype(float)
    y = (rng.uniform(size=n) < s).astype(float)
    _, _, v = fit_logistic(x, y)
    formula = 1 / (2 * n * f * (1 - f) * s * (1 - s))
    assert abs(v / formula - 1) < 0.02, (v, formula)
    print(f"[Sec inputs] case-control V ratio Wald/formula = {v / formula:.4f}")


def check_sd_y(branch: ModuleType) -> None:
    """Exact recovery, EAF/MAF symmetry and scale invariance of the quantitative ABF."""
    rng = np.random.default_rng(9)
    f = rng.uniform(0.02, 0.5, 500)
    n = rng.uniform(5e4, 2e5, 500)
    sd_y = 3.7
    se = sd_y / np.sqrt(2 * n * f * (1 - f))
    assert math.isclose(branch.estimate_sd_y(se=se, maf=f, n=n), sd_y, rel_tol=1e-12)
    assert math.isclose(
        branch.estimate_sd_y(se=se, maf=1 - f, n=n), sd_y, rel_tol=1e-12
    )
    beta = rng.normal(0, 0.02, 500)
    scale = 42.0
    a = log_abf(beta, se**2, OMEGA_QUANT * sd_y)
    b = log_abf(scale * beta, (scale * se) ** 2, OMEGA_QUANT * sd_y * scale)
    assert np.max(np.abs(a - b)) < 1e-9
    print("[Sec inputs] sdY exact recovery, f <-> 1-f symmetry, scale invariance OK")


def check_prior_scale_derivative() -> None:
    """d lABF / d log omega = -r + r (1 - r) z^2, against a finite difference."""
    rng = np.random.default_rng(12)
    v = rng.uniform(1e-5, 1e-2, 200)
    z = rng.normal(0, 6, 200)
    beta = z * np.sqrt(v)
    h = 1e-6
    difference = (
        log_abf(beta, v, OMEGA_QUANT * math.exp(h))
        - log_abf(beta, v, OMEGA_QUANT * math.exp(-h))
    ) / (2 * h)
    r = OMEGA_QUANT**2 / (OMEGA_QUANT**2 + v)
    assert np.max(np.abs(difference - (-r + r * (1 - r) * z**2))) < 1e-6
    print("[Sec inputs] d lABF / d log omega formula confirmed")


def check_z_distribution() -> None:
    """z ~ N(lambda R[:, c], R) under the Gaussian model, by direct matrix algebra."""
    rng = np.random.default_rng(10)
    n, q, sigma, beta, c = 200, 10, 1.0, 0.4, 3
    h = correlated_genotypes(n, q, rho=0.8, rng=rng)
    h = h - h.mean(axis=0)
    h_norm = np.sqrt(np.sum(h**2, axis=0))
    transform = (h / h_norm).T / sigma
    mean = transform @ (h[:, c] * beta)
    cov = transform @ (sigma**2 * np.eye(n)) @ transform.T
    correlation = (h / h_norm).T @ (h / h_norm)
    noncentrality = beta * h_norm[c] / sigma
    assert np.max(np.abs(mean - noncentrality * correlation[:, c])) < 1e-10
    assert np.max(np.abs(cov - correlation)) < 1e-10
    print("[Sec LD] z mean = lambda * r and Cov(z) = R confirmed")


def stylized_locus(
    z1: np.ndarray, z2: np.ndarray, se: float = 0.01, omega: float = OMEGA_QUANT
) -> dict[str, float]:
    v = np.full(len(z1), se**2)
    return four_number_summary(
        log_abf(z1 * se, v, omega), log_abf(z2 * se, v, omega), P1, P2, P12
    )


def check_stylized_table() -> None:
    """Noise-free loci for the table in the note: Q = 1000, AR(1) LD, se = 0.01."""
    q = 1000
    rho = 0.98
    positions = np.arange(q)

    def tag(c: int, lam: float) -> np.ndarray:
        return lam * rho ** np.abs(positions - c)

    rows = {
        "shared variant, strong": (tag(500, 10.0), tag(500, 10.0)),
        "distinct variants, 100 apart": (tag(450, 10.0), tag(550, 10.0)),
        "distinct variants, adjacent": (tag(500, 10.0), tag(501, 10.0)),
        "trait 2 null": (tag(500, 10.0), np.zeros(q)),
        "shared variant, trait 2 weak": (tag(500, 10.0), tag(500, 3.5)),
        "two trait-1 signals, weaker shared": (
            tag(300, 10.0) + tag(700, 8.0),
            tag(700, 10.0),
        ),
    }
    print("[Table] Q=1000, r_ij=0.98^|i-j|, se=0.01, omega=0.15, default priors")
    print("  case | a1 a2 C | PP0 PP1 PP2 PP3 PP4")
    for name, (z1, z2) in rows.items():
        s = stylized_locus(z1, z2)
        print(
            f"  {name:36s} | {s['a1']:.3f} {s['a2']:.3f} {s['C']:.2e} | "
            f"{s['pp0']:.3f} {s['pp1']:.3f} {s['pp2']:.3f} {s['pp3']:.3f} {s['pp4']:.3f}"
        )


def check_clique() -> None:
    """Perfect-LD clique: C = 1 / m and H4:H3 odds = kappa / (m - 1)."""
    q, m = 1000, 40
    z = np.zeros(q)
    z[:m] = 12.0
    s = stylized_locus(z, z.copy())
    assert abs(s["C"] - 1 / m) < 1e-6
    ratio = s["pp4"] / s["pp3"]
    assert math.isclose(ratio, s["kappa"] * s["C"] / (1 - s["C"]), rel_tol=1e-9)
    assert abs(ratio - s["kappa"] / (m - 1)) / ratio < 1e-4
    print(f"[Prop clique] m={m}: C={s['C']:.5f}, PP4/PP3={ratio:.3f}")


def check_multi_causal() -> None:
    """H4:H3 odds approximately kappa exp(-(z_a^2 - z_b^2) / 2) with two trait-1 signals."""
    q = 1000
    for z_a, z_b in ((10.0, 8.0), (9.0, 8.0), (8.5, 8.0), (8.0, 8.0)):
        z1 = np.zeros(q)
        z2 = np.zeros(q)
        z1[100], z1[900] = z_a, z_b
        z2[900] = 10.0
        s = stylized_locus(z1, z2)
        predicted = 1000 * math.exp(-(z_a**2 - z_b**2) / 2)
        print(
            f"[Prop two signals] z_a={z_a} z_b={z_b}: PP4/PP3={s['pp4'] / s['pp3']:.4g}, "
            f"kappa*exp(-diff/2)={predicted:.4g}"
        )
    print(f"  2 log kappa = {2 * math.log(1000):.2f}")


def check_consistency() -> None:
    """All four truth patterns converge as N grows (fixed Q, identifiable causal variant)."""
    rng = np.random.default_rng(11)
    q, c, d, rho = 200, 50, 150, 0.9
    correlation = rho ** np.abs(np.arange(q)[:, None] - np.arange(q)[None, :])
    cholesky = np.linalg.cholesky(correlation)
    effect, residual_sd, f = 0.15, 1.0, 0.3

    def z_vector(n: float, causal: int | None) -> np.ndarray:
        noise = cholesky @ rng.standard_normal(q)
        if causal is None:
            return noise
        lam = effect * math.sqrt(n * 2 * f * (1 - f)) / residual_sd
        return lam * correlation[:, causal] + noise

    print("[Prop consistency] N, PP0..PP4 for four truths")
    for label, s1, s2 in (
        ("H0", None, None),
        ("H1", c, None),
        ("H3", c, d),
        ("H4", c, c),
    ):
        for n in (1e3, 1e4, 1e5, 1e6):
            se = residual_sd / math.sqrt(n * 2 * f * (1 - f))
            s = stylized_locus(z_vector(n, s1), z_vector(n, s2), se=se)
            print(
                f"  truth {label} N={n:9.0f}: "
                + " ".join(f"{s[k]:.3f}" for k in ("pp0", "pp1", "pp2", "pp3", "pp4"))
            )


def check_quoted_numbers() -> None:
    """Numbers quoted in passing in the note."""
    print(f"[Numbers] p12 cond prob p12/(p1+p12) = {P12 / (P1 + P12):.4f}")
    print(
        f"[Numbers] variance explained at beta=0.15, f=0.3: "
        f"{0.15**2 * 2 * 0.3 * 0.7:.5f}"
    )
    print(f"[Numbers] 95% OR interval for omega=0.2: {math.exp(-1.96 * 0.2):.3f}, "
          f"{math.exp(1.96 * 0.2):.3f}")
    print(f"[Numbers] 95% interval for omega=0.15: +-{1.96 * 0.15:.3f} sdY")
    for se in (0.05, 0.01):
        v = se**2
        r = OMEGA_QUANT**2 / (OMEGA_QUANT**2 + v)
        print(
            f"[Numbers] se={se}: r={r:.5f}, ABF>1 iff |z|>"
            f"{math.sqrt(-math.log1p(-r) / r):.3f}"
        )
    for rsq in (0.01, 0.1, 0.3):
        print(
            f"[Numbers] R2={rsq}: -log(1-R2)={-math.log(1 - rsq):.4f}, "
            f"R2/(1-R2)={rsq / (1 - rsq):.4f}"
        )
    print(f"[Numbers] Q(p1+p2+p12) for Q=1e4 etc: "
          f"{[(q, round(q * (P1 + P2 + P12) + q * (q - 1) * P1 * P2, 3)) for q in (1000, 4000, 10000)]}")
    for q in (1000, 10000):
        rho1, rho2, rho12 = P1, P2, P12
        z = 1 + q * rho1 + q * rho2 + q * (q - 1) * rho1 * rho2 + q * rho12
        print(
            f"[Numbers] prior P(H0..H4), Q={q}: "
            f"{[round(x / z, 4) for x in (1, q * rho1, q * rho2, q * (q - 1) * rho1 * rho2, q * rho12)]}"
        )


def main() -> None:
    branch = load_branch_module()
    check_theorem_exact_linear_bf()
    check_flat_nuisance_equals_projection()
    check_prior_lemma()
    check_rescaled_prior_is_exact_for_implementation()
    check_brute_force_posterior(branch)
    check_four_number_form(branch)
    check_overlap_bounds()
    check_p12_rescaling(branch)
    check_null_variance_case_control()
    check_sd_y(branch)
    check_prior_scale_derivative()
    check_z_distribution()
    check_clique()
    check_multi_causal()
    check_stylized_table()
    check_consistency()
    check_logistic_approximation()
    check_quoted_numbers()
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()

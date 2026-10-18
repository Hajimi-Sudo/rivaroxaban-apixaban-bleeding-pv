"""PATS-PV: persistent association and transient shock pharmacovigilance.

Large simulation banks use a deterministic penalized MAP plus Laplace
approximation. The final real-data analysis uses PyMC NUTS through
``fit_pats_mcmc`` and must satisfy the diagnostics frozen in the plan.
"""

from __future__ import annotations

import argparse
import json
import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, gammaln, logsumexp


DEVELOPMENT_SEEDS = list(range(110_001, 110_101))
TEST_SEEDS = list(range(910_001, 910_501))
SCENARIOS = ("null", "constant", "gradual_drift", "single_shock", "opposing_shocks")
SHOCK_THRESHOLD = math.log(1.5)


@dataclass(frozen=True)
class PatsConfig:
    engine: str = "v1_l1"
    lambda_alpha: float = 8.0
    lambda_theta: float = 12.0
    lambda_shock: float = 3.0
    shock_epsilon: float = 0.03
    shock_probability_cutoff: float = 0.80
    spike_sd: float = 0.10
    slab_sd: float = 1.00
    shock_prior_probability: float = 0.02
    ridge: float = 1e-5
    posterior_draws: int = 1000


def validate_config(config: PatsConfig) -> None:
    if config.engine not in {"v1_l1", "v2_spike_slab"}:
        raise ValueError("engine must be v1_l1 or v2_spike_slab")
    if min(config.lambda_alpha, config.lambda_theta, config.lambda_shock) <= 0:
        raise ValueError("All PATS penalty weights must be positive")
    if config.shock_epsilon <= 0:
        raise ValueError("shock_epsilon must be positive")
    if not 0 < config.shock_probability_cutoff < 1:
        raise ValueError("shock_probability_cutoff must lie strictly between 0 and 1")
    if not 0 < config.spike_sd < config.slab_sd:
        raise ValueError("Require 0 < spike_sd < slab_sd")
    if not 0 < config.shock_prior_probability < 0.5:
        raise ValueError("shock_prior_probability must lie between 0 and 0.5")
    if config.posterior_draws < 100:
        raise ValueError("posterior_draws must be at least 100")


def validate_counts(frame: pd.DataFrame) -> None:
    required = {
        "period", "rivaroxaban_reports", "rivaroxaban_aub",
        "apixaban_reports", "apixaban_aub",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing count columns: {sorted(missing)}")
    for total, events in (
        ("rivaroxaban_reports", "rivaroxaban_aub"),
        ("apixaban_reports", "apixaban_aub"),
    ):
        if (frame[total] < 0).any() or (frame[events] < 0).any() or (frame[events] > frame[total]).any():
            raise ValueError(f"Invalid binomial cells for {events}/{total}")


def _binomial_nll(y: np.ndarray, n: np.ndarray, eta: np.ndarray) -> float:
    return float(np.sum(n * np.logaddexp(0.0, eta) - y * eta))


def _binomial_log_likelihood(y: np.ndarray, n: np.ndarray, eta: np.ndarray) -> float:
    combinatorial = gammaln(n + 1) - gammaln(y + 1) - gammaln(n - y + 1)
    return float(np.sum(combinatorial + y * eta - n * np.logaddexp(0.0, eta)))


def unpack(params: np.ndarray, t: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    alpha = params[:t]
    theta = params[t : 2 * t]
    shock_raw = params[2 * t :]
    shock = shock_raw - shock_raw.mean()
    return alpha, theta, shock


def pats_objective(params: np.ndarray, y_r: np.ndarray, n_r: np.ndarray, y_a: np.ndarray, n_a: np.ndarray, config: PatsConfig) -> float:
    t = len(y_r)
    alpha, theta, shock = unpack(params, t)
    delta = theta + shock
    eta_r = alpha + delta / 2
    eta_a = alpha - delta / 2
    nll = _binomial_nll(y_r, n_r, eta_r) + _binomial_nll(y_a, n_a, eta_a)
    penalty = 0.5 * config.lambda_alpha * np.sum(np.diff(alpha) ** 2)
    penalty += 0.5 * config.lambda_theta * np.sum(np.diff(theta) ** 2)
    penalty += config.lambda_shock * np.sum(np.sqrt(shock**2 + config.shock_epsilon**2))
    penalty += 0.5 * (alpha[0] / 3.0) ** 2 + 0.5 * (theta[0] / 2.0) ** 2
    penalty += 0.5 * config.ridge * np.sum(params**2)
    return float(nll + penalty)


def _normal_logpdf(x: np.ndarray, sd: float) -> np.ndarray:
    return -0.5 * (x / sd) ** 2 - math.log(sd) - 0.5 * math.log(2 * math.pi)


def shock_slab_probability(shock: np.ndarray, config: PatsConfig) -> np.ndarray:
    """Posterior mixture responsibility of the transient slab component."""
    log_spike = math.log1p(-config.shock_prior_probability) + _normal_logpdf(shock, config.spike_sd)
    log_slab = math.log(config.shock_prior_probability) + _normal_logpdf(shock, config.slab_sd)
    denominator = logsumexp(np.vstack([log_spike, log_slab]), axis=0)
    return np.exp(log_slab - denominator)


def pats_v2_objective(
    params: np.ndarray,
    y_r: np.ndarray,
    n_r: np.ndarray,
    y_a: np.ndarray,
    n_a: np.ndarray,
    config: PatsConfig,
) -> float:
    """Spike-and-slab transient-state objective for PATS-PV v2."""
    t = len(y_r)
    alpha, theta, shock = unpack(params, t)
    delta = theta + shock
    eta_r = alpha + delta / 2
    eta_a = alpha - delta / 2
    nll = _binomial_nll(y_r, n_r, eta_r) + _binomial_nll(y_a, n_a, eta_a)
    penalty = 0.5 * config.lambda_alpha * np.sum(np.diff(alpha) ** 2)
    penalty += 0.5 * config.lambda_theta * np.sum(np.diff(theta) ** 2)
    log_spike = math.log1p(-config.shock_prior_probability) + _normal_logpdf(shock, config.spike_sd)
    log_slab = math.log(config.shock_prior_probability) + _normal_logpdf(shock, config.slab_sd)
    penalty -= np.sum(logsumexp(np.vstack([log_spike, log_slab]), axis=0))
    penalty += 0.5 * (alpha[0] / 3.0) ** 2 + 0.5 * (theta[0] / 2.0) ** 2
    penalty += 0.5 * config.ridge * np.sum(params**2)
    return float(nll + penalty)


def initial_values(y_r: np.ndarray, n_r: np.ndarray, y_a: np.ndarray, n_a: np.ndarray) -> np.ndarray:
    smooth = 0.5
    p_r = (y_r + smooth) / (n_r + 2 * smooth)
    p_a = (y_a + smooth) / (n_a + 2 * smooth)
    logit_r = np.log(p_r / (1 - p_r))
    logit_a = np.log(p_a / (1 - p_a))
    alpha = (logit_r + logit_a) / 2
    theta = logit_r - logit_a
    shock = np.zeros_like(theta)
    return np.concatenate([alpha, theta, shock])


def spike_slab_initial_values(
    y_r: np.ndarray, n_r: np.ndarray, y_a: np.ndarray, n_a: np.ndarray
) -> list[np.ndarray]:
    """Use both no-shock and robust local-residual starts for the non-convex v2 fit."""
    base = initial_values(y_r, n_r, y_a, n_a)
    t = len(y_r)
    raw_theta = base[t : 2 * t]
    padded = np.pad(raw_theta, (2, 2), mode="edge")
    smooth_theta = np.convolve(padded, np.ones(5) / 5, mode="valid")
    shock = raw_theta - smooth_theta
    shock -= shock.mean()
    residual_start = np.concatenate([base[:t], smooth_theta, shock])
    return [base, residual_start]


def fit_pats_map(frame: pd.DataFrame, config: PatsConfig = PatsConfig()) -> dict[str, object]:
    validate_config(config)
    validate_counts(frame)
    y_r = frame["rivaroxaban_aub"].to_numpy(dtype=float)
    n_r = frame["rivaroxaban_reports"].to_numpy(dtype=float)
    y_a = frame["apixaban_aub"].to_numpy(dtype=float)
    n_a = frame["apixaban_reports"].to_numpy(dtype=float)
    objective = pats_v2_objective if config.engine == "v2_spike_slab" else pats_objective
    starts = (
        spike_slab_initial_values(y_r, n_r, y_a, n_a)
        if config.engine == "v2_spike_slab"
        else [initial_values(y_r, n_r, y_a, n_a)]
    )
    candidates = [
        minimize(
            objective,
            start,
            args=(y_r, n_r, y_a, n_a, config),
            method="L-BFGS-B",
            options={"maxiter": 2500, "ftol": 1e-10, "gtol": 1e-6},
        )
        for start in starts
    ]
    successful = [candidate for candidate in candidates if candidate.success]
    result = min(successful or candidates, key=lambda candidate: candidate.fun)
    if not result.success:
        raise RuntimeError(f"PATS MAP optimization failed: {result.message}")
    alpha, theta, shock = unpack(result.x, len(frame))

    # L-BFGS supplies a positive-definite inverse-Hessian approximation.
    covariance = np.asarray(result.hess_inv.todense(), dtype=float)
    covariance = (covariance + covariance.T) / 2
    rng = np.random.default_rng(73_021)
    try:
        draws = rng.multivariate_normal(result.x, covariance, size=config.posterior_draws, check_valid="raise")
    except (ValueError, np.linalg.LinAlgError):
        covariance += np.eye(len(result.x)) * 1e-5
        draws = rng.multivariate_normal(result.x, covariance, size=config.posterior_draws, check_valid="ignore")
    persistent_draws = np.median(draws[:, len(frame) : 2 * len(frame)], axis=1)
    persistent = float(np.median(theta))
    if config.engine == "v2_spike_slab":
        shock_probability = shock_slab_probability(shock, config)
    else:
        shock_probability = np.mean(
            np.abs(draws[:, 2 * len(frame) :] - draws[:, 2 * len(frame) :].mean(axis=1, keepdims=True)) > SHOCK_THRESHOLD,
            axis=0,
        )
    return {
        "success": True,
        "objective": float(result.fun),
        "iterations": int(result.nit),
        "log_predictive_density": float(
            _binomial_log_likelihood(y_r, n_r, alpha + (theta + shock) / 2)
            + _binomial_log_likelihood(y_a, n_a, alpha - (theta + shock) / 2)
        ),
        "persistent_log_or": persistent,
        "persistent_or": float(math.exp(persistent)),
        "persistent_ci95_lower": float(math.exp(np.quantile(persistent_draws, 0.025))),
        "persistent_ci95_upper": float(math.exp(np.quantile(persistent_draws, 0.975))),
        "persistent_positive_probability": float(np.mean(persistent_draws > 0)),
        "theta": theta.tolist(),
        "shock": shock.tolist(),
        "shock_probability": shock_probability.tolist(),
        "periods": frame["period"].astype(str).tolist(),
        "inference": (
            "spike_slab_MAP_with_Laplace_persistent_interval"
            if config.engine == "v2_spike_slab"
            else "penalized_MAP_with_Laplace_approximation"
        ),
        "config": config.__dict__,
    }


def smooth_objective(params: np.ndarray, y_r: np.ndarray, n_r: np.ndarray, y_a: np.ndarray, n_a: np.ndarray, config: PatsConfig) -> float:
    t = len(y_r)
    alpha, theta = params[:t], params[t:]
    nll = _binomial_nll(y_r, n_r, alpha + theta / 2) + _binomial_nll(y_a, n_a, alpha - theta / 2)
    penalty = 0.5 * config.lambda_alpha * np.sum(np.diff(alpha) ** 2)
    penalty += 0.5 * config.lambda_theta * np.sum(np.diff(theta) ** 2)
    penalty += 0.5 * (alpha[0] / 3.0) ** 2 + 0.5 * (theta[0] / 2.0) ** 2
    penalty += 0.5 * config.ridge * np.sum(params**2)
    return float(nll + penalty)


def fit_smooth_map(frame: pd.DataFrame, config: PatsConfig = PatsConfig()) -> dict[str, float | bool]:
    """No-shock smooth random-walk baseline fit to the identical counts."""
    validate_counts(frame)
    y_r = frame["rivaroxaban_aub"].to_numpy(dtype=float)
    n_r = frame["rivaroxaban_reports"].to_numpy(dtype=float)
    y_a = frame["apixaban_aub"].to_numpy(dtype=float)
    n_a = frame["apixaban_reports"].to_numpy(dtype=float)
    start = initial_values(y_r, n_r, y_a, n_a)[: 2 * len(frame)]
    result = minimize(
        smooth_objective, start, args=(y_r, n_r, y_a, n_a, config), method="L-BFGS-B",
        options={"maxiter": 2500, "ftol": 1e-10, "gtol": 1e-6},
    )
    if not result.success:
        raise RuntimeError(f"Smooth baseline optimization failed: {result.message}")
    theta = result.x[len(frame) :]
    alpha = result.x[: len(frame)]
    return {
        "success": True,
        "persistent_log_or": float(np.median(theta)),
        "persistent_or": float(math.exp(np.median(theta))),
        "objective": float(result.fun),
        "log_predictive_density": float(
            _binomial_log_likelihood(y_r, n_r, alpha + theta / 2)
            + _binomial_log_likelihood(y_a, n_a, alpha - theta / 2)
        ),
    }


def generate_scenario(seed: int, scenario: str, template: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    rng = np.random.default_rng(seed)
    t = len(template)
    n_r = np.maximum(template["rivaroxaban_reports"].to_numpy(dtype=int), 8)
    n_a = np.maximum(template["apixaban_reports"].to_numpy(dtype=int), 8)
    x = np.linspace(-1, 1, t)
    alpha = -3.0 + 0.15 * np.sin(np.linspace(0, 3 * np.pi, t))
    persistent = np.zeros(t)
    shock = np.zeros(t)
    if scenario == "constant":
        persistent[:] = math.log(2.0)
    elif scenario == "gradual_drift":
        persistent = math.log(1.2) + (x + 1) / 2 * (math.log(2.5) - math.log(1.2))
    elif scenario == "single_shock":
        persistent[:] = math.log(1.5)
        shock[t // 2] = math.log(4.0)
    elif scenario == "opposing_shocks":
        persistent[:] = math.log(1.5)
        shock[t // 3] = math.log(4.0)
        shock[(2 * t) // 3] = -math.log(3.0)
    elif scenario != "null":
        raise ValueError(scenario)
    delta = persistent + shock
    p_r = expit(alpha + delta / 2)
    p_a = expit(alpha - delta / 2)
    simulated = pd.DataFrame(
        {
            "period": template["period"].astype(str),
            "rivaroxaban_reports": n_r,
            "rivaroxaban_aub": rng.binomial(n_r, p_r),
            "apixaban_reports": n_a,
            "apixaban_aub": rng.binomial(n_a, p_a),
        }
    )
    truth = {
        "persistent_log_or": float(np.median(persistent)),
        "shock_indices": np.flatnonzero(np.abs(shock) > 0).tolist(),
    }
    return simulated, truth


def _simulation_job(args: tuple[int, str, pd.DataFrame, PatsConfig]) -> dict[str, object]:
    seed, scenario, template, config = args
    np.random.seed(seed)
    frame, truth = generate_scenario(seed, scenario, template)
    try:
        fit = fit_pats_map(frame, config)
        smooth = fit_smooth_map(frame, config)
        detected = {
            i for i, p in enumerate(fit["shock_probability"])
            if p >= config.shock_probability_cutoff
        }
        actual = set(truth["shock_indices"])
        tp, fp, fn = len(detected & actual), len(detected - actual), len(actual - detected)
        precision = tp / (tp + fp) if tp + fp else (1.0 if not actual else 0.0)
        recall = tp / (tp + fn) if tp + fn else 1.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        ci_contains = fit["persistent_ci95_lower"] <= math.exp(truth["persistent_log_or"]) <= fit["persistent_ci95_upper"]
        positive = fit["persistent_ci95_lower"] > 1.0
        return {
            "seed": seed,
            "scenario": scenario,
            "success": True,
            "truth_log_or": truth["persistent_log_or"],
            "estimate_log_or": fit["persistent_log_or"],
            "squared_error": (fit["persistent_log_or"] - truth["persistent_log_or"]) ** 2,
            "smooth_squared_error": (smooth["persistent_log_or"] - truth["persistent_log_or"]) ** 2,
            "log_predictive_density": fit["log_predictive_density"],
            "smooth_log_predictive_density": smooth["log_predictive_density"],
            "coverage": bool(ci_contains),
            "positive_signal": bool(positive),
            "shock_precision": precision,
            "shock_recall": recall,
            "shock_f1": f1,
        }
    except Exception as exc:  # captured in the complete audit table
        return {"seed": seed, "scenario": scenario, "success": False, "error": repr(exc)}


def run_simulation_bank(
    template: pd.DataFrame, bank: str, output_dir: Path, workers: int, config: PatsConfig
) -> dict[str, object]:
    validate_config(config)
    seeds = DEVELOPMENT_SEEDS if bank == "development" else TEST_SEEDS
    jobs = [(seed, scenario, template, config) for scenario in SCENARIOS for seed in seeds]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(_simulation_job, jobs, chunksize=5))
    output_dir.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(results)
    table.to_csv(output_dir / f"{bank}_simulation_results.csv", index=False)
    completed = table[table["success"].eq(True)].copy()
    summary = (
        completed.groupby("scenario")
        .agg(
            n=("seed", "size"),
            coverage=("coverage", "mean"),
            rmse=("squared_error", lambda x: float(np.sqrt(np.mean(x)))),
            smooth_rmse=("smooth_squared_error", lambda x: float(np.sqrt(np.mean(x)))),
            positive_rate=("positive_signal", "mean"),
            shock_precision=("shock_precision", "mean"),
            shock_recall=("shock_recall", "mean"),
            shock_f1=("shock_f1", "mean"),
            mean_log_predictive_density=("log_predictive_density", "mean"),
            mean_smooth_log_predictive_density=("smooth_log_predictive_density", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(output_dir / f"{bank}_simulation_summary.csv", index=False)
    (output_dir / f"{bank}_config.json").write_text(
        json.dumps(config.__dict__, indent=2), encoding="utf-8"
    )
    gate = {
        "all_jobs_completed": bool(len(completed) == len(jobs)),
        "coverage_deviation_lte_0_05": bool((summary["coverage"] - 0.95).abs().max() <= 0.05),
        "null_type1_lte_0_07": bool(summary.loc[summary["scenario"].eq("null"), "positive_rate"].iloc[0] <= 0.07),
        "shock_f1_gte_0_70": bool(
            summary.loc[summary["scenario"].isin(["single_shock", "opposing_shocks"]), "shock_f1"].min() >= 0.70
        ),
        "persistent_rmse_within_5pct_of_smooth_baseline": bool((summary["rmse"] <= 1.05 * summary["smooth_rmse"]).all()),
    }
    gate["pass"] = bool(all(gate.values()))
    (output_dir / f"{bank}_gate.json").write_text(json.dumps(gate, indent=2), encoding="utf-8")
    return {"bank": bank, "jobs": len(jobs), "gate": gate}


def fit_pats_mcmc(
    frame: pd.DataFrame,
    output_dir: Path,
    draws: int = 1500,
    tune: int = 1500,
    chains: int = 5,
    config: PatsConfig = PatsConfig(),
) -> dict[str, object]:
    """Fit the final PATS state-space model using five NUTS chains."""
    validate_config(config)
    validate_counts(frame)
    try:
        import arviz as az
        import pymc as pm
        import pytensor.tensor as pt
    except ImportError as exc:
        raise RuntimeError("PyMC, PyTensor and ArviZ are required for final PATS-PV inference") from exc

    y_r = frame["rivaroxaban_aub"].to_numpy(dtype=int)
    n_r = frame["rivaroxaban_reports"].to_numpy(dtype=int)
    y_a = frame["apixaban_aub"].to_numpy(dtype=int)
    n_a = frame["apixaban_reports"].to_numpy(dtype=int)
    t = len(frame)
    coords = {"quarter": frame["period"].astype(str).tolist()}
    with pm.Model(coords=coords) as model:
        alpha0 = pm.Normal("alpha0", mu=-3, sigma=2)
        sigma_alpha = pm.HalfNormal("sigma_alpha", sigma=0.35)
        alpha_steps = pm.Normal("alpha_steps", 0, 1, dims="quarter")
        alpha_walk = pt.cumsum(alpha_steps) * sigma_alpha
        alpha = pm.Deterministic("alpha", alpha0 + alpha_walk - pt.mean(alpha_walk), dims="quarter")

        mu = pm.Normal("mu", 0, 1.5)
        sigma_theta = pm.HalfNormal("sigma_theta", sigma=0.25)
        theta_steps = pm.Normal("theta_steps", 0, 1, dims="quarter")
        theta_walk = pt.cumsum(theta_steps) * sigma_theta
        theta = pm.Deterministic("theta", mu + theta_walk - pt.mean(theta_walk), dims="quarter")

        if config.engine == "v2_spike_slab":
            spike = pm.Normal.dist(mu=0, sigma=config.spike_sd, shape=t)
            slab = pm.Normal.dist(mu=0, sigma=config.slab_sd, shape=t)
            z_unscaled = pm.Mixture(
                "shock_mixture",
                w=np.array([1 - config.shock_prior_probability, config.shock_prior_probability]),
                comp_dists=[spike, slab],
                dims="quarter",
            )
        else:
            tau = pm.HalfStudentT("tau", nu=3, sigma=0.25)
            local = pm.HalfStudentT("local", nu=3, sigma=1, dims="quarter")
            slab2 = pm.InverseGamma("slab2", alpha=2, beta=2)
            local_tilde = pt.sqrt(slab2 * local**2 / (slab2 + tau**2 * local**2))
            z_raw = pm.Normal("z_raw", 0, 1, dims="quarter")
            z_unscaled = z_raw * tau * local_tilde
        shock = pm.Deterministic("shock", z_unscaled - pt.mean(z_unscaled), dims="quarter")
        delta = theta + shock

        pm.Binomial("y_r", n=n_r, logit_p=alpha + delta / 2, observed=y_r, dims="quarter")
        pm.Binomial("y_a", n=n_a, logit_p=alpha - delta / 2, observed=y_a, dims="quarter")
        # PyMC 6 renamed the prior-predictive sample-count argument from
        # ``samples`` to ``draws``.  The project environment is pinned to
        # PyMC >= 6, so use the current public API explicitly.
        prior = pm.sample_prior_predictive(draws=500, random_seed=7_301)
        trace = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            cores=min(chains, 5),
            random_seed=[11, 23, 37, 53, 71][:chains],
            target_accept=0.95,
            return_inferencedata=True,
        )
        posterior_predictive = pm.sample_posterior_predictive(
            trace, var_names=["y_r", "y_a"], random_seed=81_113, return_inferencedata=True
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    trace.to_netcdf(output_dir / "pats_pv_trace.nc")
    prior.to_netcdf(output_dir / "pats_prior_predictive.nc")
    posterior_predictive.to_netcdf(output_dir / "pats_posterior_predictive.nc")
    diagnostic_vars = ["mu", "sigma_alpha", "sigma_theta", "theta", "shock"]
    if config.engine == "v1_l1":
        diagnostic_vars.append("tau")
    diagnostics = az.summary(trace, var_names=diagnostic_vars, round_to=None)
    diagnostics.to_csv(output_dir / "pats_pv_diagnostics.csv")
    divergences = int(trace.sample_stats["diverging"].sum().item())
    mu_draws = trace.posterior["mu"].values.reshape(-1)
    shock_draws = trace.posterior["shock"].values.reshape(-1, t)
    ppc_r = posterior_predictive.posterior_predictive["y_r"].values.reshape(-1, t)
    ppc_a = posterior_predictive.posterior_predictive["y_a"].values.reshape(-1, t)
    ppc_r_cover = np.mean((y_r >= np.quantile(ppc_r, 0.025, axis=0)) & (y_r <= np.quantile(ppc_r, 0.975, axis=0)))
    ppc_a_cover = np.mean((y_a >= np.quantile(ppc_a, 0.025, axis=0)) & (y_a <= np.quantile(ppc_a, 0.975, axis=0)))
    if config.engine == "v2_spike_slab":
        responsibility_draws = np.vstack(
            [shock_slab_probability(row, config) for row in shock_draws]
        )
        shock_probability = responsibility_draws.mean(axis=0)
    else:
        shock_probability = np.mean(np.abs(shock_draws) > SHOCK_THRESHOLD, axis=0)
    result = {
        "persistent_or_median": float(np.exp(np.median(mu_draws))),
        "persistent_or_ci95_lower": float(np.exp(np.quantile(mu_draws, 0.025))),
        "persistent_or_ci95_upper": float(np.exp(np.quantile(mu_draws, 0.975))),
        "persistent_positive_probability": float(np.mean(mu_draws > 0)),
        "shock_probability": shock_probability.tolist(),
        "config": config.__dict__,
        "periods": frame["period"].astype(str).tolist(),
        "max_rhat": float(diagnostics["r_hat"].max()),
        "min_ess_bulk": float(diagnostics["ess_bulk"].min()),
        "min_ess_tail": float(diagnostics["ess_tail"].min()),
        "divergences": divergences,
        "posterior_predictive_quarter_coverage_rivaroxaban": float(ppc_r_cover),
        "posterior_predictive_quarter_coverage_apixaban": float(ppc_a_cover),
    }
    result["diagnostic_gate_pass"] = bool(
        result["max_rhat"] <= 1.01
        and result["min_ess_bulk"] >= 400
        and result["min_ess_tail"] >= 400
        and divergences == 0
    )
    (output_dir / "pats_pv_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def run_approximation_check(
    template: pd.DataFrame,
    output_dir: Path,
    draws: int = 500,
    tune: int = 500,
    config: PatsConfig = PatsConfig(),
) -> dict[str, object]:
    """Bridge the scalable Laplace simulation engine to the final horseshoe MCMC."""
    rows = []
    check_seeds = [310_001, 310_002]
    for scenario in SCENARIOS:
        for seed in check_seeds:
            frame, truth = generate_scenario(seed, scenario, template)
            bridge_config = PatsConfig(**{**config.__dict__, "posterior_draws": 500})
            map_fit = fit_pats_map(frame, bridge_config)
            run_dir = output_dir / f"{scenario}_{seed}"
            mcmc = fit_pats_mcmc(
                frame, run_dir, draws=draws, tune=tune, chains=2, config=config
            )
            map_log = math.log(map_fit["persistent_or"])
            mcmc_log = math.log(mcmc["persistent_or_median"])
            interval_overlap = not (
                map_fit["persistent_ci95_upper"] < mcmc["persistent_or_ci95_lower"]
                or mcmc["persistent_or_ci95_upper"] < map_fit["persistent_ci95_lower"]
            )
            rows.append(
                {
                    "scenario": scenario,
                    "seed": seed,
                    "truth_log_or": truth["persistent_log_or"],
                    "map_log_or": map_log,
                    "mcmc_log_or": mcmc_log,
                    "absolute_log_or_difference": abs(map_log - mcmc_log),
                    "interval_overlap": interval_overlap,
                    "mcmc_diagnostic_gate": mcmc["diagnostic_gate_pass"],
                }
            )
    output_dir.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(rows)
    table.to_csv(output_dir / "approximation_check.csv", index=False)
    result = {
        "n": len(table),
        "max_absolute_log_or_difference": float(table["absolute_log_or_difference"].max()),
        "all_intervals_overlap": bool(table["interval_overlap"].all()),
        "all_mcmc_diagnostics_pass": bool(table["mcmc_diagnostic_gate"].all()),
    }
    result["pass"] = bool(
        result["max_absolute_log_or_difference"] <= 0.25
        and result["all_intervals_overlap"]
        and result["all_mcmc_diagnostics_pass"]
    )
    (output_dir / "approximation_gate.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mode", choices=("map", "development", "test", "mcmc", "approximation-check"), required=True)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--draws", type=int, default=1500)
    parser.add_argument("--tune", type=int, default=1500)
    parser.add_argument("--engine", choices=("v1_l1", "v2_spike_slab"), default="v1_l1")
    parser.add_argument("--lambda-alpha", type=float, default=8.0)
    parser.add_argument("--lambda-theta", type=float, default=12.0)
    parser.add_argument("--lambda-shock", type=float, default=3.0)
    parser.add_argument("--shock-epsilon", type=float, default=0.03)
    parser.add_argument("--shock-probability-cutoff", type=float, default=0.80)
    parser.add_argument("--spike-sd", type=float, default=0.10)
    parser.add_argument("--slab-sd", type=float, default=1.00)
    parser.add_argument("--shock-prior-probability", type=float, default=0.02)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    counts = pd.read_csv(args.counts)
    config = PatsConfig(
        engine=args.engine,
        lambda_alpha=args.lambda_alpha,
        lambda_theta=args.lambda_theta,
        lambda_shock=args.lambda_shock,
        shock_epsilon=args.shock_epsilon,
        shock_probability_cutoff=args.shock_probability_cutoff,
        spike_sd=args.spike_sd,
        slab_sd=args.slab_sd,
        shock_prior_probability=args.shock_prior_probability,
    )
    if args.mode == "map":
        args.output_dir.mkdir(parents=True, exist_ok=True)
        result = fit_pats_map(counts, config)
        (args.output_dir / "pats_map_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    elif args.mode in {"development", "test"}:
        result = run_simulation_bank(counts, args.mode, args.output_dir, args.workers, config)
    elif args.mode == "mcmc":
        result = fit_pats_mcmc(counts, args.output_dir, args.draws, args.tune, config=config)
    else:
        result = run_approximation_check(counts, args.output_dir, args.draws, args.tune, config)
    print(json.dumps(result, indent=2))

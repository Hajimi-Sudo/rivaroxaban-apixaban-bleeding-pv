"""Statistical utilities shared by the pharmacovigilance experiments."""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
from scipy.optimize import brentq
from scipy.stats import fisher_exact, norm


def corrected_cells(a: int, b: int, c: int, d: int) -> tuple[np.ndarray, bool]:
    cells = np.asarray([a, b, c, d], dtype=float)
    corrected = bool(np.any(cells == 0))
    if corrected:
        cells += 0.5
    return cells, corrected


def two_by_two_statistics(a: int, b: int, c: int, d: int) -> dict[str, float | int | bool]:
    """Return ROR, PRR, intervals and Fisher exact p for an active comparator table."""
    cells, corrected = corrected_cells(a, b, c, d)
    aa, bb, cc, dd = cells
    ror = (aa * dd) / (bb * cc)
    se_ror = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    prr = (aa / (aa + bb)) / (cc / (cc + dd))
    se_prr = math.sqrt(1 / aa - 1 / (aa + bb) + 1 / cc - 1 / (cc + dd))
    _, fisher_p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return {
        "a_exposed_event": int(a),
        "b_exposed_non_event": int(b),
        "c_comparator_event": int(c),
        "d_comparator_non_event": int(d),
        "ror": float(ror),
        "ror_ci95_lower": float(math.exp(math.log(ror) - 1.96 * se_ror)),
        "ror_ci95_upper": float(math.exp(math.log(ror) + 1.96 * se_ror)),
        "prr": float(prr),
        "prr_ci95_lower": float(math.exp(math.log(prr) - 1.96 * se_prr)),
        "prr_ci95_upper": float(math.exp(math.log(prr) + 1.96 * se_prr)),
        "fisher_exact_p": float(fisher_p),
        "jeffreys_0_5_correction": corrected,
        "exposed_event_stable_20": bool(a >= 20),
        "comparator_event_stable_20": bool(c >= 20),
        "displayable_cells_gte_3": bool(a >= 3 and c >= 3),
    }


def benjamini_hochberg(p_values: Iterable[float]) -> list[float]:
    values = np.asarray(list(p_values), dtype=float)
    if values.size == 0:
        return []
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    return restored.tolist()


def approximate_power_for_or(
    odds_ratio: float, n_exposed: int, n_comparator: int, comparator_rate: float, alpha: float = 0.05
) -> float:
    """Two-sided normal-approximation power for two independent proportions."""
    if not (odds_ratio > 0 and n_exposed > 0 and n_comparator > 0 and 0 < comparator_rate < 1):
        return float("nan")
    p0 = comparator_rate
    p1 = odds_ratio * p0 / (1 - p0 + odds_ratio * p0)
    pooled = (n_exposed * p1 + n_comparator * p0) / (n_exposed + n_comparator)
    null_se = math.sqrt(pooled * (1 - pooled) * (1 / n_exposed + 1 / n_comparator))
    alt_se = math.sqrt(p1 * (1 - p1) / n_exposed + p0 * (1 - p0) / n_comparator)
    z_alpha = norm.ppf(1 - alpha / 2)
    delta = abs(p1 - p0)
    return float(norm.cdf((delta - z_alpha * null_se) / alt_se) + norm.cdf((-delta - z_alpha * null_se) / alt_se))


def minimum_detectable_or(
    n_exposed: int, n_comparator: int, comparator_events: int, target_power: float = 0.80
) -> float:
    if n_exposed <= 0 or n_comparator <= 0 or not (0 < comparator_events < n_comparator):
        return float("nan")
    p0 = comparator_events / n_comparator

    def objective(log_or: float) -> float:
        return approximate_power_for_or(math.exp(log_or), n_exposed, n_comparator, p0) - target_power

    lower, upper = math.log(1.000001), math.log(100.0)
    if objective(upper) < 0:
        return float("inf")
    return float(math.exp(brentq(objective, lower, upper)))

"""Render revision figures from aggregate outputs only.

This module performs no case-level analysis. It reads the aggregate CSV/JSON
files produced by ``faers_full.py`` and ``external_validation.py`` and writes
publication-ready PNG figures for the revised manuscript.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RIV = "#D55E00"
API = "#0072B2"
GREEN = "#009E73"
AMBER = "#E69F00"
GREY = "#666666"
LIGHT = "#CCCCCC"


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
            "font.size": 10,
            "axes.labelsize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def _clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(LIGHT)
    ax.spines["bottom"].set_color(LIGHT)
    ax.tick_params(direction="out", color=LIGHT)


def _save(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _forest(frame: pd.DataFrame, out_dir: Path, stem: str, xlim: tuple[float, float], show_q: bool = False) -> None:
    frame = frame.reset_index(drop=True)
    y = np.arange(len(frame))[::-1]
    fig = plt.figure(figsize=(8.8, max(4.6, 0.47 * len(frame) + 1.6)))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.0, 1.55], wspace=0.03)
    ax = fig.add_subplot(gs[0, 0])
    txt = fig.add_subplot(gs[0, 1], sharey=ax)
    ax.axvline(1.0, color=GREY, linestyle="--", linewidth=1.0, zorder=0)

    for i, row in frame.iterrows():
        yi = y[i]
        if not np.isfinite(row["ror"]):
            continue
        face = "white" if row["sparse"] else row["color"]
        ax.errorbar(
            row["ror"], yi,
            xerr=[[row["ror"] - row["lower"]], [row["upper"] - row["ror"]]],
            fmt="o", markersize=6.2, markerfacecolor=face, markeredgecolor=row["color"],
            markeredgewidth=1.4, ecolor=row["color"], elinewidth=1.25, capsize=2.5, zorder=3,
        )
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_yticks(y)
    ax.set_yticklabels(frame["label"])
    ax.set_xlabel("Reporting odds ratio (log scale)", fontweight="bold")
    ax.grid(axis="x", which="major", linestyle=":", linewidth=0.7, color=LIGHT)
    _clean_axes(ax)

    txt.set_xlim(0, 1)
    txt.set_xticks([])
    txt.set_yticks(y)
    txt.tick_params(left=False, labelleft=False)
    for spine in txt.spines.values():
        spine.set_visible(False)
    txt.text(0.0, y.max() + 0.85, "ROR (95% CI)" + ("        BH q" if show_q else ""), fontsize=9.5, fontweight="bold", va="center")
    for i, row in frame.iterrows():
        if np.isfinite(row["ror"]):
            value = f"{row['ror']:.2f} ({row['lower']:.2f}–{row['upper']:.2f})"
            if show_q:
                q = row["q"]
                value += f"    {q:.3f}" if q >= 0.001 else f"    {q:.1e}"
        else:
            value = "Not estimable (2 vs 0 events)"
        txt.text(0.0, y[i], value, fontsize=8.8, va="center")
    ax.set_ylim(-0.8, len(frame) - 0.05)
    fig.subplots_adjust(left=0.29, right=0.985, bottom=0.12, top=0.96)
    _save(fig, out_dir, stem)


def render_decomposition(faers_dir: Path, out_dir: Path) -> None:
    data = pd.read_csv(faers_dir / "post_result_phenotype_decomposition.csv")
    names = {
        "post_result_pt_abnormal_uterine_bleeding": "Abnormal uterine bleeding",
        "post_result_pt_heavy_menstrual_bleeding": "Heavy menstrual bleeding",
        "post_result_pt_intermenstrual_bleeding": "Intermenstrual bleeding",
        "post_result_pt_uterine_haemorrhage": "Uterine haemorrhage",
        "post_result_pt_vaginal_haemorrhage": "Vaginal haemorrhage",
        "post_result_menstrual_specific": "Menstrual-specific grouping",
        "post_result_uterine_specific": "Uterine-specific grouping",
        "post_result_broad_excluding_vaginal_haemorrhage": "Broad excluding vaginal haemorrhage",
    }
    rows = []
    for _, row in data.iterrows():
        analysis = row["analysis"]
        if analysis not in names:
            continue
        sparse = not (bool(row["exposed_event_stable_20"]) and bool(row["comparator_event_stable_20"]))
        non_estimable = analysis == "post_result_pt_abnormal_uterine_bleeding" and int(row["c_comparator_event"]) == 0
        estimate = float(row["conditional_mle_or"]) if sparse else float(row["ror"])
        lower = float(row["conditional_exact_ci95_lower"]) if sparse else float(row["ror_ci95_lower"])
        upper = float(row["conditional_exact_ci95_upper"]) if sparse else float(row["ror_ci95_upper"])
        rows.append(
            {
                "label": names[analysis],
                "ror": np.nan if non_estimable else estimate,
                "lower": np.nan if non_estimable else lower,
                "upper": np.nan if non_estimable else upper,
                "q": float(row["bh_fdr_q"]),
                "color": RIV if float(row["ror"]) > 1 else API,
                "sparse": sparse,
            }
        )
    _forest(pd.DataFrame(rows), out_dir, "figure3_outcome_definition_forest", (0.08, 64.0), show_q=True)


def _row_from_series(label: str, row: pd.Series | dict[str, object], color: str) -> dict[str, object]:
    sparse = int(row["a_exposed_event"]) < 20 or int(row["c_comparator_event"]) < 20
    if sparse:
        lower = float(row["conditional_exact_ci95_lower"])
        upper = float(row["conditional_exact_ci95_upper"])
        estimate = float(row["conditional_mle_or"])
    else:
        lower = float(row["ror_ci95_lower"])
        upper = float(row["ror_ci95_upper"])
        estimate = float(row["ror"])
    return {
        "label": label,
        "ror": estimate,
        "lower": lower,
        "upper": upper,
        "color": color,
        "sparse": sparse,
    }


def render_main_forest(faers_dir: Path, external_dir: Path, out_dir: Path) -> None:
    primary = json.loads((faers_dir / "primary_contrast.json").read_text(encoding="utf-8"))
    sensitivity = pd.read_csv(faers_dir / "family_C_sensitivities.csv").set_index("analysis")
    external = pd.read_csv(external_dir / "family_B_external_replication.csv").set_index("database")
    specs = [
        ("FAERS primary", primary, RIV),
        ("All suspect roles", sensitivity.loc["all_suspect"], API),
        ("Age 12–55 years", sensitivity.loc["age_12_55"], API),
        ("All female reports", sensitivity.loc["all_females"], API),
        ("VTE indication", sensitivity.loc["vte_only"], API),
        ("Atrial fibrillation indication", sensitivity.loc["atrial_fibrillation_only"], API),
        ("Healthcare-professional reports", sensitivity.loc["healthcare_professional_reporter"], API),
        ("Consumer reports", sensitivity.loc["consumer_reporter"], API),
        ("Narrow AUB/HMB definition", sensitivity.loc["narrow_pt"], AMBER),
        ("Canada Vigilance", external.loc["Canada Vigilance"], GREEN),
        ("JADER (20s–40s proxy)", external.loc["JADER"], GREEN),
    ]
    rows = [_row_from_series(label, row, color) for label, row, color in specs]
    _forest(pd.DataFrame(rows), out_dir, "figure2_primary_sensitivity_contextual_forest", (0.25, 20.0))


def render_quarterly(faers_dir: Path, out_dir: Path) -> None:
    q = pd.read_csv(faers_dir / "quarterly_active_comparator_counts.csv")
    x = np.arange(len(q))
    riv_prop = 100 * q["rivaroxaban_aub"] / q["rivaroxaban_reports"]
    api_prop = 100 * q["apixaban_aub"] / q["apixaban_reports"]
    fig, axes = plt.subplots(2, 1, figsize=(8.6, 6.2), sharex=True, gridspec_kw={"hspace": 0.12})
    axes[0].plot(x, q["rivaroxaban_reports"], color=RIV, marker="o", markersize=3.3, linewidth=1.25, label="Rivaroxaban")
    axes[0].plot(x, q["apixaban_reports"], color=API, marker="s", markersize=3.0, linewidth=1.25, label="Apixaban")
    axes[0].set_ylabel("Eligible reports", fontweight="bold")
    axes[0].legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.14))
    axes[1].plot(x, riv_prop, color=RIV, marker="o", markersize=3.3, linewidth=1.25)
    axes[1].plot(x, api_prop, color=API, marker="s", markersize=3.0, linewidth=1.25)
    axes[1].set_ylabel("Broad outcome reports (%)", fontweight="bold")
    axes[1].set_xlabel("FAERS quarter", fontweight="bold")
    for ax in axes:
        for period in ("2021Q1", "2021Q2"):
            idx = int(q.index[q["period"].eq(period)][0])
            ax.axvline(idx, color=GREY, linestyle=":", linewidth=0.8, zorder=0)
        ax.grid(axis="y", linestyle=":", linewidth=0.6, color=LIGHT)
        _clean_axes(ax)
    ticks = list(range(0, len(q), 4))
    if ticks[-1] != len(q) - 1:
        ticks.append(len(q) - 1)
    axes[1].set_xticks(ticks)
    axes[1].set_xticklabels(q.loc[ticks, "period"], rotation=45, ha="right")
    axes[0].text(-0.075, 1.02, "(A)", transform=axes[0].transAxes, fontweight="bold")
    axes[1].text(-0.075, 1.02, "(B)", transform=axes[1].transAxes, fontweight="bold")
    fig.subplots_adjust(left=0.11, right=0.985, top=0.93, bottom=0.15)
    _save(fig, out_dir, "figure4_quarterly_reporting_patterns")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--faers-dir", type=Path, required=True)
    parser.add_argument("--external-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    _style()
    render_main_forest(args.faers_dir, args.external_dir, args.out_dir)
    render_decomposition(args.faers_dir, args.out_dir)
    render_quarterly(args.faers_dir, args.out_dir)
    audit = {
        "sources": [
            str(args.faers_dir / "primary_contrast.json"),
            str(args.faers_dir / "family_C_sensitivities.csv"),
            str(args.faers_dir / "post_result_phenotype_decomposition.csv"),
            str(args.faers_dir / "quarterly_active_comparator_counts.csv"),
            str(args.external_dir / "family_B_external_replication.csv"),
        ],
        "outputs": [
            "figure2_primary_sensitivity_contextual_forest.png",
            "figure3_outcome_definition_forest.png",
            "figure4_quarterly_reporting_patterns.png",
        ],
        "abnormal_uterine_bleeding_rule": "No point estimate displayed when a comparator event cell is zero.",
        "manual_values": False,
    }
    (args.out_dir / "figure_render_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

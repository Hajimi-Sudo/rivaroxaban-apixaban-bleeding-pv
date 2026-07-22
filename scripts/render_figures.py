"""Render manuscript figures directly from frozen pharmacovigilance outputs."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
AGGREGATE = ROOT / "data" / "aggregate"
FAERS = AGGREGATE
EXTERNAL = AGGREGATE
OUT = ROOT / "figures"

RIV = "#D55E00"
API = "#0072B2"
GREEN = "#009E73"
AMBER = "#E69F00"
GREY = "#666666"
LIGHT = "#CCCCCC"

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
        "font.size": 10,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
    }
)


def clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(LIGHT)
    ax.spines["bottom"].set_color(LIGHT)
    ax.tick_params(direction="out", color=LIGHT)


def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.png", format="png", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def figure1_flow() -> None:
    flow = pd.read_csv(FAERS / "cohort_flow.csv").set_index("stage")["n"].astype(int)
    primary = json.loads((FAERS / "primary_contrast.json").read_text(encoding="utf-8"))
    assert primary["a_exposed_event"] + primary["b_exposed_non_event"] == flow["exclusive_rivaroxaban_primary"]
    assert primary["c_comparator_event"] + primary["d_comparator_non_event"] == flow["exclusive_apixaban_primary"]

    fig, ax = plt.subplots(figsize=(7.2, 7.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.axis("off")

    boxes = [
        (5, 11.0, "Raw FAERS demographic rows", flow["raw_demo_rows"], GREY),
        (5, 9.4, "Deduplicated latest cases", flow["deduplicated_cases"], GREY),
        (5, 7.8, "Female reports", flow["female_cases"], GREY),
        (5, 6.2, "Female reports with convertible age", flow["female_convertible_age"], GREY),
        (5, 4.6, "Females aged 15–49 years", flow["female_age_15_49"], GREY),
        (2.6, 2.7, "Exclusive primary-suspect\nrivaroxaban reports", flow["exclusive_rivaroxaban_primary"], RIV),
        (7.4, 2.7, "Exclusive primary-suspect\napixaban reports", flow["exclusive_apixaban_primary"], API),
        (2.6, 0.8, "Broad-phenotype reports", primary["a_exposed_event"], RIV),
        (7.4, 0.8, "Broad-phenotype reports", primary["c_comparator_event"], API),
    ]
    widths = {11.0: 4.6, 9.4: 4.6, 7.8: 4.6, 6.2: 4.6, 4.6: 4.6, 2.7: 4.0, 0.8: 3.7}
    for x, y, label, value, color in boxes:
        width = widths[y]
        patch = FancyBboxPatch(
            (x - width / 2, y - 0.55), width, 1.1,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            facecolor="white", edgecolor=color, linewidth=1.6,
        )
        ax.add_patch(patch)
        ax.text(x, y + 0.13, label, ha="center", va="center", fontsize=9.5)
        ax.text(x, y - 0.25, f"n = {value:,}", ha="center", va="center", fontsize=10, fontweight="bold", color=color)

    arrows = [
        ((5, 10.43), (5, 9.97)), ((5, 8.83), (5, 8.37)),
        ((5, 7.23), (5, 6.77)), ((5, 5.63), (5, 5.17)),
        ((4.45, 4.03), (2.95, 3.27)), ((5.55, 4.03), (7.05, 3.27)),
        ((2.6, 2.13), (2.6, 1.37)), ((7.4, 2.13), (7.4, 1.37)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "-|>", "color": LIGHT, "lw": 1.3})
    save(fig, "figure1_study_flow")


def forest_plot(
    frame: pd.DataFrame,
    stem: str,
    xlim: tuple[float, float],
    show_q: bool = False,
) -> None:
    frame = frame.reset_index(drop=True)
    y = np.arange(len(frame))[::-1]
    fig = plt.figure(figsize=(8.8, max(4.6, 0.47 * len(frame) + 1.6)))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.0, 1.45], wspace=0.03)
    ax = fig.add_subplot(gs[0, 0])
    txt = fig.add_subplot(gs[0, 1], sharey=ax)

    ax.axvline(1.0, color=GREY, linestyle="--", linewidth=1.0, zorder=0)
    for i, row in frame.iterrows():
        yi = y[i]
        color = row["color"]
        face = "white" if row.get("sparse", False) else color
        ax.errorbar(
            row["ror"], yi,
            xerr=[[row["ror"] - row["lower"]], [row["upper"] - row["ror"]]],
            fmt="o", markersize=6.2, markerfacecolor=face, markeredgecolor=color,
            markeredgewidth=1.4, ecolor=color, elinewidth=1.25, capsize=2.5, zorder=3,
        )
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_yticks(y)
    ax.set_yticklabels(frame["label"])
    ax.set_xlabel("Reporting odds ratio (log scale)", fontweight="bold")
    ax.grid(axis="x", which="major", linestyle=":", linewidth=0.7, color=LIGHT)
    clean_axes(ax)

    txt.set_xlim(0, 1)
    txt.set_xticks([])
    txt.set_yticks(y)
    txt.tick_params(left=False, labelleft=False)
    for spine in txt.spines.values():
        spine.set_visible(False)
    header = "ROR (95% CI)" + ("        BH q" if show_q else "")
    txt.text(0.0, y.max() + 0.85, header, fontsize=9.5, fontweight="bold", va="center")
    for i, row in frame.iterrows():
        value = f"{row['ror']:.2f} ({row['lower']:.2f}–{row['upper']:.2f})"
        if show_q:
            q = row["q"]
            qtext = f"{q:.3f}" if q >= 0.001 else f"{q:.1e}"
            value += f"    {qtext}"
        txt.text(0.0, y[i], value, fontsize=8.8, va="center")

    ax.set_ylim(-0.8, len(frame) - 0.05)
    fig.subplots_adjust(left=0.29, right=0.985, bottom=0.12, top=0.96)
    save(fig, stem)


def figure2_main_forest() -> None:
    primary = json.loads((FAERS / "primary_contrast.json").read_text(encoding="utf-8"))
    sens = pd.read_csv(FAERS / "family_C_sensitivities.csv").set_index("analysis")
    ext = pd.read_csv(EXTERNAL / "family_B_external_replication.csv").set_index("database")
    labels = [
        ("FAERS primary", primary, RIV, False),
        ("All suspect roles", sens.loc["all_suspect"], API, False),
        ("Age 12–55 years", sens.loc["age_12_55"], API, False),
        ("All female reports", sens.loc["all_females"], API, False),
        ("VTE indication", sens.loc["vte_only"], API, False),
        ("Atrial fibrillation indication", sens.loc["atrial_fibrillation_only"], API, True),
        ("Healthcare-professional reports", sens.loc["healthcare_professional_reporter"], API, False),
        ("Consumer reports", sens.loc["consumer_reporter"], API, False),
        ("Excluding 2021 Q1–Q2", sens.loc["exclude_2021q1_q2_safety_communication"], API, False),
        ("Narrow AUB/HMB phenotype", sens.loc["narrow_pt"], AMBER, False),
        ("Canada Vigilance", ext.loc["Canada Vigilance"], GREEN, True),
        ("JADER", ext.loc["JADER"], GREEN, True),
    ]
    rows = []
    for label, obj, color, sparse in labels:
        rows.append(
            {
                "label": label,
                "ror": float(obj["ror"]),
                "lower": float(obj["ror_ci95_lower"]),
                "upper": float(obj["ror_ci95_upper"]),
                "color": color,
                "sparse": sparse,
            }
        )
    forest_plot(pd.DataFrame(rows), "figure2_primary_sensitivity_external_forest", (0.25, 16.0))


def figure3_phenotype_forest() -> None:
    phen = pd.read_csv(FAERS / "post_result_phenotype_decomposition.csv")
    names = {
        "post_result_pt_abnormal_uterine_bleeding": "Abnormal uterine bleeding",
        "post_result_pt_heavy_menstrual_bleeding": "Heavy menstrual bleeding",
        "post_result_pt_intermenstrual_bleeding": "Intermenstrual bleeding",
        "post_result_pt_uterine_haemorrhage": "Uterine haemorrhage",
        "post_result_pt_vaginal_haemorrhage": "Vaginal haemorrhage",
        "post_result_menstrual_specific": "Menstrual-specific group",
        "post_result_uterine_specific": "Uterine-specific group",
        "post_result_vaginal_haemorrhage_only": "Vaginal-haemorrhage group",
        "post_result_broad_excluding_vaginal_haemorrhage": "Broad excluding vaginal haemorrhage",
    }
    rows = []
    for _, row in phen.iterrows():
        rows.append(
            {
                "label": names[row["analysis"]],
                "ror": row["ror"],
                "lower": row["ror_ci95_lower"],
                "upper": row["ror_ci95_upper"],
                "q": row["bh_fdr_q"],
                "color": RIV if row["ror"] > 1 else API,
                "sparse": not (bool(row["exposed_event_stable_20"]) and bool(row["comparator_event_stable_20"])),
            }
        )
    forest_plot(pd.DataFrame(rows), "figure3_post_result_phenotype_forest", (0.08, 64.0), show_q=True)


def figure4_quarterly() -> None:
    q = pd.read_csv(FAERS / "quarterly_active_comparator_counts.csv")
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
    axes[1].set_ylabel("Broad phenotype (%)", fontweight="bold")
    axes[1].set_xlabel("FAERS quarter", fontweight="bold")

    for ax in axes:
        for period in ("2021Q1", "2021Q2"):
            idx = int(q.index[q["period"].eq(period)][0])
            ax.axvline(idx, color=GREY, linestyle=":", linewidth=0.8, zorder=0)
        ax.grid(axis="y", linestyle=":", linewidth=0.6, color=LIGHT)
        clean_axes(ax)

    tick_idx = list(range(0, len(q), 4))
    if tick_idx[-1] != len(q) - 1:
        tick_idx.append(len(q) - 1)
    axes[1].set_xticks(tick_idx)
    axes[1].set_xticklabels(q.loc[tick_idx, "period"], rotation=45, ha="right")
    axes[0].text(-0.075, 1.02, "(A)", transform=axes[0].transAxes, fontweight="bold")
    axes[1].text(-0.075, 1.02, "(B)", transform=axes[1].transAxes, fontweight="bold")
    fig.subplots_adjust(left=0.11, right=0.985, top=0.93, bottom=0.15)
    save(fig, "figure4_quarterly_reporting_patterns")


def write_audit() -> None:
    audit = {
        "source_files": [
            str(FAERS / "cohort_flow.csv"),
            str(FAERS / "primary_contrast.json"),
            str(FAERS / "quarterly_active_comparator_counts.csv"),
            str(FAERS / "family_C_sensitivities.csv"),
            str(FAERS / "post_result_phenotype_decomposition.csv"),
            str(EXTERNAL / "family_B_external_replication.csv"),
        ],
        "outputs": [
            "figure1_study_workflow.png",
            "figure2_primary_sensitivity_external_forest.png",
            "figure3_post_result_phenotype_forest.png",
            "figure4_quarterly_reporting_patterns.png",
        ],
        "formats": ["png"],
        "figure1_provenance": "User-provided JPG converted to PNG without content changes",
        "manual_values": False,
        "palette": "Okabe-Ito colorblind-safe",
    }
    (OUT / "figure_render_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


def main() -> None:
    # Figure 1 is a user-provided workflow graphic retained as PNG.
    figure2_main_forest()
    figure3_phenotype_forest()
    figure4_quarterly()
    write_audit()


if __name__ == "__main__":
    main()

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


OUT = Path(__file__).resolve().parents[1] / "figures" / "figure1_numeric_report_flow.png"

steps = [
    ("Raw FAERS DEMO rows", "18,159,607"),
    ("Deduplicated longitudinal cases", "15,508,934"),
    ("Cases coded as female", "7,990,400"),
    ("Female cases with convertible age", "5,077,192"),
    ("Female cases aged 15–49 years", "1,595,107"),
]

fig, ax = plt.subplots(figsize=(10.5, 7.0), dpi=240)
ax.set_xlim(0, 10.5)
ax.set_ylim(0, 10)
ax.axis("off")

blue = "#1F5A94"
light = "#EEF5FB"
dark = "#17324D"
ys = [9.0, 7.55, 6.10, 4.65, 3.20]

for index, ((label, count), y) in enumerate(zip(steps, ys)):
    box = FancyBboxPatch(
        (2.0, y - 0.45),
        6.5,
        0.9,
        boxstyle="round,pad=0.03,rounding_size=0.08",
        linewidth=1.4,
        edgecolor=blue,
        facecolor=light,
    )
    ax.add_patch(box)
    ax.text(2.25, y + 0.12, label, ha="left", va="center", fontsize=11, color=dark)
    ax.text(8.18, y + 0.12, count, ha="right", va="center", fontsize=11, fontweight="bold", color=dark)
    if index < len(ys) - 1:
        ax.add_patch(
            FancyArrowPatch(
                (5.25, y - 0.46),
                (5.25, ys[index + 1] + 0.46),
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.2,
                color=blue,
            )
        )

ax.add_patch(
    FancyArrowPatch(
        (5.25, 2.72), (3.05, 1.95), arrowstyle="-|>", mutation_scale=13, linewidth=1.2, color=blue
    )
)
ax.add_patch(
    FancyArrowPatch(
        (5.25, 2.72), (7.45, 1.95), arrowstyle="-|>", mutation_scale=13, linewidth=1.2, color=blue
    )
)

branches = [
    (0.65, "Exclusive rivaroxaban\nprimary-suspect reports", "4,193"),
    (5.15, "Exclusive apixaban\nprimary-suspect reports", "1,985"),
]
for x, label, count in branches:
    box = FancyBboxPatch(
        (x, 0.65),
        4.7,
        1.25,
        boxstyle="round,pad=0.03,rounding_size=0.08",
        linewidth=1.6,
        edgecolor=blue,
        facecolor="#DCEBFA",
    )
    ax.add_patch(box)
    ax.text(x + 0.25, 1.35, label, ha="left", va="center", fontsize=10.5, color=dark, linespacing=1.25)
    ax.text(x + 4.4, 0.94, count, ha="right", va="center", fontsize=12, fontweight="bold", color=dark)

ax.text(
    5.25,
    0.16,
    "Dual primary-suspect index-drug reports excluded from the active comparison: 0.",
    ha="center",
    va="center",
    fontsize=9.2,
    color="#465A6E",
)

fig.tight_layout(pad=0.4)
fig.savefig(OUT, bbox_inches="tight", facecolor="white")
plt.close(fig)

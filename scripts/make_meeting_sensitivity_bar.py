"""Top-10 method-parameter effect sizes (Phase 5a) as a bar chart for the 2026-11-13 meeting deck.

Read-only on results/sensitivity/ranking.csv. Effect = max |delta crop-macro mean F1@0.50| over the
perturbations, BA-TOS_sel route (12 dev crops, seed-ignore). Evaluation switches (kind == 'eval') are
excluded: they are protocol choices, not method parameters.

Usage: python -I scripts/make_meeting_sensitivity_bar.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results" / "sensitivity" / "ranking.csv"
OUT = ROOT / "paper" / "slides_meeting" / "presentations" / "felipe" / "batos" / "images" / "sens_top10.png"

# Categorical slots 1-3 of the dataviz reference palette (validated: light surface, CVD-separated).
COLORS = {"selector": "#2a78d6", "openrec_flatten": "#eb6834", "markers": "#1baf7a"}
LABELS = {"selector": "ToS selector", "openrec_flatten": "openrec / flatten", "markers": "markers"}
TEXT, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"

df = pd.read_csv(SRC)
df = df[(df["kind"] != "eval") & (df["complete"] == True)]  # noqa: E712
top = df.sort_values("sel_max_abs_delta_mean", ascending=False).head(10).iloc[::-1]
missing = set(top["group"]) - set(COLORS)
if missing:
    raise SystemExit(f"no colour for group(s) {missing}")

names = [p.replace("selector_", "").replace("openrec_", "openrec ").replace("flatten_", "flatten ") for p in top["param"]]
vals = top["sel_max_abs_delta_mean"].to_numpy()

fig, ax = plt.subplots(figsize=(6.4, 3.9), dpi=300)
ax.barh(names, vals, height=0.62, color=[COLORS[g] for g in top["group"]], edgecolor="white", linewidth=1.5)
for y, v in enumerate(vals):
    ax.text(v + 0.012, y, f"{v:.2f}", va="center", ha="left", fontsize=8, color=TEXT)
ax.set_xlim(0, 1.0)
ax.set_xlabel("max |Δ crop-macro F1@0.50|, one parameter perturbed (×0.5, ×2 or ±1)", fontsize=8, color=MUTED)
ax.tick_params(axis="y", labelsize=8, colors=TEXT, length=0)
ax.tick_params(axis="x", labelsize=7, colors=MUTED)
ax.grid(axis="x", color=GRID, linewidth=0.6)
ax.set_axisbelow(True)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(GRID)
handles = [plt.Rectangle((0, 0), 1, 1, color=COLORS[g]) for g in COLORS]
ax.legend(handles, [LABELS[g] for g in COLORS], loc="lower right", fontsize=7.5, frameon=False)
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, facecolor="white")
print(OUT)
print(top[["param", "group", "frozen_value", "sel_max_abs_delta_mean"]].iloc[::-1].to_string(index=False))

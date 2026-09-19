#!/usr/bin/env python3
"""figures/fig2_fd_vs_gamma.png -- FD vs gamma under three box-counting estimators,
with the one-sided boundary-pixel count as context bars and the zoom runs marked
separately.  Reads the CSVs produced by fd_estimators.py."""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#77766f"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"     # validated categorical slots 1-3
BAR, BAR_EDGE = "#e4e3de", "#cfcec8"

OUT = "figures/fig2_fd_vs_gamma.png"


def load(path):
    with open(path) as fh:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(fh)]


def main():
    full = sorted(load("results/fd_estimators.csv"), key=lambda r: r["gamma"])
    zoom = sorted(load("results/zoom/fd_estimators_zoom.csv"), key=lambda r: r["gamma"])
    g = [r["gamma"] for r in full]

    fig, ax = plt.subplots(figsize=(9.6, 6.0), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    # context bars (secondary axis, deliberately recessive: they are scale, not a series)
    axb = ax.twinx()
    axb.bar(g, [r["n_boundary_px"] for r in full], width=0.016,
            color=BAR, edgecolor=BAR_EDGE, linewidth=0.8, zorder=1)
    axb.set_ylabel("one-sided boundary pixels  (bars)", color=INK3, fontsize=10)
    axb.tick_params(axis="y", colors=INK3, labelsize=9)
    axb.set_ylim(0, 5200)        # keep bars in the lower third: context, not a series
    for s in axb.spines.values():
        s.set_visible(False)
    for r in full:
        axb.annotate(f"{int(r['n_boundary_px'])}", (r["gamma"], r["n_boundary_px"]),
                     textcoords="offset points", xytext=(0, 3), ha="center",
                     fontsize=8, color=INK3, zorder=2)

    # three estimator lines
    ax.plot(g, [r["FD_two_sided"] for r in full], "-o", color=S2, lw=2, ms=8,
            mec=SURFACE, mew=1.5, zorder=4, label="two-sided, ε ≥ 1")
    ax.plot(g, [r["FD_one_sided_eps1"] for r in full], "-o", color=S3, lw=2, ms=8,
            mec=SURFACE, mew=1.5, zorder=4, label="one-sided, ε ≥ 1")
    ax.errorbar(g, [r["FD_headline"] for r in full], yerr=[r["FD_unc"] for r in full],
                fmt="-o", color=S1, lw=2.4, ms=9, mec=SURFACE, mew=1.5,
                ecolor=S1, elinewidth=1.6, capsize=5, capthick=1.6, zorder=5,
                label="headline: one-sided, ε ≥ 2")

    # zoom runs, offset so they do not sit on top of the full-range markers
    zx = [r["gamma"] + 0.012 for r in zoom]
    ax.errorbar(zx, [r["FD_headline"] for r in zoom], yerr=[r["FD_unc"] for r in zoom],
                fmt="D", color=SURFACE, mec=S1, mew=2.2, ms=9,
                ecolor=S1, elinewidth=1.6, capsize=5, capthick=1.6, zorder=6,
                label="headline, zoomed window")
    # the hollow diamond plus the legend entry carries the identity; no per-point text

    # direct labels (relief rule: aqua is under 3:1 on this surface)
    for r, col, txt, dy in ((full[-1], S1, "one-sided, ε≥2", 16),
                            (full[-1], S2, "two-sided, ε≥1", 10),
                            (full[-1], S3, "one-sided, ε≥1", -16)):
        key = {S1: "FD_headline", S2: "FD_two_sided", S3: "FD_one_sided_eps1"}[col]
        ax.annotate(txt, (r["gamma"], r[key]), textcoords="offset points",
                    xytext=(14, dy), ha="left", fontsize=9, color=col, zorder=7)

    ax.set_xlabel(r"damping  $\gamma$", color=INK, fontsize=11)
    ax.set_ylabel("box-counting dimension  FD", color=INK, fontsize=11)
    ax.set_title("Fractal dimension of the inter-well boundary depends on the estimator\n"
                 "as much as on damping", color=INK, fontsize=12.5, pad=12, loc="left")
    ax.set_xlim(-0.022, 0.515)
    ax.set_ylim(0.95, 1.42)
    ax.set_xticks(g)
    ax.set_xticklabels([f"{v:g}" for v in g])
    ax.tick_params(colors=INK2, labelsize=9.5)
    ax.grid(axis="y", color="#e8e7e2", lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#d8d7d2")
    ax.set_zorder(axb.get_zorder() + 1)
    ax.patch.set_visible(False)

    leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.115), ncol=4,
                    fontsize=9, frameon=False)
    for t in leg.get_texts():
        t.set_color(INK2)

    fig.text(0.013, 0.015, "Bars show how much boundary there is to measure: the "
             "γ = 0.25 and 0.30 full-range grids resolve few pixels, which is why "
             "their error bars are wide.",
             fontsize=8.5, color=INK3)
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    os.makedirs("figures", exist_ok=True)
    fig.savefig(OUT, facecolor=SURFACE)
    print("wrote", OUT)


if __name__ == "__main__":
    main()

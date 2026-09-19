#!/usr/bin/env python3
"""figures/fig3_static_bias.png -- category maps under a static bias F0, plus the
measured G_min(F0) against the quasi-static threshold 2/(3*sqrt(3)) - F0."""
import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from duffing_parameter_space import INTRA, SWITCH, REVERT, VACIL, LABELS, COLORS

SURFACE = "#fcfcfb"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#77766f"
S1 = "#2a78d6"
F_STATIC = 2.0 / (3.0 * np.sqrt(3.0))
OUT = "figures/fig3_static_bias.png"


def main():
    paths = sorted(glob.glob("results/static_bias/static_bias_F*_n256.npz"))
    if not paths:
        raise SystemExit("no static-bias .npz files found")
    runs = sorted((np.load(p) for p in paths), key=lambda z: float(z["F0"]))

    cmap = ListedColormap([COLORS[INTRA], COLORS[SWITCH], COLORS[REVERT], COLORS[VACIL]])
    fig = plt.figure(figsize=(13.4, 7.4), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.25], hspace=0.34, wspace=0.28,
                          left=0.055, right=0.975, top=0.86, bottom=0.145)

    for idx, z in enumerate(runs[:4]):
        ax = fig.add_subplot(gs[idx // 2, idx % 2])
        Om, Gs = z["Omega"], z["G"]
        ax.imshow(z["cat"], origin="lower", extent=[Om[0], Om[-1], Gs[0], Gs[-1]],
                  aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5, interpolation="nearest")
        ax.set_title(rf"$F_0$ = {float(z['F0']):.2f}   (wells {float(z['left']):+.3f}, "
                     rf"{float(z['right']):+.3f};  hilltop {float(z['hilltop']):+.3f})",
                     fontsize=9.5, color=INK, pad=5)
        ax.tick_params(colors=INK2, labelsize=8.5)
        if idx % 2 == 0:
            ax.set_ylabel(r"$G$", fontsize=10.5, color=INK)
        if idx // 2 == 1:
            ax.set_xlabel(r"$\Omega$", fontsize=10.5, color=INK)
        for s in ax.spines.values():
            s.set_color("#d8d7d2")

    F0 = np.array([float(z["F0"]) for z in runs])
    gmin = np.array([float(z["gmin"]) for z in runs])
    static = F_STATIC - F0

    axr = fig.add_subplot(gs[:, 2])
    axr.set_facecolor(SURFACE)
    axr.plot(F0, static, ls=(0, (5, 2)), color=INK, lw=2, zorder=3,
             label=r"quasi-static  $2/(3\sqrt{3}) - F_0$")
    axr.plot(F0, gmin, "-o", color=S1, lw=2.4, ms=9, mec=SURFACE, mew=1.5, zorder=4,
             label=r"measured  $G_{\min}(F_0)$")

    for x, y, s in zip(F0, gmin, static):
        axr.annotate(f"{y:.3f}", (x, y), textcoords="offset points", xytext=(0, -16),
                     ha="center", fontsize=8.5, color=S1)
        axr.annotate(f"{100*y/s:.0f}%", (x, (y + s) / 2), textcoords="offset points",
                     xytext=(7, 0), ha="left", va="center", fontsize=8, color=INK3)
        axr.plot([x, x], [y, s], color="#d8d7d2", lw=1, zorder=1)

    axr.text(F0[-1], static[-1], "  quasi-static", ha="left", va="center",
             fontsize=9, color=INK)
    axr.text(F0[-1], gmin[-1], "  measured", ha="left", va="center",
             fontsize=9, color=S1)
    axr.set_xlabel(r"static bias  $F_0$", fontsize=10.5, color=INK)
    axr.set_ylabel(r"forcing amplitude  $G$", fontsize=10.5, color=INK)
    axr.set_title("Harmonic forcing escapes far below the\nquasi-static threshold, "
                  "at every bias", fontsize=10.5, color=INK, pad=8, loc="left")
    axr.set_xlim(-0.02, 0.215)
    axr.set_ylim(0, 0.42)
    axr.tick_params(colors=INK2, labelsize=9)
    axr.grid(axis="y", color="#e8e7e2", lw=0.8)
    axr.set_axisbelow(True)
    for side in ("top", "right"):
        axr.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axr.spines[side].set_color("#d8d7d2")
    legr = axr.legend(loc="lower left", fontsize=8.5, frameon=False)
    for t in legr.get_texts():
        t.set_color(INK2)

    handles = [Patch(facecolor=COLORS[k], edgecolor="#cfcec8", label=LABELS[k])
               for k in (SWITCH, REVERT, VACIL, INTRA)]
    leg = fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9.5,
                     frameon=False, bbox_to_anchor=(0.34, 0.045))
    for t in leg.get_texts():
        t.set_color(INK2)

    fig.suptitle(r"Static bias $F_0$ added to the drive:  "
                 r"$u'' + \gamma u' - u + u^3 = F_0 + G\cos(\Omega\tau)$,  $\gamma$ = 0.07",
                 fontsize=12.5, color=INK, x=0.012, ha="left", y=0.965)
    fig.text(0.012, 0.015, "256 x 256, 500 forcing cycles, started from the left-well "
             "root; hilltop crossing is defined against the middle root of "
             "$u^3 - u = F_0$. Exploration, not a published figure.",
             fontsize=8.5, color=INK3)

    os.makedirs("figures", exist_ok=True)
    fig.savefig(OUT, facecolor=SURFACE)
    print("wrote", OUT)
    for x, y, s in zip(F0, gmin, static):
        print(f"  F0={x:.2f}  G_min={y:.4f}  static={s:.4f}  ratio={y/s:.3f}")


if __name__ == "__main__":
    main()

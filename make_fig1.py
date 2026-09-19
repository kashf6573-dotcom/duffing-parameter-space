#!/usr/bin/env python3
"""figures/fig1_gamma_panels.png -- category maps for the five gammas plus a zoom
of the gamma = 0.07 tongue, cropped from the same full-range grid."""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle

from duffing_parameter_space import INTRA, SWITCH, REVERT, VACIL, LABELS, COLORS

SURFACE = "#fcfcfb"
INK, INK2 = "#0b0b0b", "#52514e"
GAMMAS = [0.001, 0.07, 0.15, 0.25, 0.30]
ZOOM_OF = 0.07
ZOOM_OM, ZOOM_G = (1.05, 1.55), (0.08, 0.22)
OUT = "figures/fig1_gamma_panels.png"


def load(gamma):
    z = np.load(f"results/duffing_gamma{gamma:.3f}_n256.npz")
    return z["cat"], z["Omega"], z["G"]


def draw(ax, cat, Om, Gs, cmap, title):
    ax.imshow(cat, origin="lower", extent=[Om[0], Om[-1], Gs[0], Gs[-1]],
              aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5, interpolation="nearest")
    ax.set_title(title, fontsize=10.5, color=INK, pad=6)
    ax.tick_params(colors=INK2, labelsize=8.5)
    for s in ax.spines.values():
        s.set_color("#d8d7d2")


def main():
    cmap = ListedColormap([COLORS[INTRA], COLORS[SWITCH], COLORS[REVERT], COLORS[VACIL]])
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.6), dpi=200)
    fig.patch.set_facecolor(SURFACE)

    for ax, gamma in zip(axes.ravel(), GAMMAS):
        cat, Om, Gs = load(gamma)
        draw(ax, cat, Om, Gs, cmap, rf"$\gamma$ = {gamma:g}")

    # sixth panel: the gamma = 0.07 tongue, cropped from the same grid
    cat, Om, Gs = load(ZOOM_OF)
    ci = slice(*np.searchsorted(Om, ZOOM_OM))
    ri = slice(*np.searchsorted(Gs, ZOOM_G))
    axz = axes.ravel()[5]
    draw(axz, cat[ri, ci], Om[ci], Gs[ri], cmap,
         rf"$\gamma$ = {ZOOM_OF:g}, tongue detail")

    # outline the crop on the gamma = 0.07 panel it came from
    src = axes.ravel()[GAMMAS.index(ZOOM_OF)]
    src.add_patch(Rectangle((Om[ci][0], Gs[ri][0]),
                            Om[ci][-1] - Om[ci][0], Gs[ri][-1] - Gs[ri][0],
                            fill=False, ec=INK, lw=1.2, ls=(0, (4, 2)), zorder=3))

    for ax in axes[:, 0]:
        ax.set_ylabel(r"$G$", fontsize=11, color=INK)
    for ax in axes[1, :]:
        ax.set_xlabel(r"$\Omega$", fontsize=11, color=INK)

    handles = [Patch(facecolor=COLORS[k], edgecolor="#cfcec8", label=LABELS[k])
               for k in (SWITCH, REVERT, VACIL, INTRA)]
    leg = fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=10,
                     frameon=False, bbox_to_anchor=(0.5, 0.005))
    for t in leg.get_texts():
        t.set_color(INK2)

    fig.suptitle("Behaviour of the bistable Duffing oscillator across the "
                 r"($\Omega$, $G$) plane, by damping",
                 fontsize=13, color=INK, x=0.012, ha="left", y=0.985)
    fig.text(0.012, 0.055, "256 x 256 grids, 500 forcing cycles, initial condition "
             "(-1, 0). Dashed box on the γ = 0.07 panel marks the detail crop.",
             fontsize=8.5, color="#77766f")
    fig.tight_layout(rect=(0, 0.075, 1, 0.965))
    os.makedirs("figures", exist_ok=True)
    fig.savefig(OUT, facecolor=SURFACE)
    print("wrote", OUT)


if __name__ == "__main__":
    main()

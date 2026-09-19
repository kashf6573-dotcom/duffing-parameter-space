#!/usr/bin/env python3
"""Build results/summary.csv from the production .npz files written by duffing_parameter_space.py.

G_min and the category percentages are recomputed from the stored `cat` grid
exactly as main() computes them, so the CSV is derived from the saved data
rather than scraped from the logs.
"""
import csv
import glob
import os
import re

import numpy as np

from duffing_parameter_space import INTRA, SWITCH, REVERT, VACIL
from fd_estimators import analyse

OUT = "results/summary.csv"
PATTERN = "results/duffing_gamma*_n256.npz"


def row(path):
    z = np.load(path)
    cat = z["cat"]
    Gs = z["G"]
    Om = z["Omega"]
    OM, GG = np.meshgrid(Om, Gs)            # same orientation as main()
    inter = cat != INTRA
    gmin = GG[inter].min() if inter.any() else float("nan")
    pct = {k: float(np.mean(cat == k) * 100.0) for k in (INTRA, SWITCH, REVERT, VACIL)}
    fd = analyse(path)                      # headline + disclosure estimators
    return {
        "gamma": float(z["gamma"]),
        "G_min": float(gmin),
        "FD_headline": fd["FD_headline"],
        "FD_window_halfrange": fd["FD_window_halfrange"],
        "FD_two_sided": fd["FD_two_sided"],
        "FD_one_sided_eps1": fd["FD_one_sided_eps1"],
        "n_boundary_px": fd["n_boundary_px"],
        "pct_switching": pct[SWITCH],
        "pct_reverting": pct[REVERT],
        "pct_vacillating": pct[VACIL],
        "pct_intrawell": pct[INTRA],
    }


def main():
    paths = sorted(glob.glob(PATTERN))
    if not paths:
        raise SystemExit(f"no files matching {PATTERN}")
    rows = sorted((row(p) for p in paths), key=lambda r: r["gamma"])
    fields = ["gamma", "G_min", "FD_headline", "FD_window_halfrange", "FD_two_sided",
              "FD_one_sided_eps1", "n_boundary_px", "pct_switching", "pct_reverting",
              "pct_vacillating", "pct_intrawell"]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({
                "gamma": f"{r['gamma']:g}",
                "G_min": f"{r['G_min']:.6f}",
                "FD_headline": f"{r['FD_headline']:.4f}",
                "FD_window_halfrange": f"{r['FD_window_halfrange']:.4f}",
                "FD_two_sided": f"{r['FD_two_sided']:.4f}",
                "FD_one_sided_eps1": f"{r['FD_one_sided_eps1']:.4f}",
                "n_boundary_px": r["n_boundary_px"],
                "pct_switching": f"{r['pct_switching']:.4f}",
                "pct_reverting": f"{r['pct_reverting']:.4f}",
                "pct_vacillating": f"{r['pct_vacillating']:.4f}",
                "pct_intrawell": f"{r['pct_intrawell']:.4f}",
            })
    print(f"wrote {OUT} ({len(rows)} rows)")
    for r in rows:
        print(f"  gamma={r['gamma']:<6g} G_min={r['G_min']:.3f}  "
              f"FD={r['FD_headline']:.3f} +/- {r['FD_window_halfrange']:.3f}")


if __name__ == "__main__":
    main()

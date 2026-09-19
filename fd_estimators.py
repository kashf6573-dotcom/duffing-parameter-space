#!/usr/bin/env python3
"""
Box-counting dimension of the inter-well / intra-well boundary, recomputed from
the saved .npz grids.  Read-only: nothing here re-runs a sweep or rewrites an .npz.

Headline estimator: one-sided boundary, fit over 2 <= eps <= n/4.  The one-sided
mask avoids the factor-of-two inflation of the two-sided mask, and dropping
eps = 1 avoids the saturation of a 1-pixel-thick set at the finest box size.

Fit-window sensitivity: the headline slope is refit over every contiguous window of >= 4
box sizes inside the fit range; FD is reported as mean +/- half-range over those
windows, which measures how much the answer depends on where the scaling window
is placed.
"""
import argparse
import csv
import glob
import os

import numpy as np

from duffing_parameter_space import INTRA, boundary_mask, box_counting_dimension

DEFAULT_PATTERN = "results/duffing_gamma*_n256.npz"
OUT = "results/fd_estimators.csv"
MIN_WINDOW = 4


def _slope(sizes, counts):
    return np.polyfit(np.log(1.0 / sizes), np.log(counts), 1)[0]


def window_spread(sizes, counts, min_eps, min_window=MIN_WINDOW):
    """FD as mean +/- half-range over every contiguous window of >= min_window sizes."""
    keep = (counts > 0) & (sizes >= min_eps)
    s, c = sizes[keep], counts[keep]
    slopes = [_slope(s[i:j], c[i:j])
              for i in range(len(s))
              for j in range(i + min_window, len(s) + 1)]
    if not slopes:
        return float("nan"), float("nan"), 0
    slopes = np.array(slopes)
    return float(slopes.mean()), float((slopes.max() - slopes.min()) / 2.0), len(slopes)


def analyse(path, boundary="one-sided", fit_min_eps=2):
    z = np.load(path)
    cat = z["cat"]
    inter = cat != INTRA

    edge = boundary_mask(inter, boundary)
    fd, sizes, counts = box_counting_dimension(edge, fit_min_eps)
    mean, unc, nwin = window_spread(sizes, counts, fit_min_eps)

    edge2 = boundary_mask(inter, "two-sided")
    fd_two, _, _ = box_counting_dimension(edge2, 1)
    fd_one1, _, _ = box_counting_dimension(edge, 1)

    return {
        "gamma": float(z["gamma"]),
        "FD_headline": mean,
        "FD_window_halfrange": unc,
        "FD_fullrange": float(fd),
        "FD_two_sided": float(fd_two),
        "FD_one_sided_eps1": float(fd_one1),
        "n_boundary_px": int(edge.sum()),
        "n_boundary_px_two_sided": int(edge2.sum()),
        "n_windows": nwin,
        "path": path,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boundary", choices=("two-sided", "one-sided"), default="one-sided")
    ap.add_argument("--fit-min-eps", type=int, default=2)
    ap.add_argument("--glob", default=DEFAULT_PATTERN, help="npz files to analyse")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()

    paths = sorted(glob.glob(a.glob))
    if not paths:
        raise SystemExit(f"no files matching {a.glob}")
    rows = sorted((analyse(p, a.boundary, a.fit_min_eps) for p in paths),
                  key=lambda r: (r["gamma"], r["path"]))

    fields = ["gamma", "FD_headline", "FD_window_halfrange", "FD_fullrange", "FD_two_sided",
              "FD_one_sided_eps1", "n_boundary_px", "n_boundary_px_two_sided", "n_windows"]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{r[k]:g}" if k == "gamma" else
                            f"{r[k]:.4f}" if k.startswith("FD") else r[k]) for k in fields})

    hdr = (f"{'gamma':>7} | {'headline FD':>17} | {'2-sided':>8} {'1s,e>=1':>8} | "
           f"{'px 1s':>6} {'px 2s':>6}")
    print(f"headline: {a.boundary}, eps >= {a.fit_min_eps}, windows of >= {MIN_WINDOW} sizes")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['gamma']:>7g} | {r['FD_headline']:>9.4f} +/- {r['FD_window_halfrange']:<5.4f} | "
              f"{r['FD_two_sided']:>8.4f} {r['FD_one_sided_eps1']:>8.4f} | "
              f"{r['n_boundary_px']:>6d} {r['n_boundary_px_two_sided']:>6d}")
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()

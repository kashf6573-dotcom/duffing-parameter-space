#!/usr/bin/env python3
"""
Reproduction of the forcing-amplitude / driving-frequency parameter-space map of the
dissipative bistable Duffing oscillator, following

    M. N. Hasan, T. E. Greenwood, R. G. Parker, Y. L. Kong, P. Wang,
    "Fractal patterns in the parameter space of a bistable Duffing oscillator",
    Phys. Rev. E 108, L022201 (2023); arXiv:2301.13113.

Governing equation (their Eq. 2, mu = 1):
    u'' + gamma*u' - u + u^3 = G cos(Omega * tau)
Initial condition fixed at (u, u') = (-1, 0) for every simulation.

Protocol (as described in the paper):
  * integrate 500 forcing cycles with fourth-order Runge-Kutta;
  * "intra-well" if u never exceeds the hilltop u = 0 during all 500 cycles;
  * otherwise, using the last 50 cycles as steady state, with
        xi_max = max u, xi_min = min u over those cycles:
        switching  : xi_max > 0 and xi_min > 0   (settled in the other well)
        reverting  : xi_max < 0 and xi_min < 0   (crossed, then returned)
        vacillating: xi_max > 0 and xi_min < 0   (never settles)
  * grid 256 x 256 over 0.80 <= Omega <= 1.8, 0.03 <= G <= 0.30;
  * Hausdorff (box-counting) dimension of the inter-well / intra-well boundary.

The whole grid is integrated simultaneously (vectorized over parameter points), so the
full 256 x 256 x 500-cycle run takes a few minutes in NumPy on a laptop.

Usage:
    python duffing_parameter_space.py                    # gamma = 0.07, 256x256, 500 cycles
    python duffing_parameter_space.py --gamma 0.15
    python duffing_parameter_space.py --n 128 --cycles 300 --spc 120   # quick look
"""
import argparse
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

try:
    from numba import njit, prange
    HAVE_NUMBA = True
except ImportError:                                # pragma: no cover
    HAVE_NUMBA = False

INTRA, SWITCH, REVERT, VACIL = 0, 1, 2, 3
LABELS = {INTRA: "intra-well", SWITCH: "switching", REVERT: "reverting", VACIL: "vacillating"}
# colour convention of the paper's Fig. 3: switching red, reverting blue, vacillating green, intra-well yellow
COLORS = {INTRA: "#F2D64B", SWITCH: "#D7263D", REVERT: "#1F5FBF", VACIL: "#2A9D4B"}


def rhs(u, v, tau, gamma, G, Omega):
    """State derivative for u'' + gamma u' - u + u^3 = G cos(Omega tau)."""
    return v, -gamma * v + u - u**3 + G * np.cos(Omega * tau)


def sweep(gamma, Omega, G, n_cycles=500, steps_per_cycle=200, n_steady=50, u0=-1.0, v0=0.0):
    """
    Vectorized RK4 integration over every (Omega, G) pair.
    Omega, G: flat arrays of equal length. Each point uses its own dt = T / steps_per_cycle
    so that 'cycle' counts are exact for every driving frequency.
    Returns (category, ever_crossed, xi_max, xi_min).
    """
    Omega = np.asarray(Omega, float)
    G = np.asarray(G, float)
    n = Omega.size
    dt = (2.0 * np.pi / Omega) / steps_per_cycle
    u = np.full(n, u0)
    v = np.full(n, v0)
    tau = np.zeros(n)

    ever_crossed = np.zeros(n, bool)          # u > 0 at any time during the run
    xi_max = np.full(n, -np.inf)
    xi_min = np.full(n, np.inf)
    steady_start = n_cycles - n_steady

    total_steps = n_cycles * steps_per_cycle
    h = dt
    for step in range(total_steps):
        k1u, k1v = rhs(u, v, tau, gamma, G, Omega)
        k2u, k2v = rhs(u + 0.5 * h * k1u, v + 0.5 * h * k1v, tau + 0.5 * h, gamma, G, Omega)
        k3u, k3v = rhs(u + 0.5 * h * k2u, v + 0.5 * h * k2v, tau + 0.5 * h, gamma, G, Omega)
        k4u, k4v = rhs(u + h * k3u, v + h * k3v, tau + h, gamma, G, Omega)
        u = u + (h / 6.0) * (k1u + 2 * k2u + 2 * k3u + k4u)
        v = v + (h / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
        tau = tau + h

        ever_crossed |= u > 0.0
        if step >= steady_start * steps_per_cycle:
            np.maximum(xi_max, u, out=xi_max)
            np.minimum(xi_min, u, out=xi_min)

    cat = np.full(n, INTRA, np.int8)
    inter = ever_crossed
    cat[inter & (xi_max > 0) & (xi_min > 0)] = SWITCH
    cat[inter & (xi_max < 0) & (xi_min < 0)] = REVERT
    cat[inter & (xi_max > 0) & (xi_min < 0)] = VACIL
    return cat, ever_crossed, xi_max, xi_min


# ---------------------------------------------------------------------------
# Optional numba fast path.
#
# The scalar kernel below mirrors the NumPy path in `sweep` operation for
# operation, so it is bit-for-bit identical rather than merely close -- which
# matters because the trajectories are chaotic and a 1-ULP difference would be
# amplified into a different category over 100 000 steps.  Two details make
# that work:
#   * `u ** 3.0` (float exponent) compiles to a libm pow() call, matching what
#     NumPy's `u**3` does on a float64 array; numba turns the integer form
#     `u ** 3` into u*u*u, which differs from pow() for ~25% of inputs.
#   * the xi_max / xi_min updates reproduce np.maximum/np.minimum NaN
#     propagation, which a bare `if u > xmax` would not.
# ---------------------------------------------------------------------------
if HAVE_NUMBA:

    @njit(cache=True, parallel=True, fastmath=False, nogil=True)
    def _sweep_kernel(gamma, Omega, G, n_cycles, steps_per_cycle, n_steady, u0, v0):
        n = Omega.size
        cat = np.zeros(n, np.int8)
        ever_crossed = np.zeros(n, np.bool_)
        xi_max = np.empty(n, np.float64)
        xi_min = np.empty(n, np.float64)
        steady_start = n_cycles - n_steady
        thresh = steady_start * steps_per_cycle
        total_steps = n_cycles * steps_per_cycle

        for i in prange(n):
            Om = Omega[i]
            Gi = G[i]
            h = (2.0 * np.pi / Om) / steps_per_cycle
            u = u0
            v = v0
            tau = 0.0
            crossed = False
            xmax = -np.inf
            xmin = np.inf

            for step in range(total_steps):
                k1u = v
                k1v = -gamma * v + u - u ** 3.0 + Gi * np.cos(Om * tau)
                ua = u + 0.5 * h * k1u
                va = v + 0.5 * h * k1v
                ta = tau + 0.5 * h
                k2u = va
                k2v = -gamma * va + ua - ua ** 3.0 + Gi * np.cos(Om * ta)
                ub = u + 0.5 * h * k2u
                vb = v + 0.5 * h * k2v
                k3u = vb
                k3v = -gamma * vb + ub - ub ** 3.0 + Gi * np.cos(Om * ta)
                uc = u + h * k3u
                vc = v + h * k3v
                tc = tau + h
                k4u = vc
                k4v = -gamma * vc + uc - uc ** 3.0 + Gi * np.cos(Om * tc)

                u = u + (h / 6.0) * (k1u + 2 * k2u + 2 * k3u + k4u)
                v = v + (h / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
                tau = tau + h

                if u > 0.0:
                    crossed = True
                if step >= thresh:
                    if u > xmax or np.isnan(u):
                        xmax = u
                    if u < xmin or np.isnan(u):
                        xmin = u

            ever_crossed[i] = crossed
            xi_max[i] = xmax
            xi_min[i] = xmin
            c = INTRA
            if crossed:
                if xmax > 0 and xmin > 0:
                    c = SWITCH
                elif xmax < 0 and xmin < 0:
                    c = REVERT
                elif xmax > 0 and xmin < 0:
                    c = VACIL
            cat[i] = c

        return cat, ever_crossed, xi_max, xi_min


def sweep_numba(gamma, Omega, G, n_cycles=500, steps_per_cycle=200, n_steady=50,
                u0=-1.0, v0=0.0):
    """numba-compiled equivalent of `sweep`; same signature, same return tuple."""
    if not HAVE_NUMBA:
        raise RuntimeError("numba is not available")
    Omega = np.ascontiguousarray(np.asarray(Omega, float))
    G = np.ascontiguousarray(np.asarray(G, float))
    return _sweep_kernel(float(gamma), Omega, G, int(n_cycles), int(steps_per_cycle),
                         int(n_steady), float(u0), float(v0))


def run_sweep(backend, *args, **kwargs):
    """Dispatch to the numba or the reference NumPy implementation."""
    if backend == "numba" or (backend == "auto" and HAVE_NUMBA):
        return sweep_numba(*args, **kwargs)
    return sweep(*args, **kwargs)


def boundary_mask(binary, mode="two-sided"):
    """
    Boundary pixels of a binary image, using 4-neighbourhoods.

    two-sided (default, the original behaviour): every pixel that has a
        neighbour of the other class, so the boundary is 2 pixels thick and
        straddles the interface.
    one-sided: only True (inter-well) pixels that have at least one False
        (intra-well) 4-neighbour, giving a 1-pixel-thick boundary.

    Pixels outside the grid are not treated as neighbours in either mode.
    """
    b = binary.astype(bool)
    edge = np.zeros_like(b)
    if mode == "two-sided":
        edge[:, 1:] |= b[:, 1:] != b[:, :-1]
        edge[:, :-1] |= b[:, 1:] != b[:, :-1]
        edge[1:, :] |= b[1:, :] != b[:-1, :]
        edge[:-1, :] |= b[1:, :] != b[:-1, :]
    elif mode == "one-sided":
        edge[:, 1:] |= b[:, 1:] & ~b[:, :-1]
        edge[:, :-1] |= b[:, :-1] & ~b[:, 1:]
        edge[1:, :] |= b[1:, :] & ~b[:-1, :]
        edge[:-1, :] |= b[:-1, :] & ~b[1:, :]
    else:
        raise ValueError(f"unknown boundary mode {mode!r}")
    return edge


def box_counting_dimension(edge, fit_min_eps=1):
    """
    Least-squares slope of log N(eps) vs log(1/eps) over dyadic box sizes.

    fit_min_eps excludes box sizes below N from the least-squares fit; the
    counts for those sizes are still returned so they can be plotted.  The
    default of 1 keeps every size, i.e. the original behaviour.
    """
    n = edge.shape[0]
    sizes, counts = [], []
    eps = 1
    while eps <= n // 4:
        m = n // eps
        blocks = edge[: m * eps, : m * eps].reshape(m, eps, m, eps).any(axis=(1, 3))
        sizes.append(eps)
        counts.append(blocks.sum())
        eps *= 2
    sizes = np.array(sizes, float)
    counts = np.array(counts, float)
    keep = (counts > 0) & (sizes >= fit_min_eps)
    if keep.sum() < 2:
        raise ValueError(f"fit_min_eps={fit_min_eps} leaves fewer than 2 usable box sizes")
    slope, intercept = np.polyfit(np.log(1.0 / sizes[keep]), np.log(counts[keep]), 1)
    return slope, sizes, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gamma", type=float, default=0.07)
    ap.add_argument("--n", type=int, default=256, help="grid points per axis")
    ap.add_argument("--cycles", type=int, default=500)
    ap.add_argument("--spc", type=int, default=200, help="RK4 steps per forcing cycle")
    ap.add_argument("--omega", type=float, nargs=2, default=(0.80, 1.80))
    ap.add_argument("--G", type=float, nargs=2, default=(0.03, 0.30))
    ap.add_argument("--out", default=None)
    ap.add_argument("--backend", choices=("numpy", "numba", "auto"), default="numpy",
                    help="RK4 implementation; 'numba' is a bit-identical fast path")
    ap.add_argument("--boundary", choices=("two-sided", "one-sided"), default="two-sided",
                    help="boundary extraction for the box count")
    ap.add_argument("--fit-min-eps", type=int, default=1,
                    help="exclude box sizes below this from the least-squares fit")
    a = ap.parse_args()

    Om = np.linspace(a.omega[0], a.omega[1], a.n)
    Gs = np.linspace(a.G[0], a.G[1], a.n)
    OM, GG = np.meshgrid(Om, Gs)            # rows: G, cols: Omega

    t0 = time.time()
    cat, crossed, xi_max, xi_min = run_sweep(a.backend, a.gamma, OM.ravel(), GG.ravel(),
                                             a.cycles, a.spc)
    cat = cat.reshape(a.n, a.n)
    print(f"sweep: {a.n}x{a.n} points, {a.cycles} cycles, gamma={a.gamma}, "
          f"backend={a.backend}: {time.time()-t0:.1f} s")
    for k, name in LABELS.items():
        print(f"  {name:12s} {np.mean(cat == k)*100:5.1f} %")

    inter = cat != INTRA
    edge = boundary_mask(inter, a.boundary)
    fd, sizes, counts = box_counting_dimension(edge, a.fit_min_eps)
    print(f"box-counting dimension of inter/intra boundary: FD = {fd:.3f} "
          f"[{a.boundary}, eps >= {a.fit_min_eps}, {int(edge.sum())} boundary pixels]")
    gmin = GG[inter].min() if inter.any() else float("nan")
    print(f"minimum G with inter-well behaviour (tongue tip): G_min = {gmin:.3f}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    cmap = ListedColormap([COLORS[INTRA], COLORS[SWITCH], COLORS[REVERT], COLORS[VACIL]])
    ext = [a.omega[0], a.omega[1], a.G[0], a.G[1]]
    axes[0].imshow(cat, origin="lower", extent=ext, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5,
                   interpolation="nearest")
    axes[0].set_xlabel(r"$\Omega$"); axes[0].set_ylabel(r"$G$")
    axes[0].set_title(f"behaviour categories, $\\gamma$ = {a.gamma}")
    for k in (SWITCH, REVERT, VACIL, INTRA):
        axes[0].plot([], [], "s", color=COLORS[k], label=LABELS[k])
    axes[0].legend(loc="upper left", fontsize=8, frameon=True)

    axes[1].imshow(edge, origin="lower", extent=ext, aspect="auto", cmap="Greys", interpolation="nearest")
    axes[1].set_xlabel(r"$\Omega$"); axes[1].set_ylabel(r"$G$")
    axes[1].set_title("inter-well / intra-well boundary")

    axes[2].plot(np.log(1.0 / sizes), np.log(counts), "o-")
    axes[2].set_xlabel(r"$\log(1/\epsilon)$"); axes[2].set_ylabel(r"$\log N(\epsilon)$")
    axes[2].set_title(f"box counting: FD = {fd:.3f}")
    fig.tight_layout()

    out = a.out or f"duffing_gamma{a.gamma:.3f}_n{a.n}.png"
    fig.savefig(out, dpi=160)
    np.savez(out.replace(".png", ".npz"), cat=cat, Omega=Om, G=Gs, gamma=a.gamma, fd=fd)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

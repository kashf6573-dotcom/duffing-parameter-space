#!/usr/bin/env python3
"""
Static-bias extension of the Duffing parameter-space sweep:

    u'' + gamma*u' - u + u^3 = F0 + G cos(Omega*tau)

An exploration, not a reproduction of a published figure -- see the
"Extension: static bias" section of the README.

The RK4 kernels below are the ones in duffing_parameter_space.py with two
changes and no others:
  * the constant F0 is added to the right-hand side;
  * the hilltop is the middle equilibrium root of u^3 - u = F0 rather than
    u = 0, and the initial condition is the left-well root rather than -1.
At F0 = 0 the added term is an exact +0.0 and the roots are exactly
(-1, 0, 1), so this file reproduces the baseline sweep bit-for-bit; that is
checked by --selftest.
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from duffing_parameter_space import (INTRA, SWITCH, REVERT, VACIL, LABELS,
                                     HAVE_NUMBA)

if HAVE_NUMBA:
    from numba import njit, prange

F_STATIC = 2.0 / (3.0 * np.sqrt(3.0))       # quasi-static escape threshold, 0.3849


def well_roots(F0):
    """
    The three real roots of u^3 - u = F0, ascending: (left well, hilltop, right well).

    np.roots is polished with Newton so the values are exact to machine precision;
    at F0 = 0 the polish drives the middle root to exactly 0.0, which is what makes
    the F0 = 0 case reduce to the baseline sweep bit-for-bit.
    """
    if abs(F0) >= F_STATIC:
        raise ValueError(f"F0 = {F0} is at or past the static threshold {F_STATIC:.4f}; "
                         "the potential is no longer bistable")
    r = np.sort(np.roots([1.0, 0.0, -1.0, -F0]).real).astype(float)
    for _ in range(3):
        r = r - (r**3 - r - F0) / (3.0 * r**2 - 1.0)
    return r[0], r[1], r[2]


def sweep_bias(gamma, Omega, G, F0, hilltop, u0, n_cycles=500, steps_per_cycle=200,
               n_steady=50, v0=0.0):
    """Vectorized NumPy reference path (structure identical to `sweep`)."""
    Omega = np.asarray(Omega, float)
    G = np.asarray(G, float)
    n = Omega.size
    dt = (2.0 * np.pi / Omega) / steps_per_cycle
    u = np.full(n, u0)
    v = np.full(n, v0)
    tau = np.zeros(n)

    ever_crossed = np.zeros(n, bool)
    xi_max = np.full(n, -np.inf)
    xi_min = np.full(n, np.inf)
    steady_start = n_cycles - n_steady

    total_steps = n_cycles * steps_per_cycle
    h = dt
    for step in range(total_steps):
        k1u = v
        k1v = -gamma * v + u - u**3 + F0 + G * np.cos(Omega * tau)
        ua, va, ta = u + 0.5 * h * k1u, v + 0.5 * h * k1v, tau + 0.5 * h
        k2u = va
        k2v = -gamma * va + ua - ua**3 + F0 + G * np.cos(Omega * ta)
        ub, vb = u + 0.5 * h * k2u, v + 0.5 * h * k2v
        k3u = vb
        k3v = -gamma * vb + ub - ub**3 + F0 + G * np.cos(Omega * ta)
        uc, vc, tc = u + h * k3u, v + h * k3v, tau + h
        k4u = vc
        k4v = -gamma * vc + uc - uc**3 + F0 + G * np.cos(Omega * tc)

        u = u + (h / 6.0) * (k1u + 2 * k2u + 2 * k3u + k4u)
        v = v + (h / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
        tau = tau + h

        ever_crossed |= u > hilltop
        if step >= steady_start * steps_per_cycle:
            np.maximum(xi_max, u, out=xi_max)
            np.minimum(xi_min, u, out=xi_min)

    cat = np.full(n, INTRA, np.int8)
    inter = ever_crossed
    cat[inter & (xi_max > hilltop) & (xi_min > hilltop)] = SWITCH
    cat[inter & (xi_max < hilltop) & (xi_min < hilltop)] = REVERT
    cat[inter & (xi_max > hilltop) & (xi_min < hilltop)] = VACIL
    return cat, ever_crossed, xi_max, xi_min


if HAVE_NUMBA:

    @njit(cache=True, parallel=True, fastmath=False, nogil=True)
    def _sweep_bias_kernel(gamma, Omega, G, F0, hilltop, u0, v0,
                           n_cycles, steps_per_cycle, n_steady):
        n = Omega.size
        cat = np.zeros(n, np.int8)
        ever_crossed = np.zeros(n, np.bool_)
        xi_max = np.empty(n, np.float64)
        xi_min = np.empty(n, np.float64)
        thresh = (n_cycles - n_steady) * steps_per_cycle
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
                k1v = -gamma * v + u - u ** 3.0 + F0 + Gi * np.cos(Om * tau)
                ua = u + 0.5 * h * k1u
                va = v + 0.5 * h * k1v
                ta = tau + 0.5 * h
                k2u = va
                k2v = -gamma * va + ua - ua ** 3.0 + F0 + Gi * np.cos(Om * ta)
                ub = u + 0.5 * h * k2u
                vb = v + 0.5 * h * k2v
                k3u = vb
                k3v = -gamma * vb + ub - ub ** 3.0 + F0 + Gi * np.cos(Om * ta)
                uc = u + h * k3u
                vc = v + h * k3v
                tc = tau + h
                k4u = vc
                k4v = -gamma * vc + uc - uc ** 3.0 + F0 + Gi * np.cos(Om * tc)

                u = u + (h / 6.0) * (k1u + 2 * k2u + 2 * k3u + k4u)
                v = v + (h / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
                tau = tau + h

                if u > hilltop:
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
                if xmax > hilltop and xmin > hilltop:
                    c = SWITCH
                elif xmax < hilltop and xmin < hilltop:
                    c = REVERT
                elif xmax > hilltop and xmin < hilltop:
                    c = VACIL
            cat[i] = c

        return cat, ever_crossed, xi_max, xi_min


def sweep_bias_numba(gamma, Omega, G, F0, hilltop, u0, n_cycles=500,
                     steps_per_cycle=200, n_steady=50, v0=0.0):
    if not HAVE_NUMBA:
        raise RuntimeError("numba is not available")
    Omega = np.ascontiguousarray(np.asarray(Omega, float))
    G = np.ascontiguousarray(np.asarray(G, float))
    return _sweep_bias_kernel(float(gamma), Omega, G, float(F0), float(hilltop),
                              float(u0), float(v0), int(n_cycles),
                              int(steps_per_cycle), int(n_steady))


def run_sweep_bias(backend, *args, **kwargs):
    if backend == "numba" or (backend == "auto" and HAVE_NUMBA):
        return sweep_bias_numba(*args, **kwargs)
    return sweep_bias(*args, **kwargs)


def selftest():
    """F0 = 0 must reproduce the unbiased kernels exactly, and both backends must agree."""
    from duffing_parameter_space import sweep as sweep_base

    left, mid, right = well_roots(0.0)
    ok = (left, mid, right) == (-1.0, 0.0, 1.0)
    print(f"roots at F0=0: {left!r}, {mid!r}, {right!r}  exact: {ok}")

    n = 48
    Om = np.linspace(0.80, 1.80, n)
    Gs = np.linspace(0.03, 0.30, n)
    OM, GG = np.meshgrid(Om, Gs)
    o, g = OM.ravel(), GG.ravel()

    base = sweep_base(0.07, o, g, 120, 200)
    bias_np = sweep_bias(0.07, o, g, 0.0, mid, left, 120, 200)
    bias_nb = run_sweep_bias("numba", 0.07, o, g, 0.0, mid, left, 120, 200)

    def same(a, b):
        return (np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])
                and (a[2].view(np.int64) == b[2].view(np.int64)).all()
                and (a[3].view(np.int64) == b[3].view(np.int64)).all())

    print(f"F0=0 numpy   vs baseline sweep : bit-identical = {same(base, bias_np)}")
    print(f"F0=0 numba   vs baseline sweep : bit-identical = {same(base, bias_nb)}")

    # a nonzero bias must agree between backends too
    l2, m2, _ = well_roots(0.10)
    a = sweep_bias(0.07, o, g, 0.10, m2, l2, 60, 200)
    b = run_sweep_bias("numba", 0.07, o, g, 0.10, m2, l2, 60, 200)
    print(f"F0=0.10 numpy vs numba         : bit-identical = {same(a, b)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--F0", type=float, nargs="+", default=[0.0, 0.05, 0.10, 0.15])
    ap.add_argument("--gamma", type=float, default=0.07)
    ap.add_argument("--n", type=int, default=256)
    ap.add_argument("--cycles", type=int, default=500)
    ap.add_argument("--spc", type=int, default=200)
    ap.add_argument("--omega", type=float, nargs=2, default=(0.80, 1.80))
    ap.add_argument("--G", type=float, nargs=2, default=(0.03, 0.30))
    ap.add_argument("--backend", choices=("numpy", "numba", "auto"), default="numba")
    ap.add_argument("--outdir", default="results/static_bias")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()

    if a.selftest:
        selftest()
        return

    os.makedirs(a.outdir, exist_ok=True)
    Om = np.linspace(a.omega[0], a.omega[1], a.n)
    Gs = np.linspace(a.G[0], a.G[1], a.n)
    OM, GG = np.meshgrid(Om, Gs)

    for F0 in a.F0:
        left, mid, right = well_roots(F0)
        t0 = time.time()
        cat, crossed, xi_max, xi_min = run_sweep_bias(
            a.backend, a.gamma, OM.ravel(), GG.ravel(), F0, mid, left,
            a.cycles, a.spc)
        cat = cat.reshape(a.n, a.n)
        inter = cat != INTRA
        gmin = GG[inter].min() if inter.any() else float("nan")
        static = F_STATIC - F0

        print(f"F0 = {F0:.2f}  wells ({left:+.4f}, {right:+.4f})  hilltop {mid:+.4f}  "
              f"[{time.time()-t0:.1f} s]")
        for k in (SWITCH, REVERT, VACIL, INTRA):
            print(f"    {LABELS[k]:12s} {np.mean(cat == k)*100:6.2f} %")
        print(f"    G_min = {gmin:.4f}   quasi-static 2/(3sqrt3) - F0 = {static:.4f}   "
              f"ratio = {gmin/static:.3f}")

        np.savez(os.path.join(a.outdir, f"static_bias_F{F0:.2f}_n{a.n}.npz"),
                 cat=cat, Omega=Om, G=Gs, gamma=a.gamma, F0=F0,
                 left=left, hilltop=mid, right=right, gmin=gmin, static=static)


if __name__ == "__main__":
    main()

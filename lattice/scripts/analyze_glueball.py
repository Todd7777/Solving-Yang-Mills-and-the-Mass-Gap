#!/usr/bin/env python3
"""Analyse the output of ``glueball.py``: connected 0++ correlators and effective masses for every
smearing level (jackknife errors), Creutz ratios, plaquette and Polyakov loop.

Example:
    python scripts/analyze_glueball.py results/gb_su2_b2.3.npz --md results/gb_su2_b2.3.md
"""
import argparse
import sys

import numpy as np

from ymlattice import exact
from ymlattice import observables as obs
from ymlattice.stats import autocorrelation_time, binned_mean_err, jackknife


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--nblocks", type=int, default=20)
    ap.add_argument("--md", default=None, help="write a markdown summary here")
    args = ap.parse_args()
    d = np.load(args.npz, allow_pickle=True)
    O, W, plaq, poly = d["O"], d["W"], d["plaq"], d["poly"]
    levels, beta, shape, group = d["levels"], float(d["beta"]), tuple(int(x) for x in d["shape"]), str(d["group"])
    N = 2 if group == "su2" else 3
    nmeas, L0 = O.shape[0], O.shape[2]
    lines = []
    out = lines.append
    out(f"# 0++ glueball / Wilson loop analysis: {group.upper()} beta={beta} lattice={shape}")
    out("")
    pm, pe = binned_mean_err(plaq, args.nblocks)
    out(f"- measurements: {nmeas}; plaquette = {pm:.5f} +- {pe:.5f} (tau_int = {autocorrelation_time(plaq):.1f} updates)")
    out(f"- strong-coupling u(beta) = {exact.plaquette_strong_leading(beta, N):.5f}, weak-coupling 1-loop = {exact.plaquette_weak(beta, N):.5f}")
    pr, prerr = binned_mean_err(poly.real, args.nblocks)
    pa, paerr = binned_mean_err(np.abs(poly), args.nblocks)
    out(f"- Polyakov loop: Re <L> = {pr:.4f} +- {prerr:.4f}, <|L|> = {pa:.4f} +- {paerr:.4f} "
        f"(confined phase: Re<L> consistent with 0, |L| ~ 1/sqrt(volume))")
    out(f"- two-loop asymptotic scaling: a Lambda_L = {exact.a_lambda_lattice(beta, N):.3e}")
    out("")
    # ---- Creutz ratios
    out("## Planar Wilson loops and Creutz ratios")
    out("")
    maxloop = W.shape[1]
    Wm = {}
    out("| R x T | W(R,T) |")
    out("|---|---|")
    for R in range(1, maxloop + 1):
        for T in range(R, maxloop + 1):
            m, e = binned_mean_err(W[:, R - 1, T - 1], args.nblocks)
            Wm[(R, T)] = m
            out(f"| {R}x{T} | {m:.5f} +- {e:.5f} |")
    out("")
    out("| chi(R,T) | value |")
    out("|---|---|")
    for R in range(2, maxloop + 1):
        for T in range(R, maxloop + 1):
            def f(samples, R=R, T=T):
                Wd = {k: np.mean(samples[:, k[0] - 1, k[1] - 1]) for k in
                      [(R, T), (R - 1, T - 1), (R, T - 1), (R - 1, T)]}
                return obs.creutz_ratio(Wd, R, T)
            try:
                val, err = jackknife(W, f, args.nblocks)
                out(f"| chi({R},{T}) | {val:.4f} +- {err:.4f} |")
            except Exception:
                out(f"| chi({R},{T}) | n/a |")
    out("")
    # ---- glueball
    out("## 0++ correlator and effective masses")
    out("")
    out("C(tau) = <O(t) O(t+tau)> - <O>^2 with O(t) the zero-momentum sum of spatial plaquettes; "
        "m_eff from the periodic (cosh) ratio C(tau)/C(tau+1).  Jackknife errors.")
    out("")
    for li, lev in enumerate(levels):
        Ol = O[:, li, :]
        C, Cerr = jackknife(Ol, obs.timeslice_correlator, args.nblocks)
        m, merr = jackknife(Ol, lambda s: obs.effective_mass(obs.timeslice_correlator(s)), args.nblocks)
        out(f"### smearing level {lev}")
        out("")
        out("| tau | C(tau) | m_eff(tau) |")
        out("|---|---|---|")
        for tau in range(min(L0 // 2 + 1, 6)):
            me = f"{m[tau]:.3f} +- {merr[tau]:.3f}" if tau < len(m) and np.isfinite(m[tau]) else "-"
            out(f"| {tau} | {C[tau]:.4g} +- {Cerr[tau]:.2g} | {me} |")
        out("")
    text = "\n".join(lines)
    print(text)
    if args.md:
        with open(args.md, "w") as f:
            f.write(text + "\n")


if __name__ == "__main__":
    sys.exit(main())

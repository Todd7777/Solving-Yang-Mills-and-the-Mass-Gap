#!/usr/bin/env python3
"""Variational (GEVP) analysis of the 0++ glueball from one or more ``glueball.py`` output files.

The time-slice operators O_i(t) at the different APE-smearing levels i form a basis.  With the
connected, symmetrised correlation matrix

    C_ij(tau) = < (1/L0) sum_t O_i(t) O_j(t + tau) > - <O_i><O_j>,

the generalised eigenvalue problem C(tau) v = lambda(tau, tau0) C(tau0) v gives eigenvalues whose
decay isolates the lightest state far better than any single operator.  Effective masses are
obtained from lambda_0(tau)/lambda_0(tau+1) with the periodic (cosh) form.  Errors: block jackknife.

Interpretation (see theory/lattice-foundations.tex, reflection positivity and the spectral
representation): for a single operator, C(tau) = sum_n |<n|O|Omega>|^2 (e^{-E_n tau} + e^{-E_n (L0 - tau)})
is log-convex, so effective masses decrease towards the true gap E_1 from ABOVE.  Monte Carlo
effective masses are therefore estimates of upper bounds on the gap, never lower bounds.

Several files with the same beta, lattice shape and smearing levels can be given (independent
Markov chains, e.g. a run restarted with a new seed); their measurements are pooled.

Example:
    python scripts/analyze_gevp.py ../results/gb_su2_b2.3_8x16.npz ../results/gb_su2_b2.3_8x16_s23.npz \
        --md ../results/gb_su2_b2.3_8x16.md
"""
import argparse
import sys

import numpy as np

from ymlattice import exact
from ymlattice import observables as obs
from ymlattice.stats import autocorrelation_time, binned_mean_err, jackknife


from ymlattice.analysis import corr_matrix, gevp_lambda0, load_runs as load  # noqa: E402


def meff_from(series, L0):
    return obs.effective_mass(series, periodic=True)[: L0 // 2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz", nargs="+")
    ap.add_argument("--tau0", type=int, default=0)
    ap.add_argument("--nblocks", type=int, default=20)
    ap.add_argument("--md", default=None)
    args = ap.parse_args()
    (beta, shape, levels, group), O, plaq, poly, W, sizes = load(args.npz)
    N = 2 if group == "su2" else 3
    L0 = shape[0]
    lines = []
    out = lines.append
    out(f"# 0++ glueball, variational analysis: {group.upper()} beta={beta} lattice {shape[1]}^3 x {L0}")
    out("")
    out(f"- runs pooled: {len(sizes)} ({' + '.join(str(s) for s in sizes)} = {sum(sizes)} measurements)")
    pm, pe = binned_mean_err(plaq, args.nblocks)
    out(f"- plaquette = {pm:.5f} +- {pe:.5f}; tau_int(plaquette) = {autocorrelation_time(plaq):.1f} measurements")
    pr, prerr = binned_mean_err(poly.real, args.nblocks)
    out(f"- Polyakov loop Re<L> = {pr:.4f} +- {prerr:.4f} (zero in the confined phase)")
    aL = exact.a_lambda_lattice(beta, N)
    out(f"- two-loop asymptotic scaling a*Lambda_L = {aL:.4e}")
    out("")

    # Creutz ratios from rows where Wilson loops were measured
    ok = ~np.isnan(W).any(axis=(1, 2))
    Wm = W[ok]
    maxloop = Wm.shape[1]
    out(f"## Creutz ratios ({ok.sum()} configurations with Wilson loops)")
    out("")
    out("| chi(R,R) | value |")
    out("|---|---|")
    chis = {}
    for R in range(2, maxloop + 1):
        def f(s, R=R):
            Wd = {k: np.mean(s[:, k[0] - 1, k[1] - 1]) for k in [(R, R), (R - 1, R - 1), (R, R - 1), (R - 1, R)]}
            return obs.creutz_ratio(Wd, R, R)
        try:
            v, e = jackknife(Wm, f, args.nblocks)
            chis[R] = (v, e)
            out(f"| chi({R},{R}) | {v:.4f} +- {e:.4f} |")
        except Exception:
            out(f"| chi({R},{R}) | n/a |")
    out("")

    # single-operator and GEVP effective masses
    out("## Effective masses")
    out("")
    out(f"GEVP with tau0 = {args.tau0} over smearing levels {list(levels)}; m_eff from the cosh ratio; block jackknife ({args.nblocks} blocks).")
    out("")
    cols = [f"lvl {lv}" for lv in levels] + ["GEVP"]
    out("| tau | " + " | ".join(cols) + " |")
    out("|---|" + "---|" * len(cols))
    res = {}
    for li in range(len(levels)):
        res[li] = jackknife(O[:, li, :], lambda s, li=li: meff_from(corr_matrix(s[:, None, :])[:, 0, 0], L0), args.nblocks)
    res["g"] = jackknife(O, lambda s: meff_from(gevp_lambda0(corr_matrix(s), args.tau0), L0), args.nblocks)
    for tau in range(L0 // 2):
        row = []
        for k in list(range(len(levels))) + ["g"]:
            v, e = res[k]
            row.append(f"{v[tau]:.3f} +- {e[tau]:.3f}" if np.isfinite(v[tau]) and np.isfinite(e[tau]) else "-")
        out(f"| {tau} | " + " | ".join(row) + " |")
    out("")
    v, e = res["g"]
    if L0 // 2 > 1 and np.isfinite(v[1]):
        out(f"GEVP effective mass at tau = 1: a*m = {v[1]:.3f} +- {e[1]:.3f}  (an estimate of an upper bound on the gap); "
            f"in units of the asymptotic-scaling scale: m/Lambda_L = {v[1]/aL:.0f} +- {e[1]/aL:.0f}.")
        if 4 in chis or 3 in chis:
            R = max(k for k in chis)
            s, se = chis[R]
            if s > 0:
                out(f"Ratio to the largest Creutz ratio: a*m / sqrt(chi({R},{R})) = {v[1]/np.sqrt(s):.2f} "
                    f"(chi({R},{R}) over-estimates sigma at small R, so this under-estimates m/sqrt(sigma)).")
    text = "\n".join(lines)
    print(text)
    if args.md:
        with open(args.md, "w") as fh:
            fh.write(text + "\n")


if __name__ == "__main__":
    sys.exit(main())

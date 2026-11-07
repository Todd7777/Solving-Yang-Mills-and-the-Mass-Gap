#!/usr/bin/env python3
"""Scan the average plaquette versus beta and compare with the exact strong- and weak-coupling
expansions.  Writes a CSV with columns beta, plaquette, error, u(beta) [leading strong coupling],
weak-coupling one-loop value, and tau_int.

Example:
    python scripts/plaquette_scan.py --group su2 --L 8 --betas 0.5 1 2 2.3 3 6 12 --out results/su2.csv
"""
import argparse
import csv
import sys
import time

import numpy as np

from ymlattice import Lattice, exact
from ymlattice.stats import autocorrelation_time, binned_mean_err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", default="su2")
    ap.add_argument("--L", type=int, default=8)
    ap.add_argument("--d", type=int, default=4)
    ap.add_argument("--betas", type=float, nargs="+", required=True)
    ap.add_argument("--ntherm", type=int, default=200)
    ap.add_argument("--nmeas", type=int, default=400)
    ap.add_argument("--nor", type=int, default=2, help="overrelaxation sweeps per heat-bath sweep")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    N = 2 if args.group == "su2" else 3
    rows = []
    for beta in args.betas:
        t0 = time.time()
        lat = Lattice((args.L,) * args.d, beta, args.group, seed=args.seed, start="hot")
        lat.thermalize(args.ntherm, n_overrelax=args.nor)
        vals = []
        for _ in range(args.nmeas):
            lat.update(1, args.nor)
            vals.append(lat.plaquette())
        vals = np.array(vals)
        mean, err = binned_mean_err(vals, 20)
        tau = autocorrelation_time(vals)
        u = exact.plaquette_strong_leading(beta, N)
        weak = exact.plaquette_weak(beta, N, args.d)
        rows.append(dict(beta=beta, plaquette=mean, error=err, strong_u=u, weak_1loop=weak, tau_int=tau))
        print(f"beta={beta:6.3f}  <P>={mean:.5f} +- {err:.5f}   u={u:.5f}  weak={weak:.5f}  tau={tau:.2f}  ({time.time()-t0:.0f}s)",
              flush=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote", args.out)


if __name__ == "__main__":
    sys.exit(main())

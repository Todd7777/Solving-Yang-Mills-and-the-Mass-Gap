#!/usr/bin/env python3
"""Scan the average plaquette versus beta and compare with the exact strong- and weak-coupling
expansions.  Writes a CSV with columns beta, plaquette, error, u(beta) [leading strong coupling],
weak-coupling one-loop value, and tau_int.  Rows are appended as each beta finishes; rerunning
the same command skips the betas already present (restart-safe).

Example:
    python scripts/plaquette_scan.py --group su2 --L 8 --betas 0.5 1 2 2.3 3 6 12 --out results/su2.csv
"""
import argparse
import csv
import os
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
    fields = ["beta", "plaquette", "error", "strong_u", "weak_1loop", "tau_int"]
    done = set()
    if os.path.exists(args.out):  # resume: skip betas already in the file
        with open(args.out) as f:
            done = {float(r["beta"]) for r in csv.DictReader(f)}
    else:
        with open(args.out, "w", newline="") as f:
            csv.DictWriter(f, fieldnames=fields).writeheader()
    for beta in args.betas:
        if beta in done:
            print(f"beta={beta:6.3f} already in {args.out}, skipping", flush=True)
            continue
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
        with open(args.out, "a", newline="") as f:
            csv.DictWriter(f, fieldnames=fields).writerow(
                dict(beta=beta, plaquette=mean, error=err, strong_u=u, weak_1loop=weak, tau_int=tau))
        print(f"beta={beta:6.3f}  <P>={mean:.5f} +- {err:.5f}   u={u:.5f}  weak={weak:.5f}  tau={tau:.2f}  ({time.time()-t0:.0f}s)",
              flush=True)
    print("wrote", args.out)


if __name__ == "__main__":
    sys.exit(main())

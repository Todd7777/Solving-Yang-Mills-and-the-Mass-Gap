#!/usr/bin/env python3
"""Measure zero-momentum 0++ time-slice operators (several APE-smearing levels), planar Wilson loops
and the plaquette on a sequence of configurations.  Raw per-configuration data is stored in an .npz
file for later analysis with ``analyze_glueball.py``.  A checkpoint (``<out>.ckpt.npz``: configuration,
RNG state, measurements so far) is written every ``--checkpoint-every`` measurements; rerunning the same
command resumes from it.

Example:
    python scripts/glueball.py --group su2 --beta 2.3 --shape 16 8 8 8 --nmeas 4000 --out results/gb_su2_b2.3.npz
"""
import argparse
import json
import os
import sys
import time

import numpy as np

from ymlattice import Lattice
from ymlattice import observables as obs


def _save(path, arrays, meta):
    tmp = path + ".tmp.npz"
    np.savez(tmp, **arrays, **meta)
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", default="su2")
    ap.add_argument("--beta", type=float, required=True)
    ap.add_argument("--shape", type=int, nargs="+", required=True, help="L0 (time) L1 L2 L3")
    ap.add_argument("--ntherm", type=int, default=500)
    ap.add_argument("--nmeas", type=int, default=2000)
    ap.add_argument("--nor", type=int, default=3)
    ap.add_argument("--nskip", type=int, default=1, help="compound updates between measurements")
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--levels", type=int, nargs="+", default=[0, 2, 4, 8])
    ap.add_argument("--maxloop", type=int, default=4)
    ap.add_argument("--wloop-every", type=int, default=1, help="measure Wilson loops every k-th measurement")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--checkpoint-every", type=int, default=50)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ckpt = args.out + ".ckpt.npz"
    L0 = args.shape[0]
    levels = sorted(set(args.levels))
    meta = dict(levels=np.array(levels), beta=args.beta, shape=np.array(args.shape), group=args.group,
                alpha=args.alpha, nor=args.nor, nskip=args.nskip, seed=args.seed)
    lat = Lattice(tuple(args.shape), args.beta, args.group, seed=args.seed, start="hot")
    if os.path.exists(ckpt):
        c = np.load(ckpt, allow_pickle=True)
        lat.U = c["U"]
        lat.rng.bit_generator.state = json.loads(str(c["rng"]))
        start = int(c["n"])
        O, W, plaq, poly = c["O"], c["W"], c["plaq"], c["poly"]
        print(f"resumed from checkpoint at measurement {start}", flush=True)
    else:
        t0 = time.time()
        lat.thermalize(args.ntherm, n_overrelax=args.nor)
        print(f"thermalised {args.ntherm} updates in {time.time()-t0:.0f}s, <P>={lat.plaquette():.5f}", flush=True)
        start = 0
        O = np.zeros((args.nmeas, len(levels), L0))
        W = np.full((args.nmeas, args.maxloop, args.maxloop), np.nan)
        plaq = np.zeros(args.nmeas)
        poly = np.zeros(args.nmeas, dtype=complex)

    t0 = time.time()
    for n in range(start, args.nmeas):
        for _ in range(args.nskip):
            lat.update(1, args.nor)
        plaq[n] = lat.plaquette()
        poly[n] = obs.polyakov_loop(lat)
        U = lat.U
        done = 0
        for li, lev in enumerate(levels):
            if lev > done:
                U = obs.ape_smear(lat, args.alpha, lev - done, U=U)
                done = lev
            O[n, li] = obs.spatial_plaquette_timeslices(lat, U)
        if n % args.wloop_every == 0:
            for R in range(1, args.maxloop + 1):
                for T in range(1, args.maxloop + 1):
                    W[n, R - 1, T - 1] = obs.wilson_loop(lat, R, T)
        if (n + 1) % args.checkpoint_every == 0 or n + 1 == args.nmeas:
            arrays = dict(O=O[:n + 1], W=W[:n + 1], plaq=plaq[:n + 1], poly=poly[:n + 1])
            _save(args.out, arrays, meta)
            _save(ckpt, dict(O=O, W=W, plaq=plaq, poly=poly, U=lat.U,
                             rng=json.dumps(lat.rng.bit_generator.state), n=n + 1), {})
        if (n + 1) % 100 == 0:
            el = time.time() - t0
            print(f"{n+1}/{args.nmeas}  <P>={plaq[:n+1].mean():.5f}  {el:.0f}s elapsed, "
                  f"~{el/(n+1-start)*(args.nmeas-n-1):.0f}s left", flush=True)
    print("wrote", args.out)


if __name__ == "__main__":
    sys.exit(main())

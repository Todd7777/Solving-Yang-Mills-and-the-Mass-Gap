#!/usr/bin/env python3
"""Rerun the deficiency search with other choices of the deficient link (largest; random, three seeds); the completeness
argument holds for any choice, so the sets must coincide.  Usage: python3 variant_orders.py D L"""
import sys
from lattice16b import Lattice
from fronts16b import FrontCache
from enum_junctions import enumerate_junctions

D, L = int(sys.argv[1]), int(sys.argv[2])
lat = Lattice(D, L)
fc = FrontCache(lat, [4, 6, 7])
q = lat.plaq(lat.origin(), 0, 1)
for n, sizes in ((3, (4, 4)), (2, (4, 6)), (2, (4, 7))):
    for lam in (0, 1):
        ref = set(enumerate_junctions(lat, fc, q, lam, n, sizes, order='min'))
        res = [len(ref)]
        same = True
        for order, seed in (('max', 0), ('random', 1), ('random', 2), ('random', 3)):
            s = set(enumerate_junctions(lat, fc, q, lam, n, sizes, order=order, seed=seed))
            res.append(len(s))
            same &= (s == ref)
        print(f"D={D} L={L} n={n} sizes={sizes} lam={lam}: counts (min, max, rnd1, rnd2, rnd3) = {res}; identical: {same}", flush=True)

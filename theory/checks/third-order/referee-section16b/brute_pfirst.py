#!/usr/bin/env python3
"""A second, structurally different enumeration (plaquettes first), to cross-check enum_junctions.py.

Loop over all sets P of n plaquettes of the torus that contain a plaquette sharing a vertex with q (necessary: p is touched
by a marker).  Given P, the links of T := {e in E(P) : exactly one plaquette of P contains e, e not in dq} u (dq \\ E(P)) need a
front other than p (condition (N)), so the other front-pairs cover T; every front-pair touches V(P) (condition (T)).  The
fronts are chosen to cover the smallest uncovered link of T, or among all fronts touching V(P) once T is covered.  Every
candidate is then checked with is_junction (all conditions).  No result of the paper is used.
Usage: python3 brute_pfirst.py D L n sizes lam
"""
import sys
import time
import itertools
from collections import Counter
from lattice16b import Lattice
from fronts16b import FrontCache
from enum_junctions import is_junction, enumerate_junctions


def brute(lat, fc, q, lam, n, sizes):
    p = (frozenset(lat.plaq_edges(q)), lam)
    dq = p[0]
    allplaq = range(lat.nv * lat.NP)
    nearq = set()
    for v in lat.plaq_verts(q):
        nearq.update(lat.plaqs_at_vertex(v))
    Psets = set()
    for c1 in nearq:
        for rest in itertools.combinations([c for c in allplaq if c != c1], n - 1):
            Psets.add(frozenset((c1,) + rest))
    found = set()
    sizes = tuple(sizes)
    for P in Psets:
        pi = Counter()
        for c in P:
            for e in lat.plaq_edges(c):
                pi[e] += 1
        T = {e for e, m in pi.items() if m == 1 and e not in dq} | {e for e in dq if e not in pi}
        if len(T) > sum(sizes):
            continue
        VP = lat.verts_of_plaqs(P)

        def touching(k):
            s = set()
            for v in VP:
                s.update(fc.through(v, k))
            return s

        def place(chosen, remaining_sizes, Tleft):
            if not remaining_sizes:
                if not Tleft:
                    Phi = frozenset([p] + chosen)
                    if len(Phi) == 1 + len(sizes) and is_junction(lat, P, Phi, p=p, nfp=1 + len(sizes)):
                        found.add((P, Phi))
                return
            if Tleft:
                e = min(Tleft)
                for idx, k in enumerate(remaining_sizes):
                    if k in remaining_sizes[:idx]:
                        continue
                    rs = remaining_sizes[:idx] + remaining_sizes[idx + 1:]
                    for F in fc.containing_edge(e, k):
                        if not (lat.verts_of_edges(F) & VP):
                            continue
                        for mu in (0, 1):
                            place(chosen + [(F, mu)], rs, Tleft - F)
            else:
                k = remaining_sizes[0]
                rs = remaining_sizes[1:]
                for F in touching(k):
                    for mu in (0, 1):
                        place(chosen + [(F, mu)], rs, Tleft)

        place([], sizes, T)
    return found


def main():
    D, L, n = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    sizes = tuple(int(x) for x in sys.argv[4].split(','))
    lam = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    lat = Lattice(D, L)
    fc = FrontCache(lat, sorted(set(sizes)))
    q = lat.plaq(lat.origin(), 0, 1)
    t0 = time.time()
    b = brute(lat, fc, q, lam, n, sizes)
    t1 = time.time()
    s = set(enumerate_junctions(lat, fc, q, lam, n, sizes))
    print(f"D={D} L={L} n={n} sizes={sizes} lam={lam}: plaquettes-first {len(b)} ({t1 - t0:.0f}s), "
          f"deficiency search {len(s)}, identical sets: {b == s}", flush=True)


if __name__ == "__main__":
    main()

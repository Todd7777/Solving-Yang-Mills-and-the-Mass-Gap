#!/usr/bin/env python3
"""Plaquettes-first brute force for case (i) of Lemma mg:lem:threelow (|P| = 3, other front-pairs with 4 and 4 edges),
made feasible in d = 4 by one more necessary condition.

As in brute_pfirst.py: P = {c1, c2, c3} with c1 sharing a vertex with q; T := links in exactly one plaquette of P and not
in dq, together with the links of dq outside E(P); the two other fronts (4 edges each) cover T, so |T| <= 8.
Pruning (justified): given c1, c2, let T2 be T for {c1, c2}.  A third plaquette c3 sharing no link with
E(c1) u E(c2) u dq adds its four links to T (and removes none), so if |T2| > 4 then c3 contains a link of
E(c1) u E(c2) u dq.  All surviving P are completed as in brute_pfirst.py and checked with is_junction.
Usage: python3 brute_pfirst_i.py D L lam
"""
import sys
import time
from collections import Counter
from lattice16b import Lattice
from fronts16b import FrontCache
from enum_junctions import is_junction, enumerate_junctions


def main():
    D, L, lam = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    lat = Lattice(D, L)
    fc = FrontCache(lat, [4])
    q = lat.plaq(lat.origin(), 0, 1)
    p = (frozenset(lat.plaq_edges(q)), lam)
    dq = p[0]
    allplaq = list(range(lat.nv * lat.NP))
    PE = [frozenset(lat.plaq_edges(c)) for c in allplaq]
    nearq = set()
    for v in lat.plaq_verts(q):
        nearq.update(lat.plaqs_at_vertex(v))

    def Tof(P):
        pi = Counter()
        for c in P:
            for e in PE[c]:
                pi[e] += 1
        return {e for e, m in pi.items() if m == 1 and e not in dq} | {e for e in dq if e not in pi}

    t0 = time.time()
    Psets = set()
    for c1 in nearq:
        for c2 in allplaq:
            if c2 == c1:
                continue
            T2 = Tof((c1, c2))
            if len(T2) > 4:
                S12 = PE[c1] | PE[c2] | dq
                c3s = set()
                for e in S12:
                    c3s.update(lat.plaqs_at_edge(e))
            else:
                c3s = allplaq
            for c3 in c3s:
                if c3 in (c1, c2):
                    continue
                P = frozenset((c1, c2, c3))
                if P in Psets:
                    continue
                if len(Tof(P)) <= 8:
                    Psets.add(P)
    t1 = time.time()
    found = set()
    for P in Psets:
        T = Tof(P)
        VP = lat.verts_of_plaqs(P)

        def place(chosen, left, Tleft):
            if left == 0:
                if not Tleft:
                    Phi = frozenset([p] + chosen)
                    if len(Phi) == 3 and is_junction(lat, P, Phi, p=p, nfp=3):
                        found.add((P, Phi))
                return
            if Tleft:
                cands = [F for F in fc.containing_edge(min(Tleft), 4) if lat.verts_of_edges(F) & VP]
            else:
                s = set()
                for v in VP:
                    s.update(fc.through(v, 4))
                cands = list(s)
            for F in cands:
                for mu in (0, 1):
                    place(chosen + [(F, mu)], left - 1, Tleft - F)

        place([], 2, T)
    t2 = time.time()
    s = set(enumerate_junctions(lat, fc, q, lam, 3, (4, 4)))
    print(f"D={D} L={L} lam={lam}: {len(Psets)} plaquette sets with |T|<=8 ({t1 - t0:.0f}s); plaquettes-first {len(found)} "
          f"({t2 - t1:.0f}s), deficiency search {len(s)}, identical sets: {found == s}", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exhaustive enumeration of candidate junctions containing a given four-edge front-pair p = (dq, lam).

A candidate junction at slice s (layers 0 = s-1 and 1 = s) is a pair (P, Phi): P a finite set of plaquettes of the slice
(the plaquettes of the clusters) and Phi a set of front-pairs (phi, mu), mu in {0,1}, such that (Definitions 8.1, 9.1 and
Lemma 9.2 of the paper, nothing else):
  (F)  every phi is a front: non-empty, connected, every vertex of degree >= 2;
  (L)  distinct fronts of the same layer are vertex-disjoint (distinct components of T_mu);
  (N)  n_J(e) := [e in E_0] + [e in E_1] + #{c in P : e in c} >= 2 for every link e of E(J) (no support in J);
  (T)  every front-pair shares a vertex with E(P) (it is touched by a marker; there is no support);
  (C)  J is connected: the graph on P u Phi, adjacency = a common vertex, is connected (equivalent to the
       connectivity of H_s, clusters being link-connected).
The clusters are the link-components of P; they impose nothing further.

Search: start from p.  While some link of E(current) has n(e) = 1, branch over all items (a plaquette, or a front of a
still allowed size in a layer where (L) permits it) that contain the smallest such link: every completion contains one.
When no link is deficient and items remain to be placed, branch over all items sharing a vertex with the current ones:
by (C) every completion contains one.  Complete; duplicates are removed at the end.  Every result is re-checked by the
independent function `is_junction`.

Usage:  python3 enum_junctions.py D L n sizes lam      e.g.  python3 enum_junctions.py 3 8 2 4,6 0
        (sizes = sizes of the front-pairs other than p, comma-separated)
"""
import sys
import time
from collections import Counter
from lattice16b import Lattice, is_front
from fronts16b import FrontCache


def is_junction(lat, P, Phi, p=None, nfp=None):
    """Independent check of (F), (L), (N), (T), (C) for a candidate (P, Phi)."""
    P = list(P)
    Phi = list(Phi)
    if p is not None and p not in Phi:
        return False
    if nfp is not None and len(Phi) != nfp:
        return False
    if len(set(P)) != len(P) or len(set(Phi)) != len(Phi):
        return False
    for (phi, mu) in Phi:
        if mu not in (0, 1) or not is_front(lat, phi):
            return False
    # (L)
    for a in range(len(Phi)):
        for b in range(a + 1, len(Phi)):
            (f1, m1), (f2, m2) = Phi[a], Phi[b]
            if m1 == m2 and (lat.verts_of_edges(f1) & lat.verts_of_edges(f2)):
                return False
    # (N)
    n = Counter()
    for (phi, mu) in Phi:
        for e in phi:
            n[e] += 1
    for c in P:
        for e in lat.plaq_edges(c):
            n[e] += 1
    if any(v < 2 for v in n.values()):
        return False
    # (T)
    VP = lat.verts_of_plaqs(P)
    for (phi, mu) in Phi:
        if not (lat.verts_of_edges(phi) & VP):
            return False
    # (C)
    items = [set(lat.plaq_verts(c)) for c in P] + [lat.verts_of_edges(phi) for (phi, mu) in Phi]
    seen = {0}
    stack = [0]
    while stack:
        i = stack.pop()
        for j in range(len(items)):
            if j not in seen and items[i] & items[j]:
                seen.add(j)
                stack.append(j)
    return len(seen) == len(items)


def enumerate_junctions(lat, fc, q, lam, n, sizes, order='min', seed=0):
    """All candidate junctions (P, Phi) with |P| = n and Phi = {p} + front-pairs of the given sizes, p = (dq, lam)."""
    p = (frozenset(lat.plaq_edges(q)), lam)
    P = []
    Pset = set()
    Phi = []
    cnt = Counter()          # n(e)
    lay_verts = [Counter(), Counter()]
    allv = Counter()
    rem = {'plaq': n, 'sizes': Counter(sizes)}
    found = set()
    vcache = {}
    import random
    rng = random.Random(seed)

    def V(F):
        r = vcache.get(F)
        if r is None:
            r = lat.verts_of_edges(F)
            vcache[F] = r
        return r

    def add_front(F, mu):
        Phi.append((F, mu))
        for e in F:
            cnt[e] += 1
        for v in V(F):
            lay_verts[mu][v] += 1
            allv[v] += 1

    def rem_front():
        F, mu = Phi.pop()
        for e in F:
            cnt[e] -= 1
            if cnt[e] == 0:
                del cnt[e]
        for v in V(F):
            lay_verts[mu][v] -= 1
            if lay_verts[mu][v] == 0:
                del lay_verts[mu][v]
            allv[v] -= 1
            if allv[v] == 0:
                del allv[v]

    def add_plaq(c):
        P.append(c)
        Pset.add(c)
        for e in lat.plaq_edges(c):
            cnt[e] += 1
        for v in lat.plaq_verts(c):
            allv[v] += 1

    def rem_plaq():
        c = P.pop()
        Pset.discard(c)
        for e in lat.plaq_edges(c):
            cnt[e] -= 1
            if cnt[e] == 0:
                del cnt[e]
        for v in lat.plaq_verts(c):
            allv[v] -= 1
            if allv[v] == 0:
                del allv[v]

    def front_ok(F, mu):
        if (F, mu) in Phi:
            return False
        lv = lay_verts[mu]
        return not any(v in lv for v in V(F))

    def rec():
        defic = [e for e, c in cnt.items() if c == 1]
        cap = 4 * rem['plaq'] + sum(k * m for k, m in rem['sizes'].items())
        if len(defic) > cap:
            return
        ks = [k for k, m in rem['sizes'].items() if m > 0]
        if defic:
            if order == 'min':
                e = min(defic)
            elif order == 'max':
                e = max(defic)
            else:
                e = rng.choice(sorted(defic))
            if rem['plaq'] > 0:
                for c in lat.plaqs_at_edge(e):
                    if c not in Pset:
                        rem['plaq'] -= 1
                        add_plaq(c)
                        rec()
                        rem_plaq()
                        rem['plaq'] += 1
            for k in ks:
                for F in fc.containing_edge(e, k):
                    for mu in (0, 1):
                        if front_ok(F, mu):
                            rem['sizes'][k] -= 1
                            add_front(F, mu)
                            rec()
                            rem_front()
                            rem['sizes'][k] += 1
            return
        if rem['plaq'] == 0 and not ks:
            found.add((frozenset(P), frozenset(Phi)))
            return
        verts = list(allv)
        if rem['plaq'] > 0:
            cands = set()
            for v in verts:
                cands.update(lat.plaqs_at_vertex(v))
            for c in sorted(cands - Pset):
                rem['plaq'] -= 1
                add_plaq(c)
                rec()
                rem_plaq()
                rem['plaq'] += 1
        for k in ks:
            cf = set()
            for v in verts:
                cf.update(fc.through(v, k))
            for F in cf:
                for mu in (0, 1):
                    if front_ok(F, mu):
                        rem['sizes'][k] -= 1
                        add_front(F, mu)
                        rec()
                        rem_front()
                        rem['sizes'][k] += 1

    add_front(p[0], lam)
    rec()
    out = []
    for (Pf, Phif) in found:
        assert is_junction(lat, Pf, Phif, p=p, nfp=1 + len(sizes)), "search produced an invalid junction"
        out.append((Pf, Phif))
    return out


def main():
    D, L, n = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    sizes = tuple(int(x) for x in sys.argv[4].split(','))
    lam = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    lat = Lattice(D, L)
    fc = FrontCache(lat, sorted(set(sizes)))
    q = lat.plaq(lat.origin(), 0, 1)
    t0 = time.time()
    res = enumerate_junctions(lat, fc, q, lam, n, sizes)
    print(f"D={D} L={L} n={n} sizes={sizes} lam={lam}: {len(res)} junctions  ({time.time() - t0:.1f}s)")


if __name__ == "__main__":
    main()

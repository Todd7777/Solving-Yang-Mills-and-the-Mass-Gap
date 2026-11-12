#!/usr/bin/env python3
"""Fronts (connected link sets in which every vertex has degree >= 2) through a vertex, from the definition.

Search: grow a link set from a root link; while some vertex has degree 1, branch over the links at the smallest
such vertex (every front containing the current set contains one of them); when no vertex has degree 1 and the
size is below k, branch over all links touching the current vertex set (a front is connected).  Complete; the
duplicates are removed by storing frozensets.  No structural fact about fronts (cycles, thetas, ...) is used.

Usage: python3 fronts16b.py D L KMAX     prints N^v_k for k = 4..KMAX on the torus (Z/L)^D.
"""
import sys
import time
from lattice16b import Lattice, is_front


def fronts_with_root(lat, e0, k, out):
    S = {e0}
    deg = {}
    a, b = lat.ends(e0)
    deg[a] = 1
    deg[b] = 1

    def add(e):
        S.add(e)
        x, y = lat.ends(e)
        deg[x] = deg.get(x, 0) + 1
        deg[y] = deg.get(y, 0) + 1

    def rem(e):
        S.discard(e)
        x, y = lat.ends(e)
        for z in (x, y):
            deg[z] -= 1
            if deg[z] == 0:
                del deg[z]

    def rec():
        ones = [v for v, d in deg.items() if d == 1]
        room = k - len(S)
        if (len(ones) + 1) // 2 > room:
            return
        if room == 0:
            if not ones:
                out.add(frozenset(S))
            return
        if ones:
            v = min(ones)
            cands = [e for e in lat.edges_at(v) if e not in S]
        else:
            cands = sorted({e for v in deg for e in lat.edges_at(v) if e not in S})
        for e in cands:
            add(e)
            rec()
            rem(e)

    rec()


def fronts_through_vertex(lat, v0, k):
    out = set()
    for e0 in lat.edges_at(v0):
        fronts_with_root(lat, e0, k, out)
    # sanity: every set found is a front of size k through v0
    for F in out:
        assert len(F) == k and is_front(lat, F) and v0 in lat.verts_of_edges(F)
    return out


class FrontCache:
    """Fronts of size k through any vertex, by translating those through the origin."""

    def __init__(self, lat, sizes):
        self.lat = lat
        o = lat.origin()
        self.base = {}
        for k in sizes:
            fr = fronts_through_vertex(lat, o, k)
            # store as lists of (vertex offset, direction)
            enc = []
            for F in fr:
                enc.append(tuple((lat.vdiff(o, e // lat.D), e % lat.D) for e in F))
            self.base[k] = enc
        self.cache = {}

    def through(self, v, k):
        key = (v, k)
        r = self.cache.get(key)
        if r is None:
            lat = self.lat
            D = lat.D
            r = []
            for enc in self.base[k]:
                r.append(frozenset(lat.vadd(v, t) * D + i for (t, i) in enc))
            self.cache[key] = r
        return r

    def containing_edge(self, e, k):
        a, _ = self.lat.ends(e)
        return [F for F in self.through(a, k) if e in F]


if __name__ == "__main__":
    D, L, KMAX = (int(x) for x in sys.argv[1:4])
    lat = Lattice(D, L)
    o = lat.origin()
    for k in range(4, KMAX + 1):
        t0 = time.time()
        fr = fronts_through_vertex(lat, o, k)
        print(f"D={D} L={L} k={k}: N^v_k = {len(fr)}   ({time.time() - t0:.1f}s)", flush=True)

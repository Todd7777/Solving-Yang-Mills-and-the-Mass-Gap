#!/usr/bin/env python3
"""Periodic cubic lattice (Z/L)^D with integer-coded vertices, links and plaquettes.

Written for the referee report on Section 16 (Lemmas mg:lem:threelow, mg:lem:jsums5) of
theory/mass-gap-strong-coupling.tex.  Independent of the authors' scripts.

Coding:
  vertex v      : index into coords (tuple of D ints mod L)
  link   e      : v*D + i   joins v and v+e_i
  plaquette p   : v*NP + k  with k the index of the pair (i<j); corners v, v+e_i, v+e_j, v+e_i+e_j
A large L (>= 2*extent of the configurations) makes the torus equivalent to Z^D for local questions;
L = 8 tests the paper's hypothesis 'all spatial sides >= 8'.
"""
import itertools


class Lattice:
    def __init__(self, D, L):
        assert L >= 3
        self.D, self.L = D, L
        self.coords = list(itertools.product(range(L), repeat=D))
        self.vid = {c: i for i, c in enumerate(self.coords)}
        self.nv = len(self.coords)
        D_ = D
        self.sh = []   # v + e_i
        self.shm = []  # v - e_i
        for c in self.coords:
            self.sh.append([self.vid[tuple((c[k] + (k == i)) % L for k in range(D_))] for i in range(D_)])
            self.shm.append([self.vid[tuple((c[k] - (k == i)) % L for k in range(D_))] for i in range(D_)])
        self.pairs = [(i, j) for i in range(D) for j in range(i + 1, D)]
        self.NP = len(self.pairs)
        self.pidx = {pr: k for k, pr in enumerate(self.pairs)}
        self._edges_at = {}
        self._plaqs_at_edge = {}
        self._plaqs_at_vertex = {}

    # ---- links -------------------------------------------------------------
    def link(self, v, i):
        return v * self.D + i

    def ends(self, e):
        v, i = divmod(e, self.D)
        return (v, self.sh[v][i])

    def direction(self, e):
        return e % self.D

    def edges_at(self, v):
        r = self._edges_at.get(v)
        if r is None:
            D = self.D
            r = tuple([v * D + i for i in range(D)] + [self.shm[v][i] * D + i for i in range(D)])
            self._edges_at[v] = r
        return r

    # ---- plaquettes --------------------------------------------------------
    def plaq(self, v, i, j):
        if i > j:
            i, j = j, i
        return v * self.NP + self.pidx[(i, j)]

    def plaq_dirs(self, p):
        return self.pairs[p % self.NP]

    def plaq_corner(self, p):
        return p // self.NP

    def plaq_edges(self, p):
        v, k = divmod(p, self.NP)
        i, j = self.pairs[k]
        D = self.D
        vi = self.sh[v][i]
        vj = self.sh[v][j]
        return (v * D + i, vi * D + j, vj * D + i, v * D + j)

    def plaq_verts(self, p):
        v, k = divmod(p, self.NP)
        i, j = self.pairs[k]
        vi = self.sh[v][i]
        return (v, vi, self.sh[vi][j], self.sh[v][j])

    def plaqs_at_edge(self, e):
        r = self._plaqs_at_edge.get(e)
        if r is None:
            v, i = divmod(e, self.D)
            out = []
            for j in range(self.D):
                if j == i:
                    continue
                out.append(self.plaq(v, i, j))
                out.append(self.plaq(self.shm[v][j], i, j))
            r = tuple(out)
            self._plaqs_at_edge[e] = r
        return r

    def plaqs_at_vertex(self, v):
        r = self._plaqs_at_vertex.get(v)
        if r is None:
            out = []
            for (i, j) in self.pairs:
                vi = self.shm[v][i]
                vj = self.shm[v][j]
                vij = self.shm[vi][j]
                out += [self.plaq(v, i, j), self.plaq(vi, i, j), self.plaq(vj, i, j), self.plaq(vij, i, j)]
            r = tuple(out)
            self._plaqs_at_vertex[v] = r
        return r

    # ---- translations ------------------------------------------------------
    def vadd(self, v, t):
        c = self.coords[v]
        L = self.L
        return self.vid[tuple((c[k] + t[k]) % L for k in range(self.D))]

    def vdiff(self, v, w):
        """w - v as a vector with entries in (-L/2, L/2]."""
        cv, cw = self.coords[v], self.coords[w]
        L = self.L
        out = []
        for a, b in zip(cv, cw):
            d = (b - a) % L
            if d > L // 2:
                d -= L
            out.append(d)
        return tuple(out)

    def tlink(self, e, t):
        v, i = divmod(e, self.D)
        return self.vadd(v, t) * self.D + i

    def tplaq(self, p, t):
        v, k = divmod(p, self.NP)
        return self.vadd(v, t) * self.NP + k

    # ---- helpers -----------------------------------------------------------
    def verts_of_edges(self, E):
        s = set()
        for e in E:
            a, b = self.ends(e)
            s.add(a)
            s.add(b)
        return s

    def verts_of_plaqs(self, P):
        s = set()
        for p in P:
            s.update(self.plaq_verts(p))
        return s

    def origin(self):
        return self.vid[tuple([0] * self.D)]


def is_front(lat, E):
    """Definition (Section 8): non-empty, connected (edges linked through common endpoints), every vertex of degree >= 2."""
    if not E:
        return False
    deg = {}
    adj = {}
    for e in E:
        a, b = lat.ends(e)
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    if min(deg.values()) < 2:
        return False
    start = next(iter(adj))
    seen = {start}
    stack = [start]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == len(adj)

# Written by the independent referee of Section 13 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs. Code written from the definitions alone (no lemma of the paper is used).
# Independent check of Lemma 13.1(i) and (iv): all fronts have four edges, so in Z^n they are plaquette boundaries, and the condition
# n(e) >= 2 (J4) says that the multiset M = {fronts} + {plaquettes of P} of plaquette boundaries covers every link of its union at least twice
# (same-layer fronts are disjoint, so [e in A] + [e in B] = number of fronts containing e).  We enumerate ALL such multisets M containing q
# (q = plaquette at the origin in the plane (0,1)) of a given size, connected through common vertices, by branching on a link covered only
# once (some other element must contain it) or, when M is already closed, on a plaquette sharing a vertex with M.  Then we assign roles
# (p = (bd q, -1), the other front-pairs with their layers, the plaquettes of P, which must be distinct) and check J1-J4.
#   (i):  |M| = 6: two front-pairs p, p' and |P| = 4.  Lemma: P = faces of a unit cube minus q, q'; Z_{4,4} = 12(d-3).
#   (iv): |M| = 5: three four-edge front-pairs and |P| = 2.  Lemma: impossible.
#   also |M| = 4, two front-pairs and |P| = 2 (shifts and jumps, 4(P_v-1)-K_S), and |M| = 5, two front-pairs and |P| = 3 (Lemma 11.2: none),
#   as a test of the code.
import sys, time, itertools
from collections import Counter, defaultdict
from lattice import *

def key(p): return tuple(sorted(p))

def closed_multisets(q, size, n):
    res = set(); visited = set()
    def rec(M):
        if M in visited: return
        visited.add(M)
        cnt = Counter(e for p in M for e in p)
        defic = sorted(e for e, c in cnt.items() if c == 1)
        if len(M) == size:
            if not defic: res.add(M)
            return
        if len(defic) > 4 * (size - len(M)): return
        if defic:
            for r in plaqs_containing_edge(defic[0], n):
                rec(tuple(sorted(M + (key(r),))))
        else:
            V = verts(cnt.keys())
            cands = set()
            for v in V: cands |= plaqs_at_vertex(v, n)
            for r in cands:
                rec(tuple(sorted(M + (key(r),))))
    rec((key(q),))
    return res

def junction_ok(P, FP):
    """P: list of plaquettes (frozensets of edges), FP: list of (front edge set, layer)."""
    if len(set(P)) != len(P): return False
    for (f1, l1), (f2, l2) in itertools.combinations(FP, 2):
        if l1 == l2 and (f1 == f2 or (verts(f1) & verts(f2))): return False   # J1
    cnt = Counter()
    for f, l in FP:
        for e in f: cnt[e] += 1
    for p in P:
        for e in p: cnt[e] += 1
    if min(cnt.values()) < 2: return False                                    # J4
    # clusters: link-components of P
    par = list(range(len(P)))
    def fd(i):
        while par[i] != i: par[i] = par[par[i]]; i = par[i]
        return i
    for i, j in itertools.combinations(range(len(P)), 2):
        if P[i] & P[j]: par[fd(i)] = fd(j)
    comps = defaultdict(set)
    for i in range(len(P)): comps[fd(i)] |= P[i]
    nodes = [verts(es) for es in comps.values()]
    ncl = len(nodes)
    VP = frozenset().union(*nodes) if nodes else frozenset()
    for f, l in FP:
        if not (verts(f) & VP): return False                                  # J2
    nodes += [verts(f) for f, l in FP]
    lay = [None] * ncl + [l for f, l in FP]
    seen = {0}; st = [0]                                                      # J3
    while st:
        a = st.pop()
        for b in range(len(nodes)):
            if b in seen: continue
            if lay[a] is not None and lay[b] is not None and lay[a] == lay[b]: continue
            if nodes[a] & nodes[b]: seen.add(b); st.append(b)
    return len(seen) == len(nodes)

def classify(q, r, n):
    if r == q: return 'q\'=q'
    if len(q & r) == 1: return 'adjacent (share a link)'
    if not (verts(q) & verts(r)):
        # translate by a unit vector perpendicular to the plane?
        vq = sorted(verts(q)); vr = sorted(verts(r))
        d = tuple(b - a for a, b in zip(vq[0], vr[0]))
        if sum(abs(x) for x in d) == 1 and frozenset(E(add(a, d), add(b, d)) for a, b in q) == r: return 'opposite (unit translate)'
        return 'disjoint, other'
    return 'share a vertex only'

def is_cube_minus(P, q, r, n):
    Y = list(P) + [q, r]
    V = verts(frozenset().union(*Y))
    if len(V) != 8: return False
    span = sorted(max(v[k] for v in V) - min(v[k] for v in V) for k in range(n))
    if span != [0] * (n - 3) + [1, 1, 1]: return False
    return len(set(Y)) == 6

def roles(M, q, nfronts, n):
    """all role assignments: p = one copy of q (layer -1); nfronts-1 further fronts with layers; the rest = P"""
    elems = [frozenset(m) for m in M]
    i0 = elems.index(q)
    rest = elems[:i0] + elems[i0 + 1:]
    out = set()
    for idx in itertools.combinations(range(len(rest)), nfronts - 1):
        fr = [rest[i] for i in idx]; P = [rest[i] for i in range(len(rest)) if i not in idx]
        for lays in itertools.product((-1, 1), repeat=nfronts - 1):
            FP = [(q, -1)] + list(zip(fr, lays))
            if len(set(FP)) != len(FP): continue
            if junction_ok(P, FP):
                out.add((frozenset(FP[1:]), frozenset(P)))
    return out

def run(n):
    o = tuple([0] * n); q = plaq(o, 0, 1, n)
    Pv = 2 * n * (n - 1); KS = 4 * (2 * n - 3); d = n + 1
    print(f"n={n} (d={d}): P_v={Pv}, K_S={KS}")
    for size, nfr, label, expect in [(4, 2, 'two front-pairs, |P|=2 (shifts/jumps)', f'4(P_v-1)-K_S = {4*(Pv-1)-KS}'),
                                     (5, 2, 'two front-pairs, |P|=3 (Lemma 11.2)', '0'),
                                     (6, 2, 'two front-pairs, |P|=4 (Lemma 13.1(i))', f'12(d-3) = {12*(d-3)}'),
                                     (5, 3, 'three front-pairs, |P|=2 (Lemma 13.1(iv))', '0')]:
        t = time.time()
        Ms = closed_multisets(q, size, n)
        sols = set()
        for M in Ms:
            if frozenset(q) not in [frozenset(m) for m in M]: continue
            sols |= roles(M, q, nfr, n)
        tal = Counter()
        for FPs, P in sols:
            if nfr == 2:
                (r, lay), = list(FPs)
                c = classify(q, r, n)
                extra = ''
                if size == 6: extra = ', cube minus q,q\'' if is_cube_minus(P, q, r, n) else ', NOT a cube'
                tal[('same layer' if lay == -1 else 'other layer', c + extra)] += 1
            else:
                tal['found'] += 1
        print(f"  {label}: {len(Ms)} closed multisets; {len(sols)} junction configurations (expected {expect})  ({time.time()-t:.1f}s)")
        for k, v in sorted(tal.items()): print(f"      {v:4d}  {k}")
        sys.stdout.flush()

if __name__ == "__main__":
    for n in [int(x) for x in sys.argv[1:]]:
        run(n)

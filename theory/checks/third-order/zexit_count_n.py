# Section 16 of Part IV (Remark mg:rem:nextorder): Z_{n,k}, the junctions with two front-pairs p = (bd q, lambda), p' with |p'| = k >= 5, and n
# plaquettes.  Usage: python3 zexit_count_n.py D NMAX
"""By Lemma mg:lem:twofront(i) (A triangle B is contained in E(P), A and B the unions of the fronts of the two layers), the front of p' is contained in
G = bd q u E(P), and every plaquette of P is joined to V(q) by a chain of plaquettes of P sharing vertices.  The sets P with this property are generated
by adding, one at a time, a plaquette sharing a vertex with q or with a plaquette already chosen (each set once); for each, the fronts in G with at least
five edges are found by a decision search restricted to G, and the conditions of Lemma mg:lem:junction are checked as in zexit_count.py.  p' may lie in
either layer: in the layer of p only if its front is vertex-disjoint from bd q (fronts of one layer are vertex-disjoint), and then the coverage
condition is the same in both layers, so such a front gives two pairs (J, p'); a front meeting V(q) gives one (the other layer).  (The sixth referee
found that the first version counted the vertex-disjoint fronts once; for n <= 3 there are none, so Z_(3,k) is unchanged; for n = 4 in d=4 there are
8, all of six edges.)"""
import sys, itertools
from collections import Counter
from fronts_decision import make
D = int(sys.argv[1]); NMAX = int(sys.argv[2])
E, add, sub, inc, ends = make(D)
O = (0,) * D
def plinks(p):
    x, i, j = p
    return frozenset(((x, i), (add(x, E[i]), j), (add(x, E[j]), i), (x, j)))
def pverts(p):
    x, i, j = p
    return frozenset((x, add(x, E[i]), add(x, E[j]), add(add(x, E[i]), E[j])))
def plaqs_at(v):
    out = []
    for i, j in itertools.combinations(range(D), 2):
        for si in (0, 1):
            for sj in (0, 1):
                y = sub(v, E[i]) if si else v
                y = sub(y, E[j]) if sj else y
                out.append((y, i, j))
    return out
q = (O, 0, 1); phi = plinks(q); Vq = pverts(q)
def fronts_in(G):
    """all fronts (connected, all degrees >= 2) contained in the edge set G, by a decision search rooted at their smallest edge"""
    G = sorted(G); Gs = set(G); out = []
    def incG(v): return [l for l in inc(v) if l in Gs]
    for idx, e0 in enumerate(G):
        allowed = set(G[idx:])
        def rec(S, X, deg, decided):
            cand = [w for w, c in deg.items() if c >= 1 and w not in decided]
            if not cand:
                out.append(frozenset(S)); return
            v = min(cand, key=lambda w: (0 if deg[w] == 1 else 1, w))
            U = [l for l in incG(v) if l in allowed and l not in S and l not in X]
            for r in range(0, len(U) + 1):
                if deg[v] + r < 2: continue
                for T in itertools.combinations(U, r):
                    deg2 = Counter(deg)
                    for l in T:
                        for w in ends(l): deg2[w] += 1
                    rec(S | set(T), X | (set(U) - set(T)), deg2, decided | {v})
        deg0 = Counter(ends(e0))
        rec({e0}, set(), deg0, frozenset())
    return out
def sets_P(nmax):
    """all sets of at most nmax plaquettes in which every plaquette is joined to V(q) by a chain sharing vertices"""
    start = {c for v in Vq for c in plaqs_at(v)}
    seen = set(); out = []
    frontier = [frozenset([c]) for c in start]
    for S in frontier: seen.add(S)
    while frontier:
        nxt = []
        for S in frontier:
            out.append(S)
            if len(S) == nmax: continue
            V = Vq | frozenset(v for c in S for v in pverts(c))
            for v in V:
                for c in plaqs_at(v):
                    if c in S: continue
                    T = S | {c}
                    if T not in seen: seen.add(T); nxt.append(T)
        frontier = nxt
    return out
found = Counter()
for P in sets_P(NMAX):
    if len(P) < 3: continue
    L = phi | frozenset(l for c in P for l in plinks(c))
    pi = Counter(l for c in P for l in plinks(c))
    VP = frozenset(v for c in P for v in pverts(c))
    Pl = list(P)
    comp = list(range(len(Pl)))
    def fd(a):
        while comp[a] != a: comp[a] = comp[comp[a]]; a = comp[a]
        return a
    for a, b in itertools.combinations(range(len(Pl)), 2):
        if plinks(Pl[a]) & plinks(Pl[b]): comp[fd(a)] = fd(b)
    mk = {}
    for a in range(len(Pl)): mk.setdefault(fd(a), set()).update(pverts(Pl[a]))
    markers = list(mk.values())
    if not any(m & Vq for m in markers): continue
    for S in fronts_in(L):
        if len(S) < 5: continue
        EJ = L | S
        if any((l in phi) + (l in S) + pi[l] < 2 for l in EJ): continue
        VS = frozenset(v for l in S for v in ends(l))
        if not (VS & VP): continue
        nodes = markers + [Vq, VS]
        par = list(range(len(nodes)))
        def f2(a):
            while par[a] != a: par[a] = par[par[a]]; a = par[a]
            return a
        for a, b in itertools.combinations(range(len(nodes)), 2):
            if nodes[a] & nodes[b]: par[f2(a)] = f2(b)
        if len({f2(a) for a in range(len(nodes))}) == 1:
            found[(len(P), len(S))] += 1 if (VS & Vq) else 2      # a front vertex-disjoint from bd q can be p' in either layer
for n in range(3, NMAX + 1):
    print(f"D={D}: Z_({n},k), k >= 5:", dict(sorted((k, c) for (m, k), c in found.items() if m == n)))

# Referee (Section 16, seventh pass): the classes "4,4,k" (further fronts of 4 and k >= 9 edges) with n <= 2 plaquettes in d = 4, fronts first.
# Usage: python3 junc_bigfront.py FRONTS_DUMP k        (FRONTS_DUMP from fronts_redel.c D=3; e.g. k = 12, 13)
"""A candidate junction J contains p0 = (bd q, layer 0) and further front-pairs H (4 edges, so H = bd h for a plaquette h: Lemma mg:lem:smallfronts in
infinite volume) and G (k edges), and n <= 2 plaquettes P.  With A, B the unions of the fronts of layers 0, 1, Lemma mg:lem:junction(iii) gives
A xor B in E(P), so |A xor B| <= 4n <= 8.  The four layer assignments of (H, G):
  H:1, G:0  |A xor B| = |bd q| + |G| + |bd h| - 2|A n bd h| >= k > 8;     H:1, G:1  >= 4 + 4 + k - 2*4 = k > 8;     H:0, G:0  B empty, |A| = 8 + k > 8;
so H lies in layer 0 (vertex-disjoint from bd q) and G in layer 1, with |(bd q u bd h) n G| >= (8 + k - 4n)/2, hence |G n bd q| >= (8 + k - 4n)/2 - 4.
Enumeration: every placement of a k-edge front G containing at least that many links of bd q (each once: by its smallest link in bd q), every plaquette
h with bd h vertex-disjoint from bd q and enough links in G, and every set P of n <= 2 plaquettes containing a plaquette through the smallest link
of A xor B (which must lie in E(P)), the other plaquette having a link in bd q u bd h u G (a plaquette of P with no such link would need its four
links in the other plaquette); then all conditions are checked from scratch: n_J(e) >= 2 on
E(J), every front shares a vertex with P, the graph of markers (link-components of P) and front-pairs is connected (common vertices).  Each J is
counted once (a set of keys).  This is neither the authors' search nor the partition search: the big front is placed first."""
import sys, itertools
from collections import Counter

dump, k = sys.argv[1], int(sys.argv[2])
D = 3
E = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
add = lambda x, y: (x[0] + y[0], x[1] + y[1], x[2] + y[2])
sub = lambda x, y: (x[0] - y[0], x[1] - y[1], x[2] - y[2])
def ends(l): return (l[0], add(l[0], E[l[1]]))
def plinks(c):
    x, i, j = c
    return frozenset(((x, i), (add(x, E[i]), j), (add(x, E[j]), i), (x, j)))
def pverts(c):
    x, i, j = c
    return frozenset((x, add(x, E[i]), add(x, E[j]), add(add(x, E[i]), E[j])))
def plaqs_through(l):
    x, d = l; out = []
    for m in range(D):
        if m == d: continue
        i, j = min(d, m), max(d, m)
        out.append((x, i, j)); out.append((sub(x, E[m]), i, j))
    return out
def verts(F): return frozenset(v for l in F for v in ends(l))
O = (0, 0, 0); q = (O, 0, 1); dq = plinks(q); Vq = pverts(q)

classes = []
with open(dump) as fh:
    for line in fh:
        t = line.split()
        if int(t[0]) != k: continue
        L = []
        for a in range(k):
            x, y, z, d = map(int, t[1 + 4*a: 5 + 4*a]); L.append(((x, y, z), d))
        classes.append(L)
print(f"k = {k}: {len(classes)} classes of fronts loaded"); sys.stdout.flush()

found = {1: set(), 2: set()}
dq_sorted = sorted(dq)
nmax = 2
need_A = -(-(8 + k - 4*nmax) // 2)               # |(bd q u bd h) n G| >= need_A for n <= 2
need_q = max(1, need_A - 4)
if True:
    for L in classes:
        for e in L:
            for l in dq_sorted:
                if l[1] != e[1]: continue
                t = sub(l[0], e[0])
                G = frozenset((add(x, t), d) for x, d in L)
                inter = G & dq
                if len(inter) < need_q or min(inter) != l: continue          # each placement once, at its smallest link of bd q
                VG = verts(G)
                need_h = need_A - len(inter)
                # plaquettes h, bd h vertex-disjoint from bd q, with at least need_h links in G (need_h may be <= 0: then any h near G or q)
                hs = set()
                if need_h > 0:
                    for g in G:
                        for h in plaqs_through(g):
                            if len(plinks(h) & G) >= need_h and not (pverts(h) & Vq): hs.add(h)
                else:
                    raise SystemExit('need_h <= 0 not handled')
                for h in hs:
                    A = dq | plinks(h); Bset = G
                    W = A | Bset
                    cand = sorted({c for l2 in W for c in plaqs_through(l2)})
                    Delta = A ^ Bset                      # must lie in E(P): one plaquette of P contains its smallest link
                    if len(Delta) > 4*nmax: continue
                    first = plaqs_through(min(Delta))
                    Ps = {(c,) for c in first} | {tuple(sorted((c, c2))) for c in first for c2 in cand if c2 != c}
                    for P in sorted(Ps):
                        n = len(P)
                        if True:
                            EP = Counter(l2 for c in P for l2 in plinks(c))
                            EJ = W | set(EP)
                            if any((l2 in A) + (l2 in Bset) + EP.get(l2, 0) < 2 for l2 in EJ): continue
                            VP = [pverts(c) for c in P]
                            fr = [Vq, pverts(h), VG]
                            if not all(any(F & v for v in VP) for F in fr): continue
                            # markers: link-components of P
                            comp = list(range(n))
                            for a, b in itertools.combinations(range(n), 2):
                                if plinks(P[a]) & plinks(P[b]): comp[b] = comp[a] = min(comp[a], comp[b])
                            mk = {}
                            for a in range(n): mk.setdefault(comp[a], set()).update(VP[a])
                            nodes = list(mk.values()) + fr
                            par = list(range(len(nodes)))
                            def f2(a):
                                while par[a] != a: a = par[a]
                                return a
                            for a, b in itertools.combinations(range(len(nodes)), 2):
                                if nodes[a] & nodes[b]: par[f2(a)] = f2(b)
                            if len({f2(a) for a in range(len(nodes))}) != 1: continue
                            found[n].add((frozenset(P), h, G))
    print(f"  placements with n <= {nmax} done: " + ', '.join(f"n={n}: {len(found[n])}" for n in (1, 2))); sys.stdout.flush()
print(f"k = {k}: junctions with further fronts 4, {k} by n: " + ', '.join(f"{n}: {len(found[n])}" for n in (1, 2)))

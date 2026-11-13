# Section 16 of Part IV: Z_{2,k}, the junctions with two front-pairs p = (bd q, lambda), p' with |p'| = k >= 5, and two plaquettes.  Usage: python3 zexit_count.py D
"""By Lemma mg:lem:twofront(i), p' lies in the other layer and its front is contained in L = bd q u E(P) (at most 12 links); one plaquette of P shares a
vertex with q, and the other one with q or with the first.  For each such P and each subset of L that is a front with k >= 5 edges, the conditions of
Lemma mg:lem:junction are checked: n_J(e) >= 2 on E(J), both fronts share a vertex with a plaquette of P, and J is connected (markers = link-components
of P; adjacency = common vertex).  Pairs (J, p') = (P, front) are counted once.
Results: D=2: {6: 4, 7: 4, 8: 22, 9: 36, 10: 42, 11: 24, 12: 18};  D=3: {6: 20, 7: 12, 8: 182, 9: 296, 10: 606, 11: 456, 12: 458} (10 s).
The values for k = 6, 7 (and 8 in D=2) agree with three_front_count.py."""
import sys, itertools
from collections import Counter
from fronts_decision import make
D = int(sys.argv[1])
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
def nbrs(c):                              # plaquettes sharing a vertex with c (c included)
    return {x for v in pverts(c) for x in plaqs_at(v)}
q = (O, 0, 1); phi = plinks(q); Vq = pverts(q)
def is_front(S):
    if not S: return False
    deg = Counter(v for l in S for v in ends(l))
    if any(c < 2 for c in deg.values()): return False
    S = list(S); seen = {S[0]}; stack = [S[0]]
    while stack:
        l = stack.pop()
        for m in S:
            if m not in seen and set(ends(m)) & set(ends(l)): seen.add(m); stack.append(m)
    return len(seen) == len(S)
found = set()
P1 = nbrs(q)                              # includes q
for c1 in P1:
    for c2 in nbrs(c1) | P1:
        if c2 == c1: continue
        P = frozenset((c1, c2))
        if not any(pverts(c) & Vq for c in P): continue
        L = sorted(phi | plinks(c1) | plinks(c2))
        VP = pverts(c1) | pverts(c2)
        pi = Counter(l for c in P for l in plinks(c))
        for k in range(5, len(L) + 1):
            for S in itertools.combinations(L, k):
                S = frozenset(S)
                EJ = phi | S | plinks(c1) | plinks(c2)
                if any((l in phi) + (l in S) + pi[l] < 2 for l in EJ): continue
                if not is_front(S): continue
                VS = frozenset(v for l in S for v in ends(l))
                if not (VS & VP): continue
                # connectivity: markers (link-components of P), p, p'
                same = bool(plinks(c1) & plinks(c2))
                markers = [VP] if same else [pverts(c1), pverts(c2)]
                nodes = markers + [Vq, VS]
                par = list(range(len(nodes)))
                def f(a):
                    while par[a] != a: a = par[a]
                    return a
                for a, b in itertools.combinations(range(len(nodes)), 2):
                    if nodes[a] & nodes[b]: par[f(a)] = f(b)
                if len({f(a) for a in range(len(nodes))}) == 1: found.add((P, S))
print(f"D={D}: Z_(2,k), k >= 5:", dict(sorted(Counter(len(S) for P, S in found).items())))

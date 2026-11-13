# Referee (Section 16, third pass): Z_{2,k}, the pairs (J, p') with Phi(J) = {p, p'}, p = (bd q, s-1), |J|_C = 2, |p'| = k >= 5.
# Usage: python3 zexit_check.py D
"""From the definitions: n_J(e) >= 2 on E(J) (Lemma mg:lem:junction(iii)) forces every link of the exit phi' outside bd q to lie in E(P)
(only one front contains it), so phi' is a subset of L := bd q u E(P), |L| <= 12; and it forces the links of L that would otherwise be
covered once into phi' (a link of bd q with no plaquette of P on it, when p' is in the other layer; a link of E(P) in exactly one plaquette
and not in bd q).  J connected (5) and phi' inside L imply that q, c1, c2 are connected by common vertices, so (up to naming) c1 shares a
vertex with q and c2 with q or c1.  Both layers of p' are tried; every candidate is re-checked against (1)-(5) of junction_search.py.
Independent of zexit_count.py (which assumes the other layer for p'), written before reading it in detail."""
import sys, itertools
from collections import Counter
D = int(sys.argv[1])
E = [tuple(1 if k == i else 0 for k in range(D)) for i in range(D)]
add = lambda x, y: tuple(a + b for a, b in zip(x, y))
sub = lambda x, y: tuple(a - b for a, b in zip(x, y))
O = (0,) * D
def ends(l): return (l[0], add(l[0], E[l[1]]))
def plinks(c):
    x, i, j = c
    return frozenset(((x, i), (add(x, E[i]), j), (add(x, E[j]), i), (x, j)))
def pverts(c):
    x, i, j = c
    return frozenset((x, add(x, E[i]), add(x, E[j]), add(add(x, E[i]), E[j])))
def plaqs_at(v):
    out = []
    for i, j in itertools.combinations(range(D), 2):
        for a in (0, 1):
            for b in (0, 1):
                y = v
                if a: y = sub(y, E[i])
                if b: y = sub(y, E[j])
                out.append((y, i, j))
    return out
def is_front(S):
    deg = Counter(v for l in S for v in ends(l))
    if any(c < 2 for c in deg.values()): return False
    S = list(S); seen = {S[0]}; st = [S[0]]
    while st:
        l = st.pop()
        for m in S:
            if m not in seen and set(ends(m)) & set(ends(l)): seen.add(m); st.append(m)
    return len(seen) == len(S)
def connected(sets):
    par = list(range(len(sets)))
    def f(a):
        while par[a] != a: a = par[a]
        return a
    for a, b in itertools.combinations(range(len(sets)), 2):
        if sets[a] & sets[b]: par[f(a)] = f(b)
    return len({f(a) for a in range(len(sets))}) == 1
q = (O, 0, 1); dq = plinks(q); Vq = pverts(q)
near_q = set(c for v in Vq for c in plaqs_at(v))
found = set()
for c1 in near_q:
    for c2 in set(c for v in pverts(c1) for c in plaqs_at(v)) | near_q:
        if c2 == c1: continue
        P = frozenset((c1, c2))
        pi = Counter(l for c in P for l in plinks(c))
        L = dq | set(pi)
        VP = pverts(c1) | pverts(c2)
        for lam in (0, 1):
            # forced links of the exit
            R = set()
            for e in L:
                inA = e in dq
                base = inA + pi[e]                                   # coverage without the exit
                if base < 2:
                    if lam == 0 and inA: R = None; break               # same layer: a link of bd q cannot lie in a second front of layer s-1
                    R.add(e)
            if R is None: continue
            free = sorted(L - R)
            for r in range(0, len(free) + 1):
                for T in itertools.combinations(free, r):
                    S = frozenset(R) | frozenset(T)
                    if len(S) < 5: continue
                    if not is_front(S): continue
                    VS = frozenset(v for l in S for v in ends(l))
                    A = set(dq) | (S if lam == 0 else set()); B = S if lam == 1 else set()
                    EJ = A | B | set(pi)
                    if any((e in A) + (e in B) + pi[e] < 2 for e in EJ): continue
                    if lam == 0 and (VS & Vq): continue                # same-layer fronts vertex-disjoint
                    if not (VS & VP) or not (Vq & VP): continue        # both front-pairs touched by a marker
                    markers = [VP] if plinks(c1) & plinks(c2) else [pverts(c1), pverts(c2)]
                    if not connected(markers + [Vq, VS]): continue
                    found.add((P, S, lam))
print(f"D={D}: Z_(2,k) for k >= 5:", dict(sorted(Counter(len(S) for P, S, lam in found).items())),
      " layers of p':", dict(Counter('same' if lam == 0 else 'other' for P, S, lam in found)))

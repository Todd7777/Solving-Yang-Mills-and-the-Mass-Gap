# Section 16 of Part IV: junctions with three or more front-pairs of the next order, by exhaustive search.
# Usage: python3 next_order_count.py D NMAX SIZES   (SIZES: the sizes of the further front-pairs, e.g. 4,8 or 4,4,4)
"""Candidate junctions at a slice s containing p = (bd q, s-1), q a fixed plaquette of Z^D (D = d-1), with at most NMAX plaquettes and further
front-pairs of exactly the given sizes.  Conditions (Lemma mg:lem:junction, Definition mg:def:junction), as in three_front_count.py:
(1) n_J(e) >= 2 on E(J); (2) fronts of one layer vertex-disjoint, each a front (connected, all degrees >= 2); (3) every front shares a vertex with a
plaquette of P; (4) J connected (markers = link-components of P, and the front-pairs; adjacency = common vertex).
Search: deficiency-driven (at the smallest link covered once, branch over the further elements containing it), with seeding of further plaquettes or
fronts at the vertices of a closed state; fronts from the decision search of fronts_decision.py, restricted to the sizes still to be placed.
Output: the number of junctions by (n, sorted sizes of the further fronts, layers relative to p).
Results (the next order S_6 of Lemma mg:lem:corners; p in layer s-1):
  d=3 (D=2):  n=4, 4,4: 0;  n=3, 4,6: 24;  n=3, 4,7: 24;  n=2, 4,8: 48;  n=2, 4,4,4: 0;  n=2, 4,9: 16   (each well under a minute)
  d=4 (D=3):  n=4, 4,4: 42; n=3, 4,6: 680; n=3, 4,7: 344; n=2, 4,8: 992; n=2, 4,4,4: 4; n=2, 4,9: 256 (up to 7 min)
  and, as checks, N_3 = 48, 876, N_2 = 24, 320, N'_2 = 24, 312 (Lemma mg:lem:threelow) with their layers."""
import sys, itertools
from collections import Counter, defaultdict
from fronts_decision import fronts_through_origin, make
D = int(sys.argv[1]); NMAX = int(sys.argv[2]); SIZES = tuple(sorted(int(x) for x in sys.argv[3].split(',')))
E, add, sub, inc, ends = make(D)
O = (0,) * D
KF = max(SIZES)
def plinks(p):
    x, i, j = p
    return ((x, i), (add(x, E[i]), j), (add(x, E[j]), i), (x, j))
def pverts(p):
    x, i, j = p
    return frozenset((x, add(x, E[i]), add(x, E[j]), add(add(x, E[i]), E[j])))
def plaqs_through(l):
    x, k = l; out = []
    for m in range(D):
        if m == k: continue
        i, j = min(k, m), max(k, m)
        out.append((x, i, j)); out.append((sub(x, E[m]), i, j))
    return out
def plaqs_at(v):
    out = []
    for i, j in itertools.combinations(range(D), 2):
        for si in (0, 1):
            for sj in (0, 1):
                y = sub(v, E[i]) if si else v
                y = sub(y, E[j]) if sj else y
                out.append((y, i, j))
    return out
def fverts(F): return frozenset(v for l in F for v in ends(l))
def translate(F, x): return frozenset((add(l[0], x), l[1]) for l in F)
FO = [F for F in fronts_through_origin(D, KF) if len(F) in SIZES]
FL = defaultdict(set)          # fronts (of the wanted sizes) containing the link (O, i)
FV = set()                     # fronts (of the wanted sizes) through the origin, all translates
for F in FO:
    for v in fverts(F):
        G = translate(F, sub(O, v)); FV.add(G)
for G in FV:
    for v in fverts(G):
        H = translate(G, sub(O, v))
        for i in range(D):
            if (O, i) in H: FL[i].add(H)
FL = {i: list(s) for i, s in FL.items()}; FV = list(FV)
q = (O, 0, 1); p0 = (frozenset(plinks(q)), -1)
def remaining(F):
    rem = Counter(SIZES)
    for G, lam in F:
        if (G, lam) != p0: rem[len(G)] -= 1
    return rem
found = set(); seen = set()
def can_add(F, G, lam):
    if (G, lam) in F: return False
    gv = fverts(G)
    return all(not (gv & fverts(H)) for H, mu in F if mu == lam)
def finish(P, F):
    if not P or sum(remaining(F).values()) != 0: return
    PV = [pverts(c) for c in P]
    if not all(any(fverts(G) & pv for pv in PV) for G, lam in F): return
    Pl = list(P); comp = list(range(len(Pl)))
    def fd(a):
        while comp[a] != a: comp[a] = comp[comp[a]]; a = comp[a]
        return a
    for a, b in itertools.combinations(range(len(Pl)), 2):
        if set(plinks(Pl[a])) & set(plinks(Pl[b])): comp[fd(a)] = fd(b)
    mk = {}
    for a in range(len(Pl)): mk.setdefault(fd(a), set()).update(pverts(Pl[a]))
    nodes = list(mk.values()) + [fverts(G) for G, lam in F]
    par = list(range(len(nodes)))
    def f2(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    for a, b in itertools.combinations(range(len(nodes)), 2):
        if nodes[a] & nodes[b]: par[f2(a)] = f2(b)
    if len({f2(a) for a in range(len(nodes))}) == 1: found.add((P, F))
def rec(cnt, P, F, excl):
    rem = remaining(F)
    defi = [l for l, c in cnt.items() if c == 1]
    if not defi:
        key = (P, F)
        if key in seen: return
        seen.add(key)
        finish(P, F)
        V = set().union(*[fverts(G) for G, lam in F]) | set().union(*[pverts(c) for c in P]) if P else set().union(*[fverts(G) for G, lam in F])
        seeds = set()
        for v in V:
            if len(P) < NMAX:
                for c in plaqs_at(v):
                    if c not in P: seeds.add(('r', c))
            if sum(rem.values()) > 0:
                for G in FV:
                    if rem[len(G)] <= 0: continue
                    H = translate(G, v)
                    for lam in (-1, 1):
                        if can_add(F, H, lam): seeds.add(('f', (H, lam)))
        for kind, obj in seeds:
            cnt2 = Counter(cnt)
            if kind == 'r':
                for m in plinks(obj): cnt2[m] += 1
                rec(cnt2, P | {obj}, F, frozenset())
            else:
                for m in obj[0]: cnt2[m] += 1
                rec(cnt2, P, F | {obj}, frozenset())
        return
    l = min(defi)
    creal = [c for c in plaqs_through(l) if c not in P and c not in excl]
    layers_with_l = {lam for G, lam in F if l in G}
    cfr = []
    if sum(rem.values()) > 0:
        for lam in (-1, 1):
            if lam in layers_with_l: continue
            for G in FL.get(l[1], []):
                if rem[len(G)] <= 0: continue
                H = translate(G, l[0])
                if can_add(F, H, lam): cfr.append((H, lam))
    room = NMAX - len(P)
    fopts = [()] + [(x,) for x in cfr]
    for x, y in itertools.combinations(cfr, 2):
        if x[1] != y[1]:
            r2 = Counter(rem); r2[len(x[0])] -= 1; r2[len(y[0])] -= 1
            if min(r2.values()) >= 0 and can_add(F | {x}, *y): fopts.append((x, y))
    for r in range(0, min(room, len(creal)) + 1):
        for T in itertools.combinations(creal, r):
            for TF in fopts:
                if not T and not TF: continue
                cnt2 = Counter(cnt)
                for c in T:
                    for m in plinks(c): cnt2[m] += 1
                for G, lam in TF:
                    for m in G: cnt2[m] += 1
                excl2 = excl | frozenset(c for c in creal if c not in T)
                rec(cnt2, P | frozenset(T), F | frozenset(TF), excl2)
rec(Counter(plinks(q)), frozenset(), frozenset([p0]), frozenset())
stat = Counter()
for P, F in found:
    others = sorted((len(G), 'same' if lam == -1 else 'other') for G, lam in F if (G, lam) != p0)
    stat[(len(P), tuple(others))] += 1
print(f"D={D}, NMAX={NMAX}, sizes={SIZES}: {len(FO)} fronts through the origin of these sizes")
for k in sorted(stat): print('  n=%d, further fronts (size, layer rel. p) = %s: %d' % (k[0], k[1], stat[k]))
by_n = Counter()
for (n, oth), c in stat.items(): by_n[n] += c
print('  totals by n:', dict(sorted(by_n.items())))

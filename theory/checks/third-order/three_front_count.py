# Section 16 of Part IV (Lemma mg:lem:threelow): junctions with three front-pairs, by exhaustive search.  Usage: python3 three_front_count.py D NMAX KF FX
"""Junctions at a slice s containing the front-pair p = (bd q, s-1), q a fixed plaquette of Z^D (D = d-1), with at most NMAX plaquettes and at most FX further
front-pairs, all fronts with at most KF edges.  The conditions are the necessary ones of Lemma mg:lem:junction and Definition mg:def:junction:
(1) every link of E(J) lies in at least two of: the fronts of layer s-1, the fronts of layer s, the plaquettes of P (n_J(e) >= 2);
(2) the fronts of one layer are vertex-disjoint (distinct components of T_lambda), and each is a front in Gamma (connected, every vertex of degree >= 2);
(3) every front shares a vertex with a plaquette of P;  (4) J is connected in H_s (markers = link-components of P, and the front-pairs; adjacency =
common vertex, including two front-pairs of different layers).
Search: at the smallest link covered once, branch over the further elements containing it (plaquettes through it; fronts through it in a layer
in which no front contains it yet); at a closed state (no link covered once) also seed a further plaquette or front at a vertex of the current
elements, since a junction may consist of closed parts touching only at vertices.  Distinct junctions are counted once (as sets (P, front-pairs)).
At the end the families of Lemma mg:lem:threelow -- (i) caps and cubes, (ii) (a), (b), (c) -- are built from their description and compared with
the junctions found (as sets, not only in number).
Results (p in layer s-1; the other front sizes listed):
  D=2, NMAX=3, KF=6, FX=2 (1 s): n=3 with fronts (4,4): 48;  n=2 with (4,6): 24;  both families equal to the search.
  D=3, NMAX=3, KF=6, FX=2 (3 min): n=3 with (4,4): 876;  n=2 with (4,6): 320;  both families equal to the search.
  D=2, 3, NMAX=2, KF=7, FX=2: n=2 with (4,7): 24 and 312, equal to the families (iii); no junction with n<=2, three front-pairs, two of four edges and
  the third of fewer than six (Lemma mg:lem:betafour(iv)), and none with n=1 and three front-pairs (fronts <= 7; Lemma mg:lem:betafour(iii)).
  Other counts printed, e.g. D=3: Z_(1,4)=1, Z_(2,4)=32, Z_(3,4)=0, K_S=12 six-edge exits with one plaquette."""
import sys, itertools
from collections import Counter, defaultdict
D = int(sys.argv[1]); NMAX = int(sys.argv[2]); KF = int(sys.argv[3]); FX = int(sys.argv[4])
E = [tuple(1 if k == i else 0 for k in range(D)) for i in range(D)]
O = (0,) * D
def add(x, y): return tuple(a + b for a, b in zip(x, y))
def sub(x, y): return tuple(a - b for a, b in zip(x, y))
def lv(l): return (l[0], add(l[0], E[l[1]]))
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
                out.append((sub(sub(v, E[i]) if si else v, E[j]) if sj else (sub(v, E[i]) if si else v), i, j))
    return out
def links_at(v): return [(v, i) for i in range(D)] + [(sub(v, E[i]), i) for i in range(D)]
def fverts(F): return frozenset(v for l in F for v in lv(l))
def is_front(S):
    deg = Counter()
    for l in S:
        for v in lv(l): deg[v] += 1
    return all(c >= 2 for c in deg.values())
# fronts through the origin with at most KF links (connected link sets with all degrees >= 2), by Redelmeier from each link at the origin
def fronts_through_origin():
    roots = sorted(links_at(O)); res = set()
    def nbrs(l): return [m for v in lv(l) for m in links_at(v) if m != l]
    for r in roots:
        S = [r]
        def rec(untried, seen):
            if is_front(S):
                fs = frozenset(S)
                if min(m for m in fs if O in lv(m)) == r: res.add(fs)
            if len(S) == KF: return
            untried = list(untried)
            while untried:
                m = untried.pop(); S.append(m)
                new = [x for x in nbrs(m) if x not in seen]
                rec(untried + new, seen | set(new))
                S.pop()
        nb = [x for x in nbrs(r)]
        rec(nb, {r} | set(nb) | set(roots[:roots.index(r)]))
    return res
FO = fronts_through_origin()
def translate(F, x): return frozenset((add(l[0], x), l[1]) for l in F)
FL = defaultdict(list)                      # fronts containing the link (O, i), for each direction i
for F in FO:
    for v in fverts(F):
        G = translate(F, sub(O, v))         # all translates through the origin
        for i in range(D):
            if (O, i) in G and G not in FL[i]: FL[i].append(G)
def fronts_through_link(l): return [translate(G, l[0]) for G in FL[l[1]]]
FV = [translate(F, sub(O, v)) for F in FO for v in fverts(F)]   # fronts through the origin, as all translates
FV = list(set(G for G in FV if O in fverts(G)))
def fronts_through_vertex(v): return [translate(G, v) for G in FV]
q = (O, 0, 1); p0 = (frozenset(plinks(q)), -1)
found = {}
def finish(P, F):
    if not P: return
    PV = [pverts(c) for c in P]
    for (G, lam) in F:
        gv = fverts(G)
        if not any(gv & pv for pv in PV): return
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
    if len({f2(a) for a in range(len(nodes))}) != 1: return
    key = (P, F)
    if key not in found: found[key] = True
seen = set()
def can_add_front(F, G, lam):
    if (G, lam) in F: return False
    gv = fverts(G)
    return all(not (gv & fverts(H)) for H, mu in F if mu == lam)
def rec(cnt, P, F, excl):
    defi = [l for l, c in cnt.items() if c == 1]
    if not defi:
        key = (P, F)
        if key in seen: return
        seen.add(key)
        finish(P, F)
        V = set(fverts(frozenset().union(*[G for G, lam in F]))) | set().union(*[pverts(c) for c in P])
        seeds = set()
        for v in V:
            if len(P) < NMAX:
                for c in plaqs_at(v):
                    if c not in P: seeds.add(('r', c))
            if len(F) - 1 < FX:
                for G in fronts_through_vertex(v):
                    for lam in (-1, 1):
                        if can_add_front(F, G, lam): seeds.add(('f', (G, lam)))
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
    if len(F) - 1 < FX:
        for lam in (-1, 1):
            if lam in layers_with_l: continue
            for G in fronts_through_link(l):
                if can_add_front(F, G, lam): cfr.append((G, lam))
    room = NMAX - len(P); froom = FX - (len(F) - 1)
    for r in range(0, min(room, len(creal)) + 1):
        for T in itertools.combinations(creal, r):
            fopts = [()]
            if froom >= 1: fopts += [(x,) for x in cfr]
            if froom >= 2: fopts += [(x, y) for x, y in itertools.combinations(cfr, 2) if x[1] != y[1] and not (fverts(x[0]) & fverts(y[0]) and False)]
            for TF in fopts:
                if not T and not TF: continue
                if len(TF) == 2 and not can_add_front(F | {TF[0]}, *TF[1]): continue
                cnt2 = Counter(cnt)
                for c in T:
                    for m in plinks(c): cnt2[m] += 1
                for G, lam in TF:
                    for m in G: cnt2[m] += 1
                excl2 = excl | frozenset(c for c in creal if c not in T)
                rec(cnt2, P | frozenset(T), F | frozenset(TF), excl2)
rec(Counter(plinks(q)), frozenset(), frozenset([p0]), frozenset())
stat = Counter()
for (P, F) in found:
    sizes = tuple(sorted(len(G) for G, lam in F if (G, lam) != p0))
    stat[(len(P), len(F), sizes)] += 1
print(f"D={D}, NMAX={NMAX}, KF={KF}, FX={FX}: {len(FO)} fronts through the origin")
for k in sorted(stat): print('  n=%d, front-pairs=%d, other front sizes=%s: %d' % (k[0], k[1], k[2], stat[k]))

# ---- the families of Lemma mg:lem:threelow, built from their description, compared with the search ----
def link_of(u, v):
    for i in range(D):
        if add(u, E[i]) == v: return (u, i)
        if add(v, E[i]) == u: return (v, i)
    return None
def nbhd(R):
    return [(x, i, j) for x in itertools.product(range(-R, R + 1), repeat=D) for i, j in itertools.combinations(range(D), 2)]
PL = nbhd(4)
bd = lambda c: frozenset(plinks(c))
touch = lambda a, c: a != c and bool(pverts(a) & pverts(c))
Nq = [c for c in PL if touch(c, q)]
def opposite_edge(c, l):
    """the edge of plaquette c opposite to its edge l"""
    (u, v) = lv(l)
    for m in plinks(c):
        if not (set(lv(m)) & {u, v}): return m
def cubes_of(c):
    x, i, j = c; out = []
    for k in range(D):
        if k in (i, j): continue
        for s in (0, 1):
            y = sub(x, E[k]) if s else x
            out.append((y, (i, j, k)))
    return out
def faces(cube):
    y, dirs = cube; out = []
    for a, b in itertools.combinations(sorted(dirs), 2):
        k = [t for t in dirs if t not in (a, b)][0]
        out.append((y, a, b)); out.append((add(y, E[k]), a, b))
    return out
def opp_face(cube, c):
    for g in faces(cube):
        if g != c and not (pverts(g) & pverts(c)): return g
fam_i, fam_ii = set(), set()
up, dn = 1, -1
for a, c in itertools.combinations(Nq, 2):                         # caps, p alone in its layer
    if not (pverts(a) & pverts(c)): fam_i.add((frozenset([q, a, c]), frozenset([p0, (bd(a), up), (bd(c), up)])))
for a in PL:                                                        # caps, p with (bd a, s-1)
    if pverts(a) & pverts(q): continue
    for b in Nq:
        if b != a and touch(b, a): fam_i.add((frozenset([q, a, b]), frozenset([p0, (bd(a), dn), (bd(b), up)])))
for cube in cubes_of(q):                                            # cubes
    F6 = faces(cube); qh = opp_face(cube, q); sides = [g for g in F6 if g not in (q, qh)]
    for sd in sides:
        fam_i.add((frozenset(g for g in F6 if g not in (q, qh, sd)), frozenset([p0, (bd(qh), dn), (bd(sd), up)])))
    for sd in sides:
        sh = opp_face(cube, sd)
        if sd < sh: fam_i.add((frozenset(g for g in F6 if g not in (q, sd, sh)), frozenset([p0, (bd(sd), up), (bd(sh), up)])))
for l in plinks(q):                                                 # (a), (b) with r = q: c shares the link l = e with q
    for c in plaqs_through(l):
        if c == q: continue
        f = opposite_edge(c, l); Sig = bd(q) ^ bd(c)
        for rp in plaqs_through(f):
            if rp != c: fam_ii.add((frozenset([c, rp]), frozenset([p0, (bd(rp), dn), (Sig, up)])))
        for v in lv(f):
            for rp in plaqs_at(v):
                if rp != c and len(pverts(rp) & pverts(c)) == 1: fam_ii.add((frozenset([c, rp]), frozenset([p0, (bd(rp), dn), (Sig, up)])))
for f in plinks(q):                                                 # (a) with r = q': q contains the edge f of c opposite to e = c n q'
    for c in plaqs_through(f):
        if c == q: continue
        e = opposite_edge(c, f)
        for qp in plaqs_through(e):
            if qp != c: fam_ii.add((frozenset([q, c]), frozenset([p0, (bd(qp), dn), (bd(qp) ^ bd(c), up)])))
for v in pverts(q):                                                 # (b) with r = q': q shares exactly the vertex v with c, v an endpoint of f
    for c in plaqs_at(v):
        if c == q or len(pverts(c) & pverts(q)) != 1: continue
        for f in plinks(c):
            if v not in lv(f): continue
            e = opposite_edge(c, f)
            for qp in plaqs_through(e):
                if qp != c: fam_ii.add((frozenset([q, c]), frozenset([p0, (bd(qp), dn), (bd(qp) ^ bd(c), up)])))
for cube in cubes_of(q):                                            # (c)
    F6 = faces(cube); qh = opp_face(cube, q); sides = [g for g in F6 if g not in (q, qh)]
    for c1, c2 in itertools.combinations(sides, 2):
        if set(plinks(c1)) & set(plinks(c2)):
            fam_ii.add((frozenset(g for g in sides if g not in (c1, c2)), frozenset([p0, (bd(qh), dn), (bd(c1) ^ bd(c2), up)])))
fam_iii = set()                                                     # (iii): fronts 4, 4, 7 with Theta = bd r u bd c in place of Sigma
for l in plinks(q):
    for c in plaqs_through(l):
        if c == q: continue
        f = opposite_edge(c, l); Th = bd(q) | bd(c)
        for rp in plaqs_through(f):
            if rp != c: fam_iii.add((frozenset([c, rp]), frozenset([p0, (bd(rp), dn), (Th, up)])))
        for v in lv(f):
            for rp in plaqs_at(v):
                if rp != c and len(pverts(rp) & pverts(c)) == 1: fam_iii.add((frozenset([c, rp]), frozenset([p0, (bd(rp), dn), (Th, up)])))
for f in plinks(q):
    for c in plaqs_through(f):
        if c == q: continue
        e = opposite_edge(c, f)
        for qp in plaqs_through(e):
            if qp != c: fam_iii.add((frozenset([q, c]), frozenset([p0, (bd(qp), dn), (bd(qp) | bd(c), up)])))
for v in pverts(q):
    for c in plaqs_at(v):
        if c == q or len(pverts(c) & pverts(q)) != 1: continue
        for f in plinks(c):
            if v not in lv(f): continue
            e = opposite_edge(c, f)
            for qp in plaqs_through(e):
                if qp != c: fam_iii.add((frozenset([q, c]), frozenset([p0, (bd(qp), dn), (bd(qp) | bd(c), up)])))
found_iii = {k for k in found if len(k[0]) == 2 and len(k[1]) == 3 and sorted(len(G) for G, lam in k[1]) == [4, 4, 7]}
found_i = {k for k in found if len(k[0]) == 3 and len(k[1]) == 3 and all(len(G) == 4 for G, lam in k[1])}
found_ii = {k for k in found if len(k[0]) == 2 and len(k[1]) == 3 and sorted(len(G) for G, lam in k[1]) == [4, 4, 6]}
if NMAX >= 3 and FX >= 2:
    print(f"families (i): built {len(fam_i)}, found {len(found_i)}, equal: {fam_i == found_i}")
if NMAX >= 2 and FX >= 2 and KF >= 6:
    print(f"families (ii): built {len(fam_ii)}, found {len(found_ii)}, equal: {fam_ii == found_ii}")
if NMAX >= 2 and FX >= 2 and KF >= 7:
    print(f"families (iii): built {len(fam_iii)}, found {len(found_iii)}, equal: {fam_iii == found_iii}")

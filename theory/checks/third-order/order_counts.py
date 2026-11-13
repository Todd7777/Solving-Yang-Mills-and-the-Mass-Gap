# Section 16 of Part IV (Lemma mg:lem:corners, Remark mg:rem:nextorder): the exact part of M_S for S(m) = {(n,E): n + E <= m + 1/2}, generated
# automatically.  Usage: python3 order_counts.py D m [out.json]
"""For every junction class (n plaquettes; further front-pairs of sizes K, |K| >= 2) with a pair (J, p') at n + E <= m + 1/2, where for p' of size k
E = (k - 4) + (sum K - k)/2, count the candidate junctions containing p = (bd q, s-1) by the exhaustive search of next_order_count.py (same
conditions), and list the pairs as terms (count, n, E, |p'|, |K| - 1): weight count b^n a^E e^{(|K|-1) c_J} (f(k) a^{4-k} if k > 4).  The corners of
S(m) are (n, m + 1 - n), 2 <= n <= m - 1, with E >= 2 (>= 3 if n = 2).  There are no fronts with five edges (sides >= 8)."""
import sys, json, itertools
from collections import Counter, defaultdict
from fractions import Fraction
from fronts_decision import fronts_through_origin, make

PRUNE = True   # the capacity prune: a necessary condition (links covered once need a later element through them), so it changes no count
def count_junctions(D, NMAX, SIZES, front_cache={}):
    """number of candidate junctions with n <= NMAX plaquettes and further fronts of exactly the sizes SIZES, by n"""
    SIZES = tuple(sorted(SIZES))
    E, add, sub, inc, ends = make(D)
    O = (0,) * D
    KF = max(SIZES)
    if (D, KF) not in front_cache: front_cache[(D, KF)] = fronts_through_origin(D, KF)
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
    fverts = lambda F: frozenset(v for l in F for v in ends(l))
    translate = lambda F, x: frozenset((add(l[0], x), l[1]) for l in F)
    FO = [F for F in front_cache[(D, KF)] if len(F) in SIZES]
    FV = set()
    for F in FO:
        for v in fverts(F): FV.add(translate(F, sub(O, v)))
    FL = defaultdict(set)
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
    def capacity(P, rem):
        # every link covered once must be covered again: by a later plaquette (at most 4 links each) or a later front (at most its size)
        return 4*(NMAX - len(P)) + sum(k*c for k, c in rem.items() if c > 0)
    def rec(cnt, P, F, excl):
        rem = remaining(F)
        defi = [l for l, c in cnt.items() if c == 1]
        if PRUNE and len(defi) > capacity(P, rem): return
        if not defi:
            key = (P, F)
            if key in seen: return
            seen.add(key)
            finish(P, F)
            V = set().union(*[fverts(G) for G, lam in F]) | (set().union(*[pverts(c) for c in P]) if P else set())
            seeds = set()
            U = set(cnt)
            capf = capacity(P, rem)
            for v in V:
                if len(P) < NMAX:
                    for c in plaqs_at(v):
                        if c not in P: seeds.add(('r', c))
                if sum(rem.values()) > 0:
                    for G in FV:
                        if rem[len(G)] <= 0: continue
                        H = translate(G, v)
                        if PRUNE and len(H - U) > capf - len(H): continue     # the new links of H could not all be covered again
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
            U = set(cnt); D1 = set(defi); capf = capacity(P, rem)
            for lam in (-1, 1):
                if lam in layers_with_l: continue
                for G in FL.get(l[1], []):
                    if rem[len(G)] <= 0: continue
                    H = translate(G, l[0])
                    # links covered once after adding H: the old ones not in H and the new ones of H; each needs a later element
                    if PRUNE and len(D1) - len(H & D1) + len(H - U) > capf - len(H): continue
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
                    rec(cnt2, P | frozenset(T), F | frozenset(TF), excl | frozenset(c for c in creal if c not in T))
    rec(Counter(plinks(q)), frozenset(), frozenset([p0]), frozenset())
    return Counter(len(P) for P, F in found)

def E_of(K, k):
    return Fraction(k - 4) + Fraction(sum(K) - k, 2)

def classes(m):
    """the classes (n, K) with some pair at n + E <= m + 1/2, and their maximal n"""
    lim = Fraction(2 * m + 1, 2)
    sizes = [4] + list(range(6, 2 * m + 1))
    out = defaultdict(int)                 # K -> largest n needed
    for r in range(2, m):
        for K in itertools.combinations_with_replacement(sizes, r):
            Emin = E_of(K, K[0])
            for n in range(2, m):
                if n + Emin <= lim and (n > 2 or Emin >= 3): out[K] = max(out[K], n)
    return out

if __name__ == "__main__":
    D, m = int(sys.argv[1]), int(sys.argv[2])
    lim = Fraction(2 * m + 1, 2)
    terms = []; table = []
    for K, nmax in sorted(classes(m).items(), key=lambda t: (len(t[0]), t[0])):
        cnt = count_junctions(D, nmax, K)
        for n in range(2, nmax + 1):
            N = cnt.get(n, 0)
            for k, mult in Counter(K).items():
                Ek = E_of(K, k)
                if n + Ek <= lim and (n > 2 or E_of(K, K[0]) >= 3):
                    table.append((n, K, k, float(Ek), N))
                    if N: terms.append((mult * N, n, float(Ek), k, len(K) - 1))
        print(f"K={K} (n<={nmax}): {dict(sorted(cnt.items()))}"); sys.stdout.flush()
    agg = Counter()
    for c, n, Ek, k, r1 in terms: agg[(n, Ek, k, r1)] += c
    terms = [(c, n, Ek, k, r1) for (n, Ek, k, r1), c in sorted(agg.items())]
    corners = [(n, m + 1 - n) for n in range(2, m) if (m + 1 - n) >= (3 if n == 2 else 2)]
    print("terms (count, n, E, |p'|, |K|-1):", terms)
    print("corners:", corners)
    if len(sys.argv) > 3: json.dump({'terms': terms, 'corners': corners}, open(sys.argv[3], 'w'))

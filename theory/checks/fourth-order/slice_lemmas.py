# Fourth order (Part IV, Section 12): exhaustive checks of the combinatorial lemmas behind the upper bound to fourth order.
#  A. (Six plaquettes) A nonempty finite set Y of plaquettes of Z^n in which no link lies in exactly one plaquette of Y has
#     |Y| >= 6, with equality only for the six faces of a unit cube.  Enumerated: all such Y with |Y| <= 7 containing a fixed
#     plaquette, n = 2, 3, 4 (backtracking over a link that lies in exactly one plaquette).
#  B. (Slice configurations, four-edge front) F = boundary of a plaquette q; q' any plaquette near q; X a set of plaquettes, q' not
#     in X, with n(e) := [e in F] + [e in q'] + #{p in X: e in p} != 1 for every link e.  Claim: |X| <= 4 only for (q' = q, X empty)
#     and (q, q' faces of a unit cube, X its four other faces).
#  C. (Slice configurations, six-edge front) F = bd(q) xor bd(c), q and c sharing a link (coplanar or perpendicular); same
#     condition.  Claim: |X| <= 2 only for q' in {q, c} and X = {the other one}.
import itertools, sys
from collections import Counter

def unit(n, i): return tuple(1 if j == i else 0 for j in range(n))
def add(x, y): return tuple(a + b for a, b in zip(x, y))

def plaq_edges(x, i, j, n):
    ei, ej = unit(n, i), unit(n, j)
    return frozenset([(x, i), (add(x, ei), j), (add(x, ej), i), (x, j)])

def plaqs_at_link(link, n):
    """the 2(n-1) plaquettes containing the link (x, i)."""
    x, i = link; out = []
    for j in range(n):
        if j == i: continue
        ej = unit(n, j)
        for base in (x, add(x, tuple(-c for c in ej))):
            a, b = min(i, j), max(i, j)
            out.append((base, a, b))
    return out

class Lat:
    def __init__(self, n): self.n = n; self.cache = {}
    def E(self, p):
        if p not in self.cache: self.cache[p] = plaq_edges(p[0], p[1], p[2], self.n)
        return self.cache[p]

def cover_count(L, plist, extra=()):
    c = Counter()
    for p in plist:
        for e in L.E(p): c[e] += 1
    for e in extra: c[e] += 1
    return c

def enumerate_sets(L, fixed_plaq, base_count, nmax, forbidden=()):
    """all sets X (|X| <= nmax) of plaquettes not in `forbidden` such that base_count + cover(X) has no link with count 1.
    Complete for nmax <= 5: at a state with a link of count 1 every valid superset contains a plaquette through it; at a valid
    state every valid proper superset with at most nmax plaquettes contains a plaquette sharing a link with the current support
    (otherwise the added plaquettes would form a nonempty set without links of count 1 by themselves, which has >= 6 plaquettes
    by part A)."""
    found = set(); seen = set()
    def rec(X):
        key = frozenset(X)
        if key in seen: return
        seen.add(key)
        c = Counter(base_count)
        for p in X:
            for e in L.E(p): c[e] += 1
        bad = sorted(e for e, v in c.items() if v == 1)
        if not bad: found.add(key)
        if len(X) >= nmax: return
        if bad:
            cand = [p for p in plaqs_at_link(bad[0], L.n)]
        else:
            cand = sorted({p for e, v in c.items() if v > 0 for p in plaqs_at_link(e, L.n)})
        for p in cand:
            if p in X or p in forbidden: continue
            rec(X | {p})
    rec(frozenset())
    return found

def is_cube(L, Y):
    if len(Y) != 6: return False
    verts = set()
    for p in Y:
        for (x, i) in L.E(p):
            verts.add(x); verts.add(add(x, unit(L.n, i)))
    if len(verts) != 8: return False
    span = [max(v[k] for v in verts) - min(v[k] for v in verts) for k in range(L.n)]
    return sorted(span) == [0]*(L.n - 3) + [1, 1, 1]

def partA(n, nmax=7):
    L = Lat(n); o = tuple([0]*n); p0 = (o, 0, 1)
    sols = enumerate_sets(L, None, cover_count(L, [p0]), nmax - 1, forbidden=(p0,))
    sizes = Counter(len(X) + 1 for X in sols)
    cubes = sum(1 for X in sols if is_cube(L, set(X) | {p0}))
    print(f"A. Z^{n}: 2-covered sets containing a fixed plaquette, size <= {nmax}: sizes {dict(sorted(sizes.items()))}; "
          f"of size 6: {sizes.get(6, 0)}, all unit cubes: {cubes == sizes.get(6, 0)} (number of unit cubes containing it: {2*(n-2)})")
    return sizes

def near_plaquettes(n, R=2):
    out = []
    for x in itertools.product(range(-R, R + 2), repeat=n):
        if any(abs(c) > R + 1 for c in x[2:]) : pass
        for i, j in itertools.combinations(range(n), 2):
            out.append((x, i, j))
    return out

def partB(n, nmax=4):
    L = Lat(n); o = tuple([0]*n); q = (o, 0, 1); F = L.E(q)
    tally = Counter(); bad = []
    for qp in near_plaquettes(n):
        base = cover_count(L, [qp], extra=F)
        sols = enumerate_sets(L, None, base, nmax, forbidden=(qp,))
        for X in sols:
            if qp == q and len(X) == 0: tally['q\'=q, X empty'] += 1; continue
            Y = set(X) | {q, qp}
            if qp != q and q not in X and len(X) == 4 and is_cube(L, Y): tally['cube partner, X = four faces'] += 1; continue
            bad.append((qp, sorted(X)))
    print(f"B. Z^{n}, four-edge front: configurations with |X| <= {nmax}: {dict(tally)}; others: {len(bad)}"
          f"  (expected cube partners: {5*2*(n-2)})")
    for b in bad[:5]: print("   unexpected:", b)

def partC(n, nmax=2):
    L = Lat(n); o = tuple([0]*n); q = (o, 0, 1)
    cs = [p for p in plaqs_at_link((o, 0), n) if p != q]           # plaquettes sharing the link (o, 0) with q
    for c in cs:
        kind = 'coplanar' if c[1:] == q[1:] else 'perpendicular'
        F = L.E(q) ^ L.E(c)
        tally = Counter(); bad = []
        for qp in near_plaquettes(n):
            base = cover_count(L, [qp], extra=F)
            sols = enumerate_sets(L, None, base, nmax, forbidden=(qp,))
            for X in sols:
                if qp in (q, c) and set(X) == ({q, c} - {qp}): tally['main'] += 1; continue
                bad.append((qp, sorted(X)))
        print(f"C. Z^{n}, six-edge front ({kind}, c={c}): configurations with |X| <= {nmax}: {dict(tally)}; others: {len(bad)}")
        for b in bad[:5]: print("   unexpected:", b)

if __name__ == "__main__":
    for n in (2, 3, 4):
        partA(n, 7 if n < 4 else 6)
    for n in (2, 3, 4):
        partB(n)
    for n in (2, 3, 4):
        partC(n)

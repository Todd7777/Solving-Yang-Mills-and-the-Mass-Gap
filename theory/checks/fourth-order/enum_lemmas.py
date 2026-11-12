# enum_lemmas.py -- written by the independent referee of Section 12 of Part IV (2026-11-12), added to the repository with credit;
# unchanged apart from this header. A compass, not part of the proofs.
# Referee code (independent algorithm): exhaustive checks of Lemma 12.5 (six plaquettes) and Lemma 12.6 (configurations at
# a slice) in Z^n, n = 2, 3, 4.
# Algorithm: growth of rooted link-connected plaquette sets (lattice animals), one neighbouring plaquette at a time, with
# the prune  c1(S) <= 4 (Nmax - |S|),  c1 = number of links of count exactly 1: each further plaquette covers at most four
# links.  Completeness: a valid set (no link of count 1) of size <= Nmax that is not link-connected to the root would contain
# a valid set disjoint from the root's links, i.e. a valid set by itself, whose size is >= 6 by part A (checked first, by
# connected growth from a fixed plaquette); every connected valid set is reached through prefixes that satisfy the prune.
import itertools, sys
from collections import Counter

def unit(n, i): return tuple(1 if k == i else 0 for k in range(n))
def add(x, y): return tuple(a + b for a, b in zip(x, y))
def sub(x, y): return tuple(a - b for a, b in zip(x, y))

class Z:
    def __init__(self, n):
        self.n = n; self._E = {}; self._P = {}
    def E(self, p):                      # links of plaquette p = (x, i, j), i < j; link = (x, i) from x to x + e_i
        if p not in self._E:
            x, i, j = p
            self._E[p] = frozenset([(x, i), (add(x, unit(self.n, i)), j), (add(x, unit(self.n, j)), i), (x, j)])
        return self._E[p]
    def P(self, l):                      # plaquettes containing link l
        if l not in self._P:
            x, i = l; out = []
            for j in range(self.n):
                if j == i: continue
                for base in (x, sub(x, unit(self.n, j))):
                    out.append((base, min(i, j), max(i, j)))
            self._P[l] = out
        return self._P[l]
    def nbrs(self, p):
        return {r for l in self.E(p) for r in self.P(l) if r != p}

def counts(Z_, plaqs, extra=()):
    c = Counter()
    for p in plaqs:
        for l in Z_.E(p): c[l] += 1
    for l in extra: c[l] += 1
    return c

def is_unit_cube(Z_, Y):
    if len(Y) != 6: return False
    V = set()
    for p in Y:
        for (x, i) in Z_.E(p): V.add(x); V.add(add(x, unit(Z_.n, i)))
    if len(V) != 8: return False
    spans = sorted(max(v[k] for v in V) - min(v[k] for v in V) for k in range(Z_.n))
    if spans != [0] * (Z_.n - 3) + [1, 1, 1]: return False
    c = counts(Z_, Y)
    return all(v == 2 for v in c.values()) and len(c) == 12

def grow(Z_, start, root_nbrs, base_count, Nmax, forbidden=frozenset()):
    """all sets X (|X| <= Nmax), X u {root} connected in the plaquette adjacency graph where the root is adjacent to
    root_nbrs, plaquettes not in forbidden, with no link of count 1 in base_count + cover(X). 'start' is an initial set."""
    found = []; seen = set([frozenset(start)])
    frontier = [frozenset(start)]
    def c1(S):
        c = Counter(base_count)
        for p in S:
            for l in Z_.E(p): c[l] += 1
        return sum(1 for v in c.values() if v == 1)
    k0 = len(start)
    for S in frontier:
        if c1(S) == 0: found.append(S)
    while frontier:
        new = []
        for S in frontier:
            if len(S) >= Nmax: continue
            cand = set(root_nbrs)
            for p in S: cand |= Z_.nbrs(p)
            for p in cand:
                if p in S or p in forbidden: continue
                T = S | {p}
                if T in seen: continue
                seen.add(T)
                cc = c1(T)
                if cc > 4 * (Nmax - len(T)): continue
                if cc == 0: found.append(T)
                new.append(T)
        frontier = new
    return found

def partA(n, Nmax):
    Z_ = Z(n); o = (0,) * n; p0 = (o, 0, 1)
    found = grow(Z_, [p0], [], Counter(), Nmax)
    sizes = Counter(len(S) for S in found)
    cubes = sum(1 for S in found if is_unit_cube(Z_, S))
    print(f'A  Z^{n}: valid connected sets containing p0 with <= {Nmax} plaquettes: sizes {dict(sizes)}; unit cubes {cubes}'
          f' (expected {2*(n-2) if n >= 3 else 0})')

def near(n, R):
    out = []
    for x in itertools.product(range(-R, R + 1), repeat=n):
        for i, j in itertools.combinations(range(n), 2):
            out.append((x, i, j))
    return out

def partB(n, R=3):
    Z_ = Z(n); o = (0,) * n; q = (o, 0, 1); F = list(Z_.E(q))
    tally = Counter(); bad = []
    for qp in near(n, R):
        base = counts(Z_, [qp], extra=F)
        rootn = {r for l in set(F) | set(Z_.E(qp)) for r in Z_.P(l)}
        for X in grow(Z_, [], rootn, base, 4, forbidden=frozenset([qp])):
            if qp == q and not X: tally['A: q\'=q, X=0'] += 1
            elif qp != q and q not in X and len(X) == 4 and is_unit_cube(Z_, set(X) | {q, qp}): tally['B: cube'] += 1
            else: bad.append((qp, sorted(X)))
    print(f'B  Z^{n} (F = bd q, |X| <= 4): {dict(tally)}; other: {len(bad)}  (expected cube configurations {5*2*(n-2) if n>=3 else 0})')
    for b in bad[:3]: print('     unexpected', b)

def partC(n, R=3):
    Z_ = Z(n); o = (0,) * n; q = (o, 0, 1)
    for c in Z_.nbrs(q):
        F = Z_.E(q) ^ Z_.E(c)
        tally = Counter(); bad = []
        for qp in near(n, R):
            base = counts(Z_, [qp], extra=F)
            rootn = {r for l in set(F) | set(Z_.E(qp)) for r in Z_.P(l)}
            for X in grow(Z_, [], rootn, base, 2, forbidden=frozenset([qp])):
                if qp in (q, c) and set(X) == {q, c} - {qp}: tally['main'] += 1
                else: bad.append((qp, sorted(X)))
        if bad or tally['main'] != 2:
            print(f'C  Z^{n}, c={c}: {dict(tally)} other {len(bad)}', bad[:3])
    print(f'C  Z^{n} (F = bd q xor bd c, all {len(Z_.nbrs(q))} neighbours c, |X| <= 2): only the two main configurations each: checked')

if __name__ == '__main__':
    import time
    for n, Nmax in ((2, 7), (3, 7), (4, 6)):
        t = time.time(); partA(n, Nmax); print(f'   [{time.time()-t:.1f}s]'); sys.stdout.flush()
    for n in (2, 3, 4):
        t = time.time(); partB(n); print(f'   [{time.time()-t:.1f}s]'); sys.stdout.flush()
    for n in (2, 3, 4):
        t = time.time(); partC(n); print(f'   [{time.time()-t:.1f}s]'); sys.stdout.flush()

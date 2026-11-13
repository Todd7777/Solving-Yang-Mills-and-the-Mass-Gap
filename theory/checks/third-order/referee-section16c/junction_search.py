# Referee (Section 16, third pass): candidate junctions containing p = (bd q, s-1), from the definitions alone.
# Usage: python3 junction_search.py FRONTS_DUMP D n SIZES [SIZES ...]
#   FRONTS_DUMP: fronts through the origin of Z^D (output of fronts_redelmeier.c, "D K dumpfile 4");
#   n: the number of plaquettes |J|_C;  SIZES: the sizes of the further front-pairs, e.g. 4,8 or 4,4,4 (several classes allowed).
"""A candidate junction at slice s is a pair (P, Phi): P a set of plaquettes of the slice Z^D (D = d-1), Phi a set of front-pairs (G, layer),
layer 0 = s-1 (the layer of p), 1 = s, with (Definition mg:def:junction, Lemma mg:lem:junction):
  (1) every G a front (connected, all degrees >= 2) -- taken from the Redelmeier list;
  (2) distinct fronts of one layer vertex-disjoint (they are distinct components of T_layer);
  (3) n_J(e) = [e in E_0] + [e in E_1] + #{plaquettes of P containing e} >= 2 for every link e of E(J) (no support);
  (4) every front-pair shares a vertex with E(P) (it is touched by a marker; no supports);
  (5) H_s connected: markers = link-components of P, front-pairs; adjacency = common vertex.
These are necessary conditions, so the counts bound the numbers of junctions from above.

Search (written independently of three_front_count.py / next_order_count.py, which decide all elements through a link at once):
one element per step.  At a state S (a subset of the final J, starting from {p}): if some link has coverage 1, take the smallest
such link l; J contains a further element through l (a plaquette through l, or a front through l in a layer having no front
through l), and the search branches over all of them.  If no link has coverage 1 and S is not the final J, J has an element
sharing a vertex with an element of S (the elements of J are connected by common vertices, by (5)); the search branches over all
plaquettes and fronts of the remaining sizes at the vertices of S.  Visited states are memoised; results are deduplicated as sets.
Prune: a link of coverage 1 needs a new element through it, and the remaining elements have at most 4*(plaquettes left) + (sum of
the remaining front sizes) links."""
import sys, itertools, time
from collections import Counter, defaultdict

def load_fronts(path, D, sizes):
    out = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            t = list(map(int, line.split()))
            k = t[0]
            if k not in sizes: continue
            links = []
            for j in range(k):
                c = t[1 + j*(D+1): 1 + (j+1)*(D+1)]
                links.append((tuple(c[:D]), c[D]))
            out[k].append(frozenset(links))
    return out

def main():
    path, D, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    classes = [tuple(sorted(int(x) for x in s.split(','))) for s in sys.argv[4:]]
    E = [tuple(1 if k == i else 0 for k in range(D)) for i in range(D)]
    add = lambda x, y: tuple(a + b for a, b in zip(x, y))
    sub = lambda x, y: tuple(a - b for a, b in zip(x, y))
    O = (0,) * D
    def ends(l): return (l[0], add(l[0], E[l[1]]))
    def plinks(c):
        x, i, j = c
        return ((x, i), (add(x, E[i]), j), (add(x, E[j]), i), (x, j))
    def pverts(c):
        x, i, j = c
        return (x, add(x, E[i]), add(x, E[j]), add(add(x, E[i]), E[j]))
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
            for a in (0, 1):
                for b in (0, 1):
                    y = v
                    if a: y = sub(y, E[i])
                    if b: y = sub(y, E[j])
                    out.append((y, i, j))
        return out
    allsizes = set(k for c in classes for k in c)
    T = load_fronts(path, D, allsizes)
    for k in allsizes:
        assert all(len(G) == k for G in T[k])
    def translate(G, x): return frozenset((add(l[0], x), l[1]) for l in G)
    VC = {}
    def fverts(G):
        r = VC.get(G)
        if r is None:
            r = frozenset(v for l in G for v in ends(l)); VC[G] = r
        return r
    # templates through the origin vertex (each front through O once) and through the link (O, i)
    FV = {k: T[k] for k in allsizes}
    FL = {k: {i: set() for i in range(D)} for k in allsizes}
    for k in allsizes:
        for G in T[k]:
            for (x, i) in G:
                FL[k][i].add(translate(G, sub(O, x)))
    FL = {k: {i: list(s) for i, s in FL[k].items()} for k in allsizes}
    for k in allsizes:
        print(f"  fronts of size {k}: {len(FV[k])} through a vertex, {len(FL[k][0])} through a link of direction 0")
    q = (O, 0, 1); dq = frozenset(plinks(q)); p0 = (dq, 0)

    for SIZES in classes:
        t0 = time.time()
        need = Counter(SIZES)
        found = set(); visited = set()
        def remaining(F):
            r = Counter(need)
            for G, lam in F:
                if (G, lam) != p0: r[len(G)] -= 1
            return r
        def ok_front(F, G, lam):
            if (G, lam) in F: return False
            gv = fverts(G)
            return all(not (gv & fverts(H)) for H, mu in F if mu == lam)
        def check(P, F):
            """the definitions, checked from scratch on the final candidate"""
            A = set(); B = set()
            for G, lam in F: (A if lam == 0 else B).update(G)
            pi = Counter(l for c in P for l in plinks(c))
            EJ = A | B | set(pi)
            if any((l in A) + (l in B) + pi[l] < 2 for l in EJ): return False
            for (G, lam), (H, mu) in itertools.combinations(F, 2):
                if lam == mu and fverts(G) & fverts(H): return False
            VP = set(v for c in P for v in pverts(c))
            if not all(fverts(G) & VP for G, lam in F): return False
            Pl = list(P); par = list(range(len(Pl)))
            def fd(a):
                while par[a] != a: a = par[a]
                return a
            for a, b in itertools.combinations(range(len(Pl)), 2):
                if set(plinks(Pl[a])) & set(plinks(Pl[b])): par[fd(a)] = fd(b)
            comps = defaultdict(set)
            for a in range(len(Pl)): comps[fd(a)].update(pverts(Pl[a]))
            nodes = list(comps.values()) + [set(fverts(G)) for G, lam in F]
            par2 = list(range(len(nodes)))
            def f2(a):
                while par2[a] != a: a = par2[a]
                return a
            for a, b in itertools.combinations(range(len(nodes)), 2):
                if nodes[a] & nodes[b]: par2[f2(a)] = f2(b)
            return len({f2(a) for a in range(len(nodes))}) == 1
        def rec(P, F, cnt):
            key = (P, F)
            if key in visited: return
            visited.add(key)
            rem = remaining(F)
            if min(rem.values(), default=0) < 0: return
            rp = n - len(P); rf = sum(rem.values())
            defi = [l for l, c in cnt.items() if c == 1]
            if len(defi) > 4*rp + sum(k*m for k, m in rem.items()): return
            if not defi:
                if rp == 0 and rf == 0:
                    if check(P, F): found.add(key)
                    return
                V = set()
                for c in P: V.update(pverts(c))
                for G, lam in F: V.update(fverts(G))
                for v in V:
                    if rp > 0:
                        for c in plaqs_at(v):
                            if c in P: continue
                            cnt2 = Counter(cnt)
                            for l in plinks(c): cnt2[l] += 1
                            rec(P | {c}, F, cnt2)
                    for k, m in rem.items():
                        if m <= 0: continue
                        for G0 in FV[k]:
                            G = translate(G0, v)
                            for lam in (0, 1):
                                if ok_front(F, G, lam):
                                    cnt2 = Counter(cnt)
                                    for l in G: cnt2[l] += 1
                                    rec(P, F | {(G, lam)}, cnt2)
                return
            l = min(defi)
            if rp > 0:
                for c in plaqs_through(l):
                    if c in P: continue
                    cnt2 = Counter(cnt)
                    for m in plinks(c): cnt2[m] += 1
                    rec(P | {c}, F, cnt2)
            for k, m in rem.items():
                if m <= 0: continue
                for G0 in FL[k][l[1]]:
                    G = translate(G0, l[0])
                    for lam in (0, 1):
                        if ok_front(F, G, lam):
                            cnt2 = Counter(cnt)
                            for mm in G: cnt2[mm] += 1
                            rec(P, F | {(G, lam)}, cnt2)
        sys.setrecursionlimit(10000)
        rec(frozenset(), frozenset([p0]), Counter(dq))
        stat = Counter()
        for P, F in found:
            others = tuple(sorted((len(G), 'same' if lam == 0 else 'other') for G, lam in F if (G, lam) != p0))
            stat[others] += 1
        print(f"D={D} n={n} sizes={SIZES}: {len(found)} junctions  {dict(sorted(stat.items()))}  "
              f"({len(visited)} states, {time.time()-t0:.1f} s)")
        sys.stdout.flush()

if __name__ == "__main__":
    main()

# Referee check for Part IV, Section 16, Lemma mg:lem:threelow (added while Section 16 was under review): the junctions without support with
# exactly three front-pairs, one of them the given p = (bd q0, lambda), and
#   (i)  three four-edge fronts and |P| = 3 plaquettes            (claimed: N_3 = 3 nu_2 + 12(d-3) = 48 (d=3), 876 (d=4));
#   (ii) fronts of 4, 4 and 6 edges (p has 4) and |P| = 2         (claimed: N_2 = K_S^2/2 + K_S J_v + 8(d-3) = 24 (d=3), 320 (d=4)).
# Written from Definitions mg:def:fronts, mg:def:junction and Lemma mg:lem:junction, independently of three_front_count.py:
#   n(e) := sum over the front-pairs of [e in front] + #{r in P: e in r} >= 2 for every e in E(J);  fronts of one layer vertex-disjoint;
#   each front-pair shares a vertex with E(P);  H_s connected (clusters = link-components of P, front-pairs; adjacency = a common vertex).
# Fronts with four edges are plaquette boundaries, with six edges six-cycles (spatial lattice Z^m, m = d-1).
# Search over states (P, set of front-pairs) with deduplication: resolve a deficient edge by adding any plaquette or allowed front-pair through it;
# without deficient edges, add any plaquette or allowed front-pair with a vertex in V(J) (complete: see ref_junctions.py). Budgets prune.
# Usage: python3 ref_three.py m case       (case = i, ii, or iii: fronts of 4, 4, 7 edges, two plaquettes)
import sys, itertools
from collections import Counter
m = int(sys.argv[1]); CASE = sys.argv[2]
NP = 3 if CASE == 'i' else 2
SIZES = {'i': Counter({4: 3}), 'ii': Counter({4: 2, 6: 1}), 'iii': Counter({4: 2, 7: 1})}[CASE]   # front sizes, p included
# case iii (not in the lemma): two plaquettes and fronts of 4, 4 and 7 edges, the pairs with (n, E) = (2, 7/2) in the proof of Lemma mg:lem:jsums5
U = [tuple(1 if k == i else 0 for k in range(m)) for i in range(m)]
def add(x, y): return tuple(a + b for a, b in zip(x, y))
def sub(x, y): return tuple(a - b for a, b in zip(x, y))
def edge(x, y): return (x, y) if x < y else (y, x)
def links(p):
    x, i, j = p
    a, b, c, d = x, add(x, U[i]), add(add(x, U[i]), U[j]), add(x, U[j])
    return frozenset((edge(a, b), edge(b, c), edge(c, d), edge(d, a)))
def plaqs_at(v):
    out = []
    for i in range(m):
        for j in range(i + 1, m):
            for a in (0, 1):
                for b in (0, 1):
                    out.append((sub(sub(v, tuple(a * c for c in U[i])), tuple(b * c for c in U[j])), i, j))
    return out
def plaqs_through(e):
    x, y = e
    return [p for p in plaqs_at(x) if e in links(p)]
def nbrs(v): return [add(v, U[k]) for k in range(m)] + [sub(v, U[k]) for k in range(m)]
def verts(E): return frozenset(v for e in E for v in e)
_six = {}
def six_through(e):
    """all six-cycles (as edge sets) containing the edge e"""
    if e in _six: return _six[e]
    x, y = e; out = set()
    def dfs(path):
        if len(path) == 6:
            if path[-1] == y: out.add(frozenset(edge(path[k], path[k + 1]) for k in range(5)) | {e})
            return
        for w in nbrs(path[-1]):
            if w in path: continue
            if len(path) < 5 and w == y: continue
            dfs(path + [w])
    dfs([x])
    _six[e] = [c for c in out if len(c) == 6]
    return _six[e]
_seven = {}
def seven_through(e):
    """all seven-edge fronts containing e: for sides >= 8 these are bd r u bd c, r, c plaquettes sharing a link (Remark mg:rem:betafour(iii))"""
    if e in _seven: return _seven[e]
    out = set()
    x, y = e
    near = set()
    for v in (x, y):
        for w in [v] + nbrs(v):
            near.update(plaqs_at(w))
    for r in near:
        for c in near:
            if r < c and len(links(r) & links(c)) == 1:
                F = links(r) | links(c)
                if e in F: out.add(F)
    _seven[e] = list(out)
    return _seven[e]
def seven_at(v):
    out = set()
    for w in nbrs(v): out.update(seven_through(edge(v, w)))
    return out
def six_at(v):
    out = set()
    for w in nbrs(v): out.update(six_through(edge(v, w)))
    return out

q0 = (tuple([0] * m), 0, 1)
P0 = (links(q0), 1)                                   # p = (bd q0, layer 1)
visited = set(); found = set()

def budget_ok(P, FR):
    sizes = Counter(len(F) for F, l in FR)
    return len(P) <= NP and all(sizes[k] <= SIZES[k] for k in sizes) and len(FR) <= 3
def front_ok(F, l, FR):
    if (F, l) in FR: return False
    VF = verts(F)
    return all(not (VF & verts(G)) for G, l2 in FR if l2 == l)     # same layer: vertex-disjoint (hence also distinct)
def ncount(P, FR):
    n = Counter()
    for r in P:
        for e in links(r): n[e] += 1
    for F, l in FR:
        for e in F: n[e] += 1
    return n
def connected(P, FR):
    nodes = []
    # clusters
    P = list(P); par = list(range(len(P)))
    def f(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    for a, b in itertools.combinations(range(len(P)), 2):
        if links(P[a]) & links(P[b]): par[f(a)] = f(b)
    cl = {}
    for a in range(len(P)): cl.setdefault(f(a), set()).update(verts(links(P[a])))
    nodes = list(cl.values()) + [set(verts(F)) for F, l in FR]
    par2 = list(range(len(nodes)))
    def g(a):
        while par2[a] != a: par2[a] = par2[par2[a]]; a = par2[a]
        return a
    for a, b in itertools.combinations(range(len(nodes)), 2):
        if nodes[a] & nodes[b]: par2[g(a)] = g(b)
    return len({g(a) for a in range(len(nodes))}) == 1

def options_through(e, P, FR):
    out = []
    if len(P) < NP:
        out += [('P', r) for r in plaqs_through(e) if r not in P]
    sizes = Counter(len(F) for F, l in FR)
    if len(FR) < 3:
        if sizes[4] < SIZES[4]:
            for r in plaqs_through(e):
                F = links(r)
                for l in (0, 1):
                    if front_ok(F, l, FR): out.append(('F', (F, l)))
        if sizes[6] < SIZES[6]:
            for F in six_through(e):
                for l in (0, 1):
                    if front_ok(F, l, FR): out.append(('F', (F, l)))
        if sizes[7] < SIZES[7]:
            for F in seven_through(e):
                for l in (0, 1):
                    if front_ok(F, l, FR): out.append(('F', (F, l)))
    return out

def rec(P, FR):
    key = (P, FR)
    if key in visited: return
    visited.add(key)
    if not budget_ok(P, FR): return
    n = ncount(P, FR)
    defi = {e: 2 - c for e, c in n.items() if c < 2}
    sizes = Counter(len(F) for F, l in FR)
    room = 4 * (NP - len(P)) + sum((SIZES[k] - sizes[k]) * k for k in SIZES)
    if sum(defi.values()) > room: return
    if defi:
        best = None
        for e in defi:
            o = options_through(e, P, FR)
            if best is None or len(o) < len(best):
                best = o
                if not o: return
        for kind, obj in best:
            if kind == 'P': rec(P | {obj}, FR)
            else: rec(P, FR | {obj})
        return
    if len(P) == NP and sizes == SIZES and len(FR) == 3:
        VP = verts(frozenset(e for r in P for e in links(r)))
        if all(verts(F) & VP for F, l in FR) and connected(P, FR):
            found.add(key)
    # seeds at the vertices of J
    VJ = verts(frozenset(e for r in P for e in links(r))) | frozenset(v for F, l in FR for v in verts(F))
    if len(P) < NP:
        for v in VJ:
            for r in plaqs_at(v):
                if r not in P: rec(P | {r}, FR)
    if len(FR) < 3:
        for v in VJ:
            if sizes[4] < SIZES[4]:
                for r in plaqs_at(v):
                    for l in (0, 1):
                        if front_ok(links(r), l, FR): rec(P, FR | {(links(r), l)})
            if sizes[6] < SIZES[6]:
                for F in six_at(v):
                    for l in (0, 1):
                        if front_ok(F, l, FR): rec(P, FR | {(F, l)})
            if sizes[7] < SIZES[7]:
                for F in seven_at(v):
                    for l in (0, 1):
                        if front_ok(F, l, FR): rec(P, FR | {(F, l)})

sys.setrecursionlimit(100000)
rec(frozenset(), frozenset([P0]))
# classify
def classify(P, FR):
    others = [(F, l) for F, l in FR if (F, l) != P0]
    if CASE == 'i':
        caps = all(any(links(r) == F for r in P) for F, l in FR)
        return 'caps' if caps else 'other (cube?)'
    big = [l for F, l in others if len(F) > 4]
    return 'big front in the layer of p' if big and big[0] == 1 else 'big front in the other layer'
cls = Counter(classify(P, FR) for P, FR in found)
print(f"d={m+1}, case ({CASE}): {len(found)} junctions given p; classes {dict(cls)}; states visited {len(visited)}")

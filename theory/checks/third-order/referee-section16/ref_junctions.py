# Referee check for Part IV, Section 16, Remark mg:rem:practical(ii): Z_{n,4}, the number of pairs (J, p') with J a junction without
# support at a slice s, Phi(J) = {p, p'}, |p| = |p'| = 4, |J|_C = n, for a given front-pair p = (bd q0, lambda) (Lemma mg:lem:twofrontcount).
# Written from the definitions (Definition mg:def:fronts, Definition mg:def:junction, Lemma mg:lem:junction), independently of
# junction_count.py. Spatial lattice Z^m, m = d - 1 (the counts of a periodic lattice with large sides are the same).
#
# A pair (J, p') is given by the set P of the plaquettes of the clusters of J (the clusters are the link-components of P) and
# p' = (bd q', lambda'), q' a plaquette (fronts with four edges in Z^m are plaquette boundaries). The local conditions used:
#   (a) n(e) := [e in bd q0] + [e in bd q'] + #{p in P : e in p} >= 2 for every e in E(J) = bd q0 u bd q' u E(P)
#       (Lemma mg:lem:junction(iii); when lambda' = lambda, E_lambda(J) = bd q0 u bd q' with the two fronts disjoint, so n(e) has the same form);
#   (b) if lambda' = lambda, the fronts are distinct components of T_lambda, so bd q0 and bd q' have no common vertex;
#   (c) each front-pair is touched by a marker: V(bd q0) and V(bd q') meet V(E(P));
#   (d) J is connected: the graph on {p, p'} u {clusters of P}, adjacent when the edge sets share a vertex, is connected;
#   (e) (p', lambda') != (p, lambda), which (b) implies for lambda' = lambda.
# So each admissible (P, q') counts 1 (lambda' the other layer) plus 1 if bd q0 and bd q' are vertex-disjoint (lambda' = lambda).
#
# Search for P, by decisions: an edge e of E(J) is deficient if n(e) < 2; take the deficient edge with fewest undecided plaquettes through it
# and branch over the subsets T of those with |T| >= 2 - n(e) (the rest excluded). Without deficient edges, branch on the smallest undecided
# plaquette with a vertex in V(J) (include / exclude); record at the end of that chain. Completeness: if P is admissible and IN is a proper
# subset with no deficient edge, P \ IN contains a plaquette with a vertex in V(J_IN) (same cluster as a plaquette of IN, or, by (d), a
# cluster of P \ IN shares a vertex with a front or a cluster of IN). Pruning by two lower bounds on the number of plaquettes still needed:
# ceil(total deficit / 4), and, if the current J is disconnected, the L-infinity distance from a component to the others.
# Usage: python3 ref_junctions.py m NMAX [list] [nosym] [nofrontadj]
import sys, itertools
from collections import Counter

m = int(sys.argv[1]); NMAX = int(sys.argv[2]); LIST = 'list' in sys.argv[3:]
U = [tuple(1 if k == i else 0 for k in range(m)) for i in range(m)]
def add(x, y): return tuple(a + b for a, b in zip(x, y))
def sub(x, y): return tuple(a - b for a, b in zip(x, y))
def links(p):
    x, i, j = p
    return ((x, i), (add(x, U[i]), j), (add(x, U[j]), i), (x, j))
def through(e):
    x, k = e; out = []
    for t in range(m):
        if t == k: continue
        i, j = min(k, t), max(k, t)
        out.append((x, i, j)); out.append((sub(x, U[t]), i, j))
    return out
def verts_of_edge(e): x, k = e; return (x, add(x, U[k]))
def pverts(p):
    x, i, j = p
    return (x, add(x, U[i]), add(x, U[j]), add(add(x, U[i]), U[j]))
def plaqs_at(v):
    out = []
    for i in range(m):
        for j in range(i + 1, m):
            for a in (0, 1):
                for b in (0, 1):
                    out.append((sub(sub(v, tuple(a * c for c in U[i])), tuple(b * c for c in U[j])), i, j))
    return out
def linf(A, B): return min(max(abs(a - b) for a, b in zip(x, y)) for x in A for y in B)

q0 = (tuple([0] * m), 0, 1)
PHI0 = frozenset(links(q0)); V0 = frozenset(pverts(q0))

NOFRONTADJ = 'nofrontadj' in sys.argv     # variant: the two front-pairs (nodes 0 and 1) are never joined directly, as in junction_count.py

def components(nodes):
    """nodes: list of vertex sets; return list of lists of node indices, grouped by shared vertices"""
    parent = list(range(len(nodes)))
    def find(i):
        while parent[i] != i: parent[i] = parent[parent[i]]; i = parent[i]
        return i
    if NOFRONTADJ:
        for a in range(len(nodes)):
            for b in range(a + 1, len(nodes)):
                if a <= 1 and b <= 1: continue
                if nodes[a] & nodes[b]: parent[find(a)] = find(b)
        groups = {}
        for i in range(len(nodes)): groups.setdefault(find(i), []).append(i)
        return list(groups.values())
    owner = {}
    for idx, vs in enumerate(nodes):
        for v in vs:
            if v in owner: parent[find(idx)] = find(owner[v])
            else: owner[v] = idx
    groups = {}
    for i in range(len(nodes)): groups.setdefault(find(i), []).append(i)
    return list(groups.values())

def clusters(P):
    P = list(P); idx = {p: k for k, p in enumerate(P)}
    parent = list(range(len(P)))
    def find(i):
        while parent[i] != i: parent[i] = parent[parent[i]]; i = parent[i]
        return i
    owner = {}
    for k, p in enumerate(P):
        for e in links(p):
            if e in owner: parent[find(k)] = find(owner[e])
            else: owner[e] = k
    groups = {}
    for k in range(len(P)): groups.setdefault(find(k), []).append(P[k])
    return list(groups.values())

def jnodes(P, V1):
    nodes = [V0, V1]
    for C in clusters(P): nodes.append(frozenset(v for p in C for v in pverts(p)))
    return nodes

def search(q1, n, results):
    PHI1 = frozenset(links(q1)); V1 = frozenset(pverts(q1))
    F = Counter(PHI0) + Counter(PHI1)
    stats = [0]
    def rec(IN, OUT, pi):
        stats[0] += 1
        # deficits
        edges = set(PHI0) | set(PHI1) | set(pi)
        defi = {e: 2 - F[e] - pi[e] for e in edges if F[e] + pi[e] < 2}
        lb = (sum(defi.values()) + 3) // 4
        nodes = jnodes(IN, V1); comps = components(nodes)
        if len(comps) > 1:
            vs = [frozenset().union(*(nodes[i] for i in c)) for c in comps]
            lbc = max(linf(vs[a], frozenset().union(*(vs[b] for b in range(len(vs)) if b != a))) for a in range(len(vs)))
            lb = max(lb, lbc)
        if len(IN) + lb > n: return
        if defi:
            best = None
            for e, need in defi.items():
                c = [p for p in through(e) if p not in IN and p not in OUT]
                if best is None or len(c) - need < len(best[1]) - best[2]:
                    best = (e, c, need)
            e, cands, need = best
            room = n - len(IN)
            for r in range(need, min(len(cands), room) + 1):
                for T in itertools.combinations(cands, r):
                    pi2 = Counter(pi)
                    for p in T:
                        for ee in links(p): pi2[ee] += 1
                    rec(IN | frozenset(T), OUT | (frozenset(cands) - frozenset(T)), pi2)
            return
        if len(IN) < n:
            VJ = set(V0) | set(V1) | {v for p in IN for v in pverts(p)}
            und = sorted({p for v in VJ for p in plaqs_at(v) if p not in IN and p not in OUT})
            if und:
                c = und[0]; pi2 = Counter(pi)
                for ee in links(c): pi2[ee] += 1
                rec(IN | {c}, OUT, pi2)
                rec(IN, OUT | {c}, pi)
                return
        if len(IN) != n: return
        # (c) touching and (d) connectivity
        VP = {v for p in IN for v in pverts(p)}
        if not (V0 & VP) or not (V1 & VP): return
        if len(components(jnodes(IN, V1))) != 1: return
        results.append((q1, IN))
    rec(frozenset(), frozenset(), Counter())
    return stats[0]

def relation(q1):
    if q1 == q0: return 'same plaquette'
    if len(set(links(q1)) & PHI0) == 1: return 'share a link'
    if set(pverts(q1)) & V0: return 'share a vertex only'
    return 'vertex-disjoint'

# The lattice symmetries fixing q0 (vertices {0,1}^2 x {0}^(m-2)): the symmetries of the unit square in the first two coordinates, times
# permutations and reflections x_k -> -x_k of the others. They map admissible (P, q') to admissible ones, so it suffices to search one q' per
# orbit and to weight it by the orbit size (option 'nosym' searches every q').
def symmetries():
    out = []
    for perm in itertools.permutations(range(m)):
        if set(perm[:2]) != {0, 1}: continue
        for signs in itertools.product((1, -1), repeat=m):
            def g(v, perm=perm, signs=signs):
                w = [v[perm[k]] for k in range(m)]
                return tuple((w[k] if signs[k] == 1 else (1 - w[k] if k < 2 else -w[k])) for k in range(m))
            out.append(g)
    return out
SYM = symmetries()
def plaq_from_verts(vs):
    x = tuple(min(v[k] for v in vs) for k in range(m))
    dirs = [k for k in range(m) if max(v[k] for v in vs) != x[k]]
    return (x, dirs[0], dirs[1])
def gplaq(g, p): return plaq_from_verts([g(v) for v in pverts(p)])
assert all(gplaq(g, q0) == q0 for g in SYM)
R = NMAX + 2
allq = []
for x in itertools.product(range(-R, R + 1), repeat=m):
    for i in range(m):
        for j in range(i + 1, m):
            q1 = (x, i, j)
            if linf(V0, pverts(q1)) <= NMAX: allq.append(q1)
NOSYM = 'nosym' in sys.argv
reps = []
if NOSYM: reps = [(q1, 1) for q1 in allq]
else:
    seen = set()
    for q1 in allq:
        if q1 in seen: continue
        orb = {gplaq(g, q1) for g in SYM}
        seen |= orb
        reps.append((min(orb), len(orb)))
print(f"d={m+1}: {len(allq)} candidates q', {len(reps)} searched")
for n in range(1, NMAX + 1):
    res = []; nodes = 0; weights = {}
    for q1, w in reps:
        k0 = len(res); nodes += search(q1, n, res)
        for t in range(k0, len(res)): weights[t] = w
    other = sum(weights[t] for t in range(len(res))); same = sum(weights[t] for t, (q1, P) in enumerate(res) if relation(q1) == 'vertex-disjoint')
    rel = Counter(); 
    for t, (q1, P) in enumerate(res):
        rel[(relation(q1), 'other layer')] += weights[t]
        if relation(q1) == 'vertex-disjoint': rel[(relation(q1), 'same layer')] += weights[t]
    print(f"d={m+1}, n={n}: Z_n,4 = {other + same}  (lambda' other layer: {other}, same layer: {same}); by relation of q' to q0: {dict(rel)}; nodes {nodes}")
    sys.stdout.flush()
    if LIST:
        def is_cube(Y):
            if len(Y) != 6: return False
            E = Counter(e for p in Y for e in links(p))
            return len(E) == 12 and all(c == 2 for c in E.values())
        def describe(q1, P):
            P = set(P); tags = []
            if is_cube(P | {q0, q1}) and q0 not in P and q1 not in P and q1 != q0: tags.append("P = Gamma minus {q0, q'}")
            if q1 != q0 and q0 in P and is_cube(P | {q1}): tags.append("P = Gamma minus {q'} (q0 added)")
            if q1 != q0 and q1 in P and is_cube(P | {q0}): tags.append("P = Gamma minus {q0} (q' added)")
            if q1 == q0 and is_cube(P | {q0}) and q0 not in P: tags.append("tube through the five other faces of a cube")
            return ', '.join(tags) or 'other'
        cls = Counter()
        for t, (q1, P) in enumerate(res):
            lay = 2 if relation(q1) == 'vertex-disjoint' else 1
            cls[(describe(q1, P), relation(q1))] += weights[t] * lay
        print(f"     classes (weighted, both layers): {dict(cls)}")
        if n <= 5:
            for t, (q1, P) in enumerate(res):
                print(f"     weight {weights[t]}: q'={q1} ({relation(q1)}), P={sorted(P)}: {describe(q1, P)}")

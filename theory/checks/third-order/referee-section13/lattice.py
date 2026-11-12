# Written by the independent referee of Section 13 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs. Code written from the definitions alone (no lemma of the paper is used).
# Independent referee code for Section 13 (Lemma 13.1).  Written from the definitions only:
#   front  = non-empty connected set of edges in which every vertex has degree >= 2 (Section 8);
#   junction J at a slice (no support) = clusters (link-components of a set P of plaquettes of the slice) and front-pairs (phi, layer),
#   layer in {-1 (= s-1), +1 (= s)}, such that
#     (J1) fronts of one layer are distinct components of T_layer: pairwise vertex-disjoint;
#     (J2) every front-pair shares a vertex with some cluster (it is "touched by a marker"; no supports here);
#     (J3) J is connected in H_s (nodes: clusters and front-pairs; adjacent when their edge sets share a vertex);
#     (J4) n(e) = [e in A] + [e in B] + #{plaquettes of P containing e} >= 2 for every e in E(J)   (no free link, Lemma 9.2(iii)),
#          A, B = unions of the fronts of layer -1, +1.
# Everything is done in Z^n (n = spatial dimension = d-1) near the origin, which is what a torus with all sides >= 8 looks like locally.
import itertools
from collections import defaultdict

def E(u, v):
    return (u, v) if u < v else (v, u)

def add(x, y):
    return tuple(a + b for a, b in zip(x, y))

def unit(n, i, s=1):
    return tuple(s if k == i else 0 for k in range(n))

def nbrs(v, n):
    for i in range(n):
        for s in (1, -1):
            yield add(v, unit(n, i, s))

def plaq(x, i, j, n):
    """plaquette spanned at corner x by e_i, e_j: frozenset of its four edges"""
    a = tuple(x); b = add(a, unit(n, i)); c = add(b, unit(n, j)); dd = add(a, unit(n, j))
    return frozenset([E(a, b), E(b, c), E(c, dd), E(dd, a)])

def verts(edges):
    return frozenset(v for e in edges for v in e)

def edge_dir(e):
    u, v = e
    return [k for k in range(len(u)) if u[k] != v[k]][0]

def plaqs_containing_edge(e, n):
    u, v = e
    i = edge_dir(e)
    out = []
    for j in range(n):
        if j == i: continue
        out.append(plaq(u, i, j, n))
        out.append(plaq(add(u, unit(n, j, -1)), i, j, n))
    return out

def plaqs_at_vertex(v, n):
    out = set()
    for i, j in itertools.combinations(range(n), 2):
        for a in (0, 1):
            for b in (0, 1):
                x = add(add(v, unit(n, i, -a)), unit(n, j, -b))
                out.add(plaq(x, i, j, n))
    return out

def plaqs_touching(p, n):
    """plaquettes sharing at least a vertex with p (including p)"""
    out = set()
    for v in verts(p):
        out |= plaqs_at_vertex(v, n)
    return out

def is_front(F):
    if not F: return False
    deg = defaultdict(int)
    for (u, v) in F:
        deg[u] += 1; deg[v] += 1
    if min(deg.values()) < 2: return False
    # connectivity (edges linked through common endpoints)
    F = list(F); seen = {F[0]}; st = [F[0]]
    inc = defaultdict(list)
    for e in F:
        inc[e[0]].append(e); inc[e[1]].append(e)
    while st:
        e = st.pop()
        for w in e:
            for f in inc[w]:
                if f not in seen: seen.add(f); st.append(f)
    return len(seen) == len(F)

def fronts_through(v0, K, n, R=None):
    """all fronts with at most K edges containing the vertex v0, by Redelmeier enumeration of connected edge sets in the line graph
    restricted to the L1-ball of radius R around v0, pruned by 'number of degree-1 vertices <= 2*(edges still to add)'."""
    if R is None: R = K // 2 + 1
    inball = lambda w: sum(abs(a - b) for a, b in zip(w, v0)) <= R
    def edges_at(w):
        return [E(w, x) for x in nbrs(w, n) if inball(x)]
    root_edges = sorted(edges_at(v0))
    out = set()
    deg = defaultdict(int)
    for idx, root in enumerate(root_edges):
        seen = set(root_edges[:idx]) | {root}
        S = []
        state = {'L': 0}
        def addE(e):
            for w in e:
                deg[w] += 1
                if deg[w] == 1: state['L'] += 1
                elif deg[w] == 2: state['L'] -= 1
        def remE(e):
            for w in e:
                deg[w] -= 1
                if deg[w] == 1: state['L'] += 1
                elif deg[w] == 0: state['L'] -= 1
        def rec(untried):
            untried = list(untried)
            while untried:
                e = untried.pop()
                S.append(e); addE(e)
                if state['L'] == 0 and len(S) >= 4:
                    out.add(frozenset(S))
                if len(S) < K and state['L'] <= 2 * (K - len(S)):
                    added = []
                    for w in e:
                        for f in edges_at(w):
                            if f not in seen:
                                seen.add(f); added.append(f)
                    rec(untried + added)
                    for f in added: seen.discard(f)
                remE(e); S.pop()
        rec([root])
    return out

# junction_enum.py -- written by the independent referee of Section 11 of Part IV (2026-11-11), added to the repository with credit;
# unchanged apart from this header. A compass, not part of the proofs.
# Referee check (independent code) for Lemmas 11.2 and 11.3 of Part IV (Section 11), and for the counts Z_{n,4} of
# Lemma 8.7 (mg:lem:twofrontcount).  Written from the definitions only:
#   * a junction J at slice s of a set X: markers = clusters (link-connected components of the spatial plaquettes of X at
#     time s; no support here), front-pairs (psi, mu) = components psi of T_mu (mu in {s-1, s}) touched by a marker, J a
#     connected component of the graph H_s (common vertex);
#   * the no-free-link condition, locally: n(e) = [e in T_{s-1}] + [e in T_s] + #{plaquettes of J containing e} is 0 or >= 2
#     for every spatial link e at time s (for e in E(J) it is >= 1, hence >= 2);
#   * fronts of a layer have minimum degree >= 2 (temporal links); distinct fronts of one layer are vertex-disjoint.
# The four-edge fronts are ALL 4-cycles of the torus graph (plaquette boundaries, and straight cycles when a side is 4),
# found by brute force over closed walks; no lemma of the paper is used to restrict them.
#
# Part A (Lemma 11.2 and Z_{n,4}): fix p = (phi, s-1), phi a 4-cycle; for every 4-cycle phi' (within graph distance 2*nmax)
#   and layer, enumerate ALL sets Q of n <= nmax plaquettes making {Q, phi, phi'} a junction with exactly these two
#   front-pairs, by a complete backtracking search: while some edge of E(J) has n(e) < 2, branch over the plaquettes containing
#   it; when none is left, record Q and also branch over every plaquette sharing an edge with E(J) (a plaquette sharing no edge
#   with E(J) has four unsatisfied edges, which need four further distinct plaquettes, so it cannot occur for n <= 4 after Q).
# Part B (Lemma 11.3): one plaquette c and front-pairs whose fronts are 4-cycles (any number, any layers): enumerate all
#   admissible junctions; and, for sides >= 6, those with fronts of 4 and 6 edges (6-cycles), as a tabulation.
import itertools, sys
from collections import defaultdict, Counter

class Torus:
    def __init__(self, L):
        self.L = tuple(L); D = self.D = len(L)
        self.verts = list(itertools.product(*[range(l) for l in L]))
        self.vidx = {v: i for i, v in enumerate(self.verts)}
        self.edges = []          # (v, i)
        self.eidx = {}
        for v in self.verts:
            for i in range(D):
                self.eidx[(v, i)] = len(self.edges); self.edges.append((v, i))
        self.ends = [frozenset((self.vidx[v], self.vidx[self.add(v, i)])) for (v, i) in self.edges]
        self.vedges = defaultdict(list)
        for e, en in enumerate(self.ends):
            for x in en: self.vedges[x].append(e)
        self.plaqs = []          # frozenset of 4 edge indices
        self.pkey = []
        for v in self.verts:
            for i in range(D):
                for j in range(i+1, D):
                    E = frozenset([self.eidx[(v, i)], self.eidx[(self.add(v, i), j)], self.eidx[(self.add(v, j), i)], self.eidx[(v, j)]])
                    assert len(E) == 4
                    self.plaqs.append(E); self.pkey.append((v, i, j))
        self.eplaqs = defaultdict(list)
        for p, E in enumerate(self.plaqs):
            for e in E: self.eplaqs[e].append(p)
        self.pverts = [frozenset(x for e in E for x in self.ends[e]) for E in self.plaqs]
        self.plaq_of_edges = {E: p for p, E in enumerate(self.plaqs)}

    def add(self, v, i, s=1):
        w = list(v); w[i] = (w[i] + s) % self.L[i]; return tuple(w)

    def nbrs(self, x):
        return [y for e in self.vedges[x] for y in self.ends[e] if y != x]

    def edge_between(self, x, y):
        for e in self.vedges[x]:
            if y in self.ends[e]: return e
        return None

    def cycles(self, k, through=None):
        """All cycles with k edges (as frozensets of edges) through the given vertices (all if None)."""
        out = set()
        starts = range(len(self.verts)) if through is None else through
        for x0 in starts:
            def walk(path):
                if len(path) == k:
                    if path[0] in self.nbrs(path[-1]):
                        es = frozenset(self.edge_between(path[i], path[(i+1) % k]) for i in range(k))
                        if len(es) == k: out.add(es)
                    return
                for y in self.nbrs(path[-1]):
                    if y in path: continue
                    walk(path + [y])
            walk([x0])
        return out

    def verts_of(self, E):
        return frozenset(x for e in E for x in self.ends[e])

    def dist_sets(self, A, B, cap):
        """graph distance between vertex sets A and B (or cap+1 if larger)."""
        seen = set(A); frontier = set(A); d = 0
        while frontier:
            if frontier & B: return d
            if d >= cap: return cap+1
            nxt = set()
            for x in frontier:
                for y in self.nbrs(x):
                    if y not in seen: seen.add(y); nxt.add(y)
            frontier = nxt; d += 1
        return cap+1

def components_by_links(T, Q):
    """clusters: link-connected components of the plaquette set Q."""
    Q = list(Q); parent = {p: p for p in Q}
    def f(p):
        while parent[p] != p: parent[p] = parent[parent[p]]; p = parent[p]
        return p
    for a, b in itertools.combinations(Q, 2):
        if T.plaqs[a] & T.plaqs[b]: parent[f(a)] = f(b)
    comps = defaultdict(list)
    for p in Q: comps[f(p)].append(p)
    return list(comps.values())

def junction_ok(T, Q, fronts):
    """fronts: list of (edge set, layer). Check: every front touched by a cluster; H_s connected (common vertices; two fronts
    of the same layer are never adjacent: they are vertex-disjoint)."""
    clusters = components_by_links(T, Q)
    nodes = [(T.verts_of(set().union(*[T.plaqs[p] for p in C])), 'C') for C in clusters] + [(T.verts_of(f), lay) for f, lay in fronts]
    VQ = T.verts_of(set().union(*[T.plaqs[p] for p in Q])) if Q else frozenset()
    for f, lay in fronts:
        if not (T.verts_of(f) & VQ): return False
    n = len(nodes); seen = {0}; st = [0]
    while st:
        a = st.pop()
        for b in range(n):
            if b in seen: continue
            if nodes[a][1] != 'C' and nodes[b][1] != 'C' and nodes[a][1] == nodes[b][1]: continue
            if nodes[a][0] & nodes[b][0]: seen.add(b); st.append(b)
    return len(seen) == n

def ncount(T, Q, layer_sets):
    """n(e) for all e in E(J)."""
    cnt = Counter()
    for S in layer_sets:
        for e in S: cnt[e] += 1
    for p in Q:
        for e in T.plaqs[p]: cnt[e] += 1
    return cnt

def enumerate_Q(T, fronts, nmax):
    """All plaquette sets Q, |Q| <= nmax, such that n(e) >= 2 on E(J) for the given fronts (list of (edges, layer))."""
    Tm = set(); Tp = set()
    for f, lay in fronts: (Tm if lay == -1 else Tp).update(f)
    found = set(); visited = set()
    def rec(Q):
        key = frozenset(Q)
        if key in visited: return
        visited.add(key)
        cnt = ncount(T, Q, [Tm, Tp])
        bad = sorted(e for e, c in cnt.items() if c < 2)
        if bad:
            if len(Q) >= nmax: return
            e = bad[0]
            for p in T.eplaqs[e]:
                if p not in Q: rec(Q | {p})
            return
        if Q: found.add(key)
        if len(Q) >= nmax: return
        EJ = set(cnt)
        cand = {p for e in EJ for p in T.eplaqs[e]} - Q
        for p in cand: rec(Q | {p})
    rec(frozenset())
    return found

def partA(L, nmax=3, phis=None, both_layers=True):
    T = Torus(L)
    o = T.vidx[tuple(0 for _ in L)]
    C4 = T.cycles(4)
    print(f"torus {L}: {len(T.plaqs)} plaquettes, {len(C4)} 4-cycles "
          f"({sum(1 for c in C4 if c in T.plaq_of_edges)} plaquette boundaries, {sum(1 for c in C4 if c not in T.plaq_of_edges)} straight)")
    # entries phi: the plaquette boundaries at the origin (each plane) and straight 4-cycles through the origin
    if phis is None:
        phis = [c for c in C4 if o in T.verts_of(c) and (c not in T.plaq_of_edges or T.pkey[T.plaq_of_edges[c]][0] == tuple(0 for _ in L))]
    res = {}
    for phi in phis:
        kind = f"plaquette {T.pkey[T.plaq_of_edges[phi]][1:]}" if phi in T.plaq_of_edges else "straight cycle"
        Vphi = T.verts_of(phi)
        tally = Counter(); examples = defaultdict(list)
        for phi2 in C4:
            if T.dist_sets(Vphi, T.verts_of(phi2), 2*nmax) > 2*nmax: continue
            for lay2 in ((-1, 0) if both_layers else (0,)):
                if lay2 == -1 and (phi2 == phi or (T.verts_of(phi2) & Vphi)): continue     # same layer: distinct components
                fronts = [(phi, -1), (phi2, lay2)]
                for Q in enumerate_Q(T, fronts, nmax):
                    if not junction_ok(T, Q, fronts): continue
                    n = len(Q)
                    sub = ('cap' if (n == 1 and phi2 == phi) else
                           ('exit plaquette' if phi2 in T.plaq_of_edges else 'exit straight'))
                    tally[(n, 'same layer' if lay2 == -1 else 'other layer', sub)] += 1
                    if len(examples[n]) < 3: examples[n].append(([T.pkey[p] for p in Q], 'exit', T.pkey[T.plaq_of_edges[phi2]] if phi2 in T.plaq_of_edges else 'straight', lay2))
        res[kind] = (tally, examples)
        print(f"  entry phi = {kind} at the origin, layer s-1: admissible pairs (J, p') with two four-edge front-pairs, |J|_C <= {nmax}:")
        for key in sorted(tally): print(f"     |J|_C={key[0]}  p' in the {key[1]:11s} {key[2]:15s}: {tally[key]}")
        Z = Counter()
        for key, v in tally.items(): Z[key[0]] += v
        print(f"     Z_n,4 = {dict(sorted(Z.items()))}")
        sys.stdout.flush()
    return res

def partB(L, sizes=(4,), maxfp=6):
    """one-plaquette junctions: c = the plaquette at the origin in plane (0,1); fronts = cycles with the given sizes through V(c)."""
    T = Torus(L)
    o = tuple(0 for _ in L); c = T.pkey.index((o, 0, 1)); Ec = T.plaqs[c]; Vc = T.pverts[c]
    fr = []
    for k in sizes:
        fr += [f for f in T.cycles(k, through=list(Vc))]
    fr = sorted(set(fr), key=lambda f: (len(f), sorted(f)))
    fr = [f for f in fr if T.verts_of(f) & Vc]
    print(f"torus {L}: plaquette c={T.pkey[c]}, {len(fr)} candidate fronts (cycles of sizes {sizes}) through V(c)")
    tally = Counter(); n_ok = 0
    # choose a set of fronts for layer s-1 (pairwise vertex-disjoint) and for layer s (pairwise vertex-disjoint)
    def disjoint_sets(start, chosen, usedV, out):
        out.append(list(chosen))
        if len(chosen) >= 4: return                           # at most one front per vertex of c... (4 vertices) -- enough
        for i in range(start, len(fr)):
            Vf = T.verts_of(fr[i])
            if Vf & usedV: continue
            chosen.append(fr[i]); disjoint_sets(i+1, chosen, usedV | Vf, out); chosen.pop()
    layer_choices = []
    disjoint_sets(0, [], frozenset(), layer_choices)
    print(f"   {len(layer_choices)} choices of vertex-disjoint front sets per layer")
    for A in layer_choices:
        TA = set().union(*A) if A else set()
        for B in layer_choices:
            if len(A) + len(B) > maxfp: continue
            TB = set().union(*B) if B else set()
            cnt = Counter()
            for e in TA: cnt[e] += 1
            for e in TB: cnt[e] += 1
            for e in Ec: cnt[e] += 1
            if any(v < 2 for v in cnt.values()): continue
            nfp = len(A) + len(B)
            prof = (tuple(sorted(len(f) for f in A)), tuple(sorted(len(f) for f in B)))
            tally[(nfp, prof)] += 1; n_ok += 1
    print(f"   admissible one-plaquette junctions (all fronts from the list): {n_ok}")
    for key in sorted(tally):
        print(f"     {key[0]} front-pairs, sizes at s-1 {key[1][0]}, at s {key[1][1]}: {tally[key]}")
    three_all4 = sum(v for (nfp, prof), v in tally.items() if nfp == 3 and all(x == 4 for x in prof[0] + prof[1]))
    print(f"   with exactly three front-pairs, all of four edges: {three_all4}")
    sys.stdout.flush()
    return tally

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else 'A'
    if which == 'A2':      # spatial dimension 2 (d = 3)
        for L in ((5, 5), (5, 6), (6, 6), (4, 5), (4, 4)):
            partA(L, nmax=4)
    if which == 'A3':      # spatial dimension 3 (d = 4)
        for L in ((5, 5, 5), (5, 6, 6), (4, 5, 5)):
            partA(L, nmax=3)
    if which == 'A3n4':
        partA((6, 6, 6), nmax=4)
    if which == 'A4':
        partA((5, 5, 5, 5), nmax=3)
    if which == 'B':
        for L in ((5, 5), (4, 4), (4, 5), (5, 5, 5), (4, 4, 4), (4, 5, 5)):
            partB(L, sizes=(4,), maxfp=8)
    if which == 'B6':
        for L in ((6, 6), (6, 6, 6)):
            partB(L, sizes=(4, 6), maxfp=4)

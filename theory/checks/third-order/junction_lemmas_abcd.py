# Section 13 of Part IV (Lemma mg:lem:betafour): exhaustive checks of the four junction lemmas (a compass; the proofs are in the paper).
# Paper numbering: A = Lemma 13.1(i), D = 13.1(ii), B = 13.1(iii), C = 13.1(iv).
# A junction J at a slice s (without support) is given by its plaquettes P (all at slice s) and its front-pairs (phi, layer), layer in {-1 = s-1, +1 = s}.
# Local conditions checked (Lemma mg:lem:junction and Definition mg:def:junction):
#   n(e) := [e in A] + [e in B] + #{p in P : e in p} >= 2 for every link e of E(J)  (A, B: union of the fronts of layer s-1, s);
#   fronts of the same layer are vertex-disjoint; every front touches (shares a vertex with) a plaquette of P;
#   J is connected: markers (link-connected components of P) and front-pairs, adjacent when they share a vertex (two fronts of the
#   same layer never do).
# Fronts here: all cycles of 4 edges (plaquette boundaries; no straight cycles in Z^n) and of 6 edges.
#  A: two front-pairs of four edges, four plaquettes: exactly the faces of a unit cube other than q, q' (other layer: 5 per cube;
#     same layer: the opposite face), i.e. Z_{4,4} = 12(d-3) = 6(n-2) per front-pair, n = d-1 the spatial dimension.
#  B: one plaquette, exactly three front-pairs, two of them with four edges: impossible with the third of <= 6 edges.
#  C: two plaquettes, exactly three front-pairs, two of them with four edges: impossible with the third of four edges.
#  D: one plaquette, two front-pairs, p with four edges, p' with <= 6 edges and p' != p: exactly the K_S fronts bd(q) xor bd(c), c ~ q, in the
#     other layer (K_S = 4(2n-3) plaquettes sharing a link with q).
import itertools, sys
from collections import defaultdict, Counter

def unit(n, i): return tuple(1 if j == i else 0 for j in range(n))
def add(x, y): return tuple(a + b for a, b in zip(x, y))
def sub(x, y): return tuple(a - b for a, b in zip(x, y))

def plaq_edges(p, n):
    x, i, j = p
    return frozenset([(x, i), (add(x, unit(n, i)), j), (add(x, unit(n, j)), i), (x, j)])

def edge_verts(e, n):
    x, i = e
    return (x, add(x, unit(n, i)))

def verts(E, n):
    return frozenset(v for e in E for v in edge_verts(e, n))

def plaqs_near(n, R):
    out = []
    for x in itertools.product(range(-R, R + 1), repeat=n):
        for i, j in itertools.combinations(range(n), 2):
            out.append((x, i, j))
    return out

def cycles(n, k, center_verts, R):
    """all simple cycles with k edges that pass through a vertex of center_verts and stay in the box [-R-1, R+2]^n."""
    lo, hi = -R - 1, R + 2
    def nbrs(v):
        for i in range(n):
            for s in (1, -1):
                w = add(v, tuple(s if j == i else 0 for j in range(n)))
                if all(lo <= c <= hi for c in w): yield w
    def edge(v, w):
        d = sub(w, v); i = [j for j in range(n) if d[j] != 0][0]
        return (v, i) if d[i] == 1 else (w, i)
    out = set()
    for v0 in center_verts:
        stack = [(v0, [v0])]
        while stack:
            v, path = stack.pop()
            if len(path) == k:
                if v0 in set(nbrs(v)):
                    E = frozenset(edge(path[t], path[(t + 1) % k]) for t in range(k))
                    if len(E) == k: out.add(E)
                continue
            for w in nbrs(v):
                if w in path: continue
                stack.append((w, path + [w]))
    return out

def is_junction(P, fronts, n):
    """P: set of plaquettes; fronts: list of (edge set, layer)."""
    A = set(); B = set()
    for E, lay in fronts:
        (A if lay < 0 else B).update(E)
    # same-layer fronts vertex-disjoint
    for (E1, l1), (E2, l2) in itertools.combinations(fronts, 2):
        if l1 == l2 and verts(E1, n) & verts(E2, n): return False
    cnt = Counter()
    for e in A: cnt[e] += 1
    for e in B: cnt[e] += 1
    for p in P:
        for e in plaq_edges(p, n): cnt[e] += 1
    if any(v < 2 for v in cnt.values()): return False
    VP = verts(set().union(*[plaq_edges(p, n) for p in P]), n) if P else frozenset()
    for E, lay in fronts:
        if not (verts(E, n) & VP): return False
    # markers: link-connected components of P
    P = list(P); parent = {p: p for p in P}
    def f(p):
        while parent[p] != p: parent[p] = parent[parent[p]]; p = parent[p]
        return p
    for p1, p2 in itertools.combinations(P, 2):
        if plaq_edges(p1, n) & plaq_edges(p2, n): parent[f(p1)] = f(p2)
    comps = defaultdict(list)
    for p in P: comps[f(p)].append(p)
    nodes = [(verts(set().union(*[plaq_edges(p, n) for p in C]), n), 'M') for C in comps.values()] + [(verts(E, n), lay) for E, lay in fronts]
    seen = {0}; st = [0]
    while st:
        a = st.pop()
        for b in range(len(nodes)):
            if b in seen: continue
            if nodes[a][1] != 'M' and nodes[b][1] != 'M' and nodes[a][1] == nodes[b][1]: continue
            if nodes[a][0] & nodes[b][0]: seen.add(b); st.append(b)
    return len(seen) == len(nodes)

def is_cube_faces(P, extra, n):
    Y = set(P) | set(extra)
    if len(Y) != 6: return False
    V = verts(set().union(*[plaq_edges(p, n) for p in Y]), n)
    if len(V) != 8: return False
    span = sorted(max(v[k] for v in V) - min(v[k] for v in V) for k in range(n))
    return span == [0]*(n - 3) + [1, 1, 1]

def run(n, R=2):
    o = tuple([0]*n); q = (o, 0, 1); phi = plaq_edges(q, n); Vq = verts(phi, n)
    near = plaqs_near(n, R)
    near1 = [p for p in near if verts(plaq_edges(p, n), n) & Vq]                       # plaquettes touching q
    cyc4 = [plaq_edges(p, n) for p in near]
    # candidate 6-cycles: through vertices within distance 2 of q
    cvs = [v for v in itertools.product(range(-1, 3), repeat=n)]
    cyc6 = list(cycles(n, 6, cvs, R))
    KS = 4*(2*n - 3)
    print(f"spatial dimension n={n} (d={n+1}): {len(near)} plaquettes in the box, {len(cyc6)} six-cycles near q; K_S={KS}")
    # D: one plaquette, two front-pairs p=(phi,-1), p'=(phi',lay') with |phi'| in {4,6}, phi' != phi or lay' != -1
    tallyD = Counter(); badD = []
    for c in near1:
        for phi2 in cyc4 + cyc6:
            for lay in (-1, 1):
                if phi2 == phi and lay == -1: continue
                if is_junction({c}, [(phi, -1), (phi2, lay)], n):
                    if len(phi2) == 4: tallyD['p\' four edges (cap if c=q)'] += 1; continue
                    ok = lay == 1 and c != q and len(plaq_edges(c, n) & phi) == 1 and phi2 == phi ^ plaq_edges(c, n)
                    tallyD['six-edge exit, domino' if ok else 'other'] += 1
                    if not ok: badD.append((c, sorted(phi2), lay))
    print(f"  D: one plaquette, two front-pairs (4 and <=6 edges): {dict(tallyD)}  (expected: 1 cap and K_S={KS} domino exits)")
    # B: one plaquette, three front-pairs: p, and two more with 4 or 6 edges, not both six
    cnt = 0; bad = []
    for c in near1:
        Vc = verts(plaq_edges(c, n), n)
        cand = [(E, lay) for E in cyc4 + cyc6 if verts(E, n) & Vc for lay in (-1, 1) if not (E == phi and lay == -1)]
        for (E1, l1), (E2, l2) in itertools.combinations(cand, 2):
            if len(E1) == 6 and len(E2) == 6: continue
            fr = [(phi, -1), (E1, l1), (E2, l2)]
            if is_junction({c}, fr, n):
                cnt += 1; bad.append((c, len(E1), l1, len(E2), l2))
    print(f"  B: one plaquette, three front-pairs with sizes (4, 4, <=6): {cnt} found (expected 0)")
    # C: two plaquettes, three front-pairs all of four edges
    cnt = 0
    plist = near
    for c1 in near1:
        V1 = verts(plaq_edges(c1, n), n)
        for c2 in plist:
            if c2 <= c1 and c2 in near1: continue
            V2 = verts(plaq_edges(c2, n), n)
            VP = V1 | V2
            cand = [(E, lay) for E in cyc4 if verts(E, n) & VP for lay in (-1, 1) if not (E == phi and lay == -1)]
            # quick filter: every edge of c1, c2 must lie in a front or in the other plaquette
            for (E1, l1), (E2, l2) in itertools.combinations(cand, 2):
                fr = [(phi, -1), (E1, l1), (E2, l2)]
                if is_junction({c1, c2}, fr, n): cnt += 1
    print(f"  C: two plaquettes, three front-pairs all of four edges: {cnt} found (expected 0)")
    # A: two front-pairs of four edges, four plaquettes: search over sets of 4 plaquettes covering the needed links
    tallyA = Counter(); badA = 0
    for phi2 in cyc4:
        for lay in (-1, 1):
            if phi2 == phi and lay == -1: continue
            if lay == -1 and verts(phi2, n) & Vq: continue
            base = Counter()
            for e in phi: base[e] += 1
            for e in phi2: base[e] += 1
            # backtracking over plaquettes at links with count 1, then supersets touching the support (complete for <= 4 plaquettes)
            found = set(); seen = set()
            def rec(X):
                key = frozenset(X)
                if key in seen: return
                seen.add(key)
                c_ = Counter(base)
                for p in X:
                    for e in plaq_edges(p, n): c_[e] += 1
                bad_links = [e for e, v in c_.items() if v == 1]
                if not bad_links and len(X) == 4: found.add(key)
                if len(X) >= 4: return
                if bad_links:
                    e = bad_links[0]; x, i = e
                    cands = []
                    for j in range(n):
                        if j == i: continue
                        for base_pt in (x, sub(x, unit(n, j))):
                            a_, b_ = min(i, j), max(i, j)
                            cands.append((base_pt, a_, b_))
                else:
                    cands = []
                    for e in c_:
                        x, i = e
                        for j in range(n):
                            if j == i: continue
                            for base_pt in (x, sub(x, unit(n, j))):
                                cands.append((base_pt, min(i, j), max(i, j)))
                for p in cands:
                    if p not in X: rec(X | {p})
            rec(frozenset())
            for X in found:
                if not is_junction(set(X), [(phi, -1), (phi2, lay)], n): continue
                qq = [p for p in near if plaq_edges(p, n) == phi2]
                ok = qq and is_cube_faces(X, [q, qq[0]], n) and q not in X and qq[0] not in X
                tallyA[('other layer' if lay == 1 else 'same layer', 'cube' if ok else 'OTHER')] += 1
    print(f"  A: two front-pairs of four edges, four plaquettes: {dict(tallyA)}  (expected 5*2(n-2) other layer, 2(n-2) same layer)")
    sys.stdout.flush()

if __name__ == "__main__":
    for n in [int(a) for a in sys.argv[1:]] or [2, 3]:
        run(n)

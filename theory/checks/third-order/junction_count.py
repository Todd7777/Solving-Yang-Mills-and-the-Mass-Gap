# Section 16 of Part IV (Remark mg:rem:practical): exact Z_{n,4}. Usage: python3 junction_count.py D NMAX (D = d-1)
"""Z_{n,4}: junctions at a slice s with exactly two front-pairs p=(bd q, s-1) (fixed) and p'=(bd q', l') with four edges each, and n plaquettes,
in the spatial lattice Z^D. Local conditions (Lemma mg:lem:junction, Definition mg:def:junction): every link of E(J) lies in at least two
of {front of layer s-1, front of layer s, plaquettes of P} (counted with multiplicity); fronts of one layer vertex-disjoint; each front shares
a vertex with a plaquette of P; J connected (markers = link-components of P, front-pairs; adjacency = common vertex).
Deficiency-driven search on the multiset {q~ (layer s-1), q~' (unassigned, layer s-1 or s), P}: at the smallest link covered once,
branch over which further elements contain it. When no link is covered once, the current multiset is closed; a junction may consist of
several closed parts that touch only at vertices (a shift: the tube ends on q and restarts on q'), so the search then also seeds a further
real plaquette or the exit front at a vertex of the current elements. Duplicates are removed; Z_(n,4) counts distinct junctions.
Results: D=2 (d=3): Z_1=1, Z_2=8, none with 3..7 plaquettes; D=3 (d=4): Z_1=1, Z_2=32, Z_3=0, Z_4=12, Z_5=26, Z_6=302, Z_7=288."""
import sys, itertools
from collections import Counter
D = int(sys.argv[1]); NMAX = int(sys.argv[2])
E = [tuple(1 if k == i else 0 for k in range(D)) for i in range(D)]
def add(x, y): return tuple(a + b for a, b in zip(x, y))
def links(p):
    x, i, j = p
    return ((x, i), (add(x, E[i]), j), (add(x, E[j]), i), (x, j))
def verts(p):
    x, i, j = p
    return {x, add(x, E[i]), add(x, E[j]), add(add(x, E[i]), E[j])}
def plaqs_through(l):
    x, k = l; out = []
    for m in range(D):
        if m == k: continue
        i, j = min(k, m), max(k, m)
        out.append((x, i, j)); out.append((add(x, tuple(-c for c in E[m])), i, j))
    return out
q = (tuple([0] * D), 0, 1)
count = Counter(); kinds = Counter(); done_j = set()
def finish(P, exit_):
    qp, lay = exit_
    if lay == -1 and verts(q) & verts(qp): return          # same layer: fronts vertex-disjoint (and distinct)
    if not P: return
    # each front touches a plaquette of P
    VP = set().union(*(verts(p) for p in P))
    if not (verts(q) & VP) or not (verts(qp) & VP): return
    # connectivity: markers (link components of P) and the two front-pairs, adjacency by common vertex
    P = list(P); comp = list(range(len(P)))
    def find(a):
        while comp[a] != a: comp[a] = comp[comp[a]]; a = comp[a]
        return a
    for a, b in itertools.combinations(range(len(P)), 2):
        if set(links(P[a])) & set(links(P[b])): comp[find(a)] = find(b)
    markers = {}
    for a in range(len(P)): markers.setdefault(find(a), set()).update(verts(P[a]))
    nodes = list(markers.values()) + [verts(q), verts(qp)]
    par = list(range(len(nodes)))
    def f2(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    for a, b in itertools.combinations(range(len(nodes)), 2):
        if a >= len(markers) and b >= len(markers): continue      # two fronts are never adjacent
        if nodes[a] & nodes[b]: par[f2(a)] = f2(b)
    if len({f2(a) for a in range(len(nodes))}) != 1: return
    keyj = (frozenset(P), exit_)
    if keyj in done_j: return
    done_j.add(keyj)
    count[len(P)] += 1
    kinds[(len(P), 'same' if lay == -1 else 'other', 'q in P' if q in P else '', "q' in P" if qp in P else '', 'q=q\'' if qp == q else '')] += 1
seen_states = set()
def plaqs_at_vertex(v):
    out = []
    for i, j in itertools.combinations(range(D), 2):
        for si in (0, 1):
            for sj in (0, 1):
                x = tuple(v[k] - (si if k == i else 0) - (sj if k == j else 0) for k in range(D))
                out.append((x, i, j))
    return out
def rec(cnt, P, exit_, excl):
    key = (P, exit_)
    defi = [l for l, c in cnt.items() if c == 1]
    if not defi:
        if key in seen_states: return
        seen_states.add(key)
        if exit_ is not None: finish(P, exit_)
        # seed a further closed part at a vertex of the current elements (duplicates removed by seen_states / the final set)
        if len(P) >= NMAX: return
        V = set(verts(q)) | set().union(*[verts(p) for p in P]) | (verts(exit_[0]) if exit_ else set())
        seeds = set()
        for v in V:
            for p in plaqs_at_vertex(v):
                if p not in P: seeds.add(('r', p))
                if exit_ is None:
                    for lay in (-1, 1):
                        if not (lay == -1 and p == q): seeds.add(('v', (p, lay)))
        for kind, obj in seeds:
            cnt2 = Counter(cnt)
            if kind == 'r':
                for ll in links(obj): cnt2[ll] += 1
                rec(cnt2, P | {obj}, exit_, frozenset())
            else:
                for ll in links(obj[0]): cnt2[ll] += 1
                rec(cnt2, P, obj, frozenset())
        return
    l = min(defi)
    cands = []
    for p in plaqs_through(l):
        if p not in P and ('r', p) not in excl: cands.append(('r', p))
    if exit_ is None:
        for p in plaqs_through(l):
            for lay in (-1, 1):
                if lay == -1 and p == q: continue
                cands.append(('v', (p, lay)))
    room = NMAX - len(P)
    for r in range(1, len(cands) + 1):
        for T in itertools.combinations(cands, r):
            reals = [c[1] for c in T if c[0] == 'r']; virt = [c[1] for c in T if c[0] == 'v']
            if len(reals) > room or len(virt) > 1: continue
            cnt2 = Counter(cnt)
            for p in reals:
                for ll in links(p): cnt2[ll] += 1
            ex2 = exit_
            if virt:
                ex2 = virt[0]
                for ll in links(ex2[0]): cnt2[ll] += 1
            # every element through l is now decided: exclude the other real candidates (virtual ones lapse once the exit is set)
            excl2 = excl | {c for c in cands if c[0] == 'r' and c not in T}
            rec(cnt2, P | set(reals), ex2, excl2)
rec(Counter(links(q)), frozenset(), None, frozenset())
print(f"D={D}, n<= {NMAX}: Z_(n,4) by n: {dict(sorted(count.items()))}")
for k, v in sorted(kinds.items()): print('  ', k, v)

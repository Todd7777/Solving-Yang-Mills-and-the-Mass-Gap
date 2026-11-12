# Written by the independent referee of Section 15 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs.
# Referee of Section 15 (2026-11-12). A compass, not part of the proofs.
# Combinatorial test of Definition 15.6 (moves and partners) and of step (1) of the proof of Theorem 15.9, from the definitions of
# Sections 7-8 (fronts, clusters, markers, front-pairs, junctions as components of H_s), on Z x Z^3 (d = 4), for explicit sets X'.
# (E1) a turn whose rest R has a lid touching y at sigma-1: the partner is not a hop in the sense of Definition 15.6;
# (E2) a hop whose rest has a lid touching v at sigma: the partner is not a turn in the sense of Definition 15.6;
# (E3) a hop (q,c,4) and a turn (q',c',5) of the rest: flipping the hop removes the second move (its entry front becomes touched),
#      so the move vertices, hence J^natural, are not invariant under the flip;
# plus, for each example: no free links in X' and in X'^* (Lemma 15.7(i)), and the size of E(X') xor E(X'^*) (the count behind c_Gamma).
import itertools
from geom import add, unit, boundary, cube_faces

def sp(s, p):            # spatial plaquette p=(x,i,j) at slice s
    return ('S', s, p)

def tp(s, link):         # temporal plaquette of layer s over the spatial link (x,i)
    return ('T', s, link)

def plinks(P):
    if P[0] == 'S':
        _, s, (x, i, j) = P
        return [('L', s, l) for l in boundary((x, i, j))]
    _, s, (x, i) = P
    y = add(x, unit(i))
    return [('L', s, (x, i)), ('L', s + 1, (x, i)), ('V', s, x), ('V', s, y)]

def verts(link):
    x, i = link
    return {x, add(x, unit(i))}

def tube(q, layers):
    return {tp(s, l) for s in layers for l in boundary(q)}

def components(items, adj):
    items = list(items); seen = set(); comps = []
    for a in items:
        if a in seen:
            continue
        comp, stack = [], [a]; seen.add(a)
        while stack:
            b = stack.pop(); comp.append(b)
            for c in items:
                if c not in seen and adj(b, c):
                    seen.add(c); stack.append(c)
        comps.append(comp)
    return comps

def free_links(X, Ssupp):
    cnt = {}
    for P in X:
        for e in plinks(P):
            cnt[e] = cnt.get(e, 0) + 1
    return [e for e, k in cnt.items() if k == 1 and e not in Ssupp]

def junctions(X, s, supports):
    """Junctions at slice s: components of H_s (markers = clusters at s and supports at s; front-pairs)."""
    T = {lam: {P[2] for P in X if P[0] == 'T' and P[1] == lam} for lam in (s - 2, s - 1, s)}
    fronts = {lam: [frozenset(cmp) for cmp in components(T[lam], lambda a, b: bool(verts(a) & verts(b)))] for lam in T}
    spat = [P[2] for P in X if P[0] == 'S' and P[1] == s]
    clusters = [frozenset(cmp) for cmp in components(spat, lambda a, b: bool(set(boundary(a)) & set(boundary(b))))]
    markers = [('cl', C, set(l for p in C for l in boundary(p))) for C in clusters]
    for (ss, E) in supports:
        if ss == s:
            markers.append(('supp', frozenset(E), set(E)))
    V = lambda E: {v for l in E for v in verts(l)}
    fps = [('fp', (phi, lam), set(phi)) for lam in (s - 1, s) for phi in fronts[lam]
           if any(V(phi) & V(m[2]) for m in markers)]
    nodes = markers + fps
    comps = components(range(len(nodes)), lambda a, b: bool(V(nodes[a][2]) & V(nodes[b][2])))
    out = []
    for cmp in comps:
        cl = frozenset(nodes[k][1] for k in cmp if nodes[k][0] == 'cl')
        fp = frozenset(nodes[k][1] for k in cmp if nodes[k][0] == 'fp')
        sup = [nodes[k][1] for k in cmp if nodes[k][0] == 'supp']
        out.append((cl, fp, sup))
    return out, fronts, markers

def has_turn(X, q, c, sig, supports):
    Gam = [F for F in (cube_faces(add(q[0], d)) for d in itertools.product((-1, 0), repeat=3)) if q in F and c in F][0]
    rest = frozenset(Gam) - {q, c}
    J, fr, mk = junctions(X, sig, supports)
    target = (frozenset({rest}), frozenset({(frozenset(boundary(q)), sig - 1), (frozenset(boundary(c)), sig)}))
    ok_sig = any((cl, fp) == target and not sup for cl, fp, sup in J)
    J1, fr1, mk1 = junctions(X, sig - 1, supports)
    dq = frozenset(boundary(q))
    V = lambda E: {v for l in E for v in verts(l)}
    untouched = dq in fr1[sig - 2] and not any(V(dq) & V(m[2]) for m in mk1)
    return ok_sig and untouched

def has_hop(X, q, c, sig, supports):
    dq, dc = frozenset(boundary(q)), frozenset(boundary(c))
    hexf = frozenset(set(dq) ^ set(dc))
    J1 = junctions(X, sig - 1, supports)[0]; J2 = junctions(X, sig, supports)[0]
    t1 = (frozenset({frozenset({c})}), frozenset({(dq, sig - 2), (hexf, sig - 1)}))
    t2 = (frozenset({frozenset({q})}), frozenset({(hexf, sig - 1), (dc, sig)}))
    return any((cl, fp) == t1 for cl, fp, _ in J1) and any((cl, fp) == t2 for cl, fp, _ in J2)

def local(kind, q, c, sig):
    Gam = [F for F in (cube_faces(add(q[0], d)) for d in itertools.product((-1, 0), repeat=3)) if q in F and c in F][0]
    if kind == 'turn':
        return {sp(sig, p) for p in Gam if p not in (q, c)} | tube(q, [sig - 1])
    hexf = set(boundary(q)) ^ set(boundary(c))
    return {sp(sig - 1, c), sp(sig, q)} | {tp(sig - 1, l) for l in hexf}

def report(name, X, q, c, sig, supports, kind):
    Ssupp = {('L', s, l) for s, E in supports for l in E}
    loc = local(kind, q, c, sig)
    assert loc <= X
    other = 'hop' if kind == 'turn' else 'turn'
    Xs = (X - loc) | local(other, q, c, sig)
    E = lambda Y: {e for P in Y for e in plinks(P)}
    print(f'{name}: free links X\' {len(free_links(X, Ssupp))}, X\'* {len(free_links(Xs, Ssupp))};'
          f' |E(X\') xor E(X\'*)| = {len(E(X) ^ E(Xs))}')
    print(f'   X\'  has turn: {has_turn(X, q, c, sig, supports)}, hop: {has_hop(X, q, c, sig, supports)}')
    print(f'   X\'* has turn: {has_turn(Xs, q, c, sig, supports)}, hop: {has_hop(Xs, q, c, sig, supports)}')
    return Xs

x0 = (0, 0, 0)
q = (x0, 0, 1); c = (x0, 0, 2); sig = 4                 # sigma(q,c) = -1; Gamma = unit cube at the origin; y = (0,0,1), v = (0,1,1)
y, v = (0, 0, 1), (0, 1, 1)
main_turn = tube(q, range(0, sig - 1)) | local('turn', q, c, sig) | tube(c, range(sig, 8))
main_hop = tube(q, range(0, sig - 1)) | local('hop', q, c, sig) | tube(c, range(sig, 8))
dq, dc = list(boundary(q)), list(boundary(c))

# E1: turn + a lid p at sigma-1 = 3 through y (p spanned by e2, e3 at y - e2; no link of Gamma)
p = ((0, -1, 1), 1, 2)
X1 = main_turn | tube(p, range(0, 3)) | {sp(3, p)}
sup1 = [(0, dq + list(boundary(p))), (8, dc)]
report('E1 (turn, lid at y at sigma-1)', X1, q, c, sig, sup1, 'turn')
J = junctions((X1 - local('turn', q, c, sig)) | local('hop', q, c, sig), 3, sup1)[0]
print('   junction of the partner at sigma-1 containing c:', [ (len(cl), len(fp)) for cl, fp, _ in J if any(c in C for C in cl)], '(clusters, front-pairs)')

# E2: hop + a lid p2 at sigma = 4 through v (p2 spanned by e2, e3 at v)
p2 = (v, 1, 2)
X2 = main_hop | tube(p2, range(0, 4)) | {sp(4, p2)}
sup2 = [(0, dq + list(boundary(p2))), (8, dc)]
report('E2 (hop, lid at v at sigma)', X2, q, c, sig, sup2, 'hop')

# E3: hop (q,c,4) and, in the rest, a tube on q' = p2 through layers 0..4 and a turn (q',c',5) through the cube Gamma' at v - e1
qq = p2; cc = ((-1, 1, 1), 0, 1)
Gp = cube_faces((-1, 1, 1)); assert qq in Gp and cc in Gp
X3 = main_hop | tube(qq, range(0, 5)) | {sp(5, f) for f in Gp if f not in (qq, cc)} | tube(cc, range(5, 8))
sup3 = [(0, dq + list(boundary(qq))), (8, dc + list(boundary(cc)))]
X3s = report('E3 (hop at (3,4); rest has a turn at (4,5))', X3, q, c, sig, sup3, 'hop')
print(f'   second move (q\',c\',5): turn in X\': {has_turn(X3, qq, cc, 5, sup3)};  turn in X\'*: {has_turn(X3s, qq, cc, 5, sup3)}')

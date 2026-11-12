#!/usr/bin/env python3
"""The families of Lemma mg:lem:threelow, built literally from their descriptions in the statement.

Junction keys are (frozenset(P), frozenset({(front, layer)})), as in enum_junctions.py.  Layers: lam and lb = 1 - lam.
"""
from lattice16b import Lattice


def PE(lat, p):
    return frozenset(lat.plaq_edges(p))


def PV(lat, p):
    return frozenset(lat.plaq_verts(p))


def near(lat, q):
    """plaquettes != q sharing a vertex with q"""
    s = set()
    for v in lat.plaq_verts(q):
        s.update(lat.plaqs_at_vertex(v))
    s.discard(q)
    return s


def link_nbrs(lat, r):
    s = set()
    for e in lat.plaq_edges(r):
        s.update(lat.plaqs_at_edge(e))
    s.discard(r)
    return s


def opposite_edge(lat, c, e):
    E = lat.plaq_edges(c)
    i = E.index(e)
    return E[(i + 2) % 4]


def cubes_containing(lat, q):
    """unit cubes containing q, as (list of six faces, dict face -> opposite face, k = direction perpendicular to q)"""
    v = lat.plaq_corner(q)
    i, j = lat.plaq_dirs(q)
    out = []
    for k in range(lat.D):
        if k in (i, j):
            continue
        for w in (v, lat.shm[v][k]):
            a, b, c = sorted((i, j, k))
            sh = lat.sh
            pairs = [(lat.plaq(w, a, b), lat.plaq(sh[w][c], a, b)),
                     (lat.plaq(w, a, c), lat.plaq(sh[w][b], a, c)),
                     (lat.plaq(w, b, c), lat.plaq(sh[w][a], b, c))]
            faces = [x for pr in pairs for x in pr]
            opp = {}
            for x, y in pairs:
                opp[x] = y
                opp[y] = x
            assert q in faces
            out.append((faces, opp, k))
    return out


def key(P, Phi):
    return (frozenset(P), frozenset(Phi))


def family_i(lat, q, lam):
    lb = 1 - lam
    fam = {}
    Nq = near(lat, q)
    dq = PE(lat, q)
    # Caps, p alone in its layer: unordered pairs {q2,q3} of vertex-disjoint plaquettes != q sharing a vertex with q
    Nl = sorted(Nq)
    for a in range(len(Nl)):
        for b in range(a + 1, len(Nl)):
            q2, q3 = Nl[a], Nl[b]
            if PV(lat, q2) & PV(lat, q3):
                continue
            fam[key({q, q2, q3}, {(dq, lam), (PE(lat, q2), lb), (PE(lat, q3), lb)})] = "caps: p alone"
    # Caps, p and (dq2, lam) in one layer: q2 vertex-disjoint from q, q3 sharing a vertex with q and q2
    for q3 in Nq:
        for q2 in near(lat, q3):
            if q2 == q or (PV(lat, q2) & PV(lat, q)):
                continue
            fam[key({q, q2, q3}, {(dq, lam), (PE(lat, q2), lam), (PE(lat, q3), lb)})] = "caps: p with q2"
    # Cubes
    for faces, opp, k in cubes_containing(lat, q):
        qo = opp[q]
        sides = [x for x in faces if x not in (q, qo)]
        for x in sides:   # q and its opposite face in lam, one side face in lb
            P = set(faces) - {q, qo, x}
            fam[key(P, {(dq, lam), (PE(lat, qo), lam), (PE(lat, x), lb)})] = "cube: q, opp(q) in lam"
        done = set()
        for x in sides:   # two opposite faces adjacent to q in lb
            y = opp[x]
            if frozenset((x, y)) in done:
                continue
            done.add(frozenset((x, y)))
            P = set(faces) - {q, x, y}
            fam[key(P, {(dq, lam), (PE(lat, x), lb), (PE(lat, y), lb)})] = "cube: opposite pair in lb"
    return fam


def _ab_configs(lat, q):
    """Configurations (r, c, e, f, r', kind) of (a)/(b) with q in {r, r'}."""
    out = []
    cand_c = set()
    for v in lat.plaq_verts(q):
        cand_c.update(lat.plaqs_at_vertex(v))
    # all (r, c) with c sharing link e with r, and {r, r'} containing q
    seen = set()
    rc_pairs = set()
    for c in cand_c | link_nbrs(lat, q):
        for e in lat.plaq_edges(c):
            for r in lat.plaqs_at_edge(e):
                if r != c:
                    rc_pairs.add((r, c, e))
    for c in link_nbrs(lat, q):
        e = next(iter(PE(lat, q) & PE(lat, c)))
        rc_pairs.add((q, c, e))
    for (r, c, e) in rc_pairs:
        f = opposite_edge(lat, c, e)
        fv = set(lat.ends(f))
        # (a): r' contains f, r' != c
        for rp in lat.plaqs_at_edge(f):
            if rp == c:
                continue
            if q in (r, rp):
                t = (r, c, rp, 'a')
                if t not in seen:
                    seen.add(t)
                    out.append((r, c, e, f, rp, 'a'))
        # (b): r' shares exactly one vertex with c, an endpoint of f
        for v in fv:
            for rp in lat.plaqs_at_vertex(v):
                if rp == c:
                    continue
                if len(PV(lat, rp) & PV(lat, c)) != 1:
                    continue
                if q in (r, rp):
                    t = (r, c, rp, 'b')
                    if t not in seen:
                        seen.add(t)
                        out.append((r, c, e, f, rp, 'b'))
    return out


def family_ii(lat, q, lam, seven=False):
    lb = 1 - lam
    fam = {}
    for (r, c, e, f, rp, kind) in _ab_configs(lat, q):
        qp = rp if r == q else r
        big = (PE(lat, r) | PE(lat, c)) if seven else (PE(lat, r) ^ PE(lat, c))
        k = key({c, rp}, {(PE(lat, q), lam), (PE(lat, qp), lam), (big, lb)})
        lab = f"({kind}) r={'q' if r == q else 'q_prime'}"
        assert k not in fam or fam[k] == lab, "a configuration has two descriptions"
        fam[k] = lab
    if not seven:
        for faces, opp, kdir in cubes_containing(lat, q):
            qo = opp[q]
            sides = [x for x in faces if x not in (q, qo)]
            # edges of the cube perpendicular to q: links in direction kdir lying in two side faces
            for a in range(4):
                for b in range(a + 1, 4):
                    c1, c2 = sides[a], sides[b]
                    common = PE(lat, c1) & PE(lat, c2)
                    if len(common) != 1:
                        continue
                    e = next(iter(common))
                    if lat.direction(e) != kdir:
                        continue
                    P = set(sides) - {c1, c2}
                    k = key(P, {(PE(lat, q), lam), (PE(lat, qo), lam), (PE(lat, c1) ^ PE(lat, c2), lb)})
                    fam[k] = "(c) cube"
    return fam

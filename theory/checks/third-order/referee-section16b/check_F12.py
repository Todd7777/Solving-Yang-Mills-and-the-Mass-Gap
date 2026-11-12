#!/usr/bin/env python3
"""Facts (F1), (F2) of Section 16 and the structure of six- and seven-edge fronts, by enumeration.

For every front with six edges through the origin of the torus (Z/L)^D (found from the definition by fronts16b):
  * it is a cycle, and exactly one of: planar domino, bent domino, skew hexagon (defined from coordinates);
  * (F1) every plaquette shares <= 3 links with it, and exactly 3 iff Sigma = d r (sym.diff.) d c with c sharing a link with r;
  * (F2) planar domino: only c1, c2 share >= 2 links; bent domino: exactly two further plaquettes share two links,
         they are the faces of the cube perpendicular to e, opposite, each containing one endpoint of e;
         skew hexagon: each face of its cube shares exactly two consecutive links, every other plaquette <= 1.
For every front with seven edges through the origin:
  * it is d r0 u d c0 for exactly one pair of plaquettes sharing a link e0;
  * a plaquette other than r0, c0 shares <= 2 links with it, and 2 only for the two faces perpendicular to e0 of the
    cube of a bent pair (none for a coplanar pair).
Usage: python3 check_F12.py D L
"""
import sys
from collections import Counter
from lattice16b import Lattice
from fronts16b import fronts_through_vertex


def cubes_containing(lat, p):
    v = lat.plaq_corner(p)
    i, j = lat.plaq_dirs(p)
    for k in range(lat.D):
        if k in (i, j):
            continue
        dirs = tuple(sorted((i, j, k)))
        yield (v, dirs, k)
        yield (lat.shm[v][k], dirs, k)


def cube_faces(lat, v, dirs):
    i, j, k = dirs
    sh = lat.sh
    return [lat.plaq(v, i, j), lat.plaq(sh[v][k], i, j), lat.plaq(v, i, k), lat.plaq(sh[v][j], i, k),
            lat.plaq(v, j, k), lat.plaq(sh[v][i], j, k)]


def cube_verts(lat, v, dirs):
    i, j, k = dirs
    sh = lat.sh
    out = []
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                w = v
                if a:
                    w = sh[w][i]
                if b:
                    w = sh[w][j]
                if c:
                    w = sh[w][k]
                out.append(w)
    return out


def all_cubes_near(lat, verts):
    """All unit cubes (corner, dirs) having a vertex in verts."""
    out = set()
    import itertools
    for w in verts:
        for dirs in itertools.combinations(range(lat.D), 3):
            # corners v with w in cube(v, dirs): v = w - (subset of dirs)
            for a in (0, 1):
                for b in (0, 1):
                    for c in (0, 1):
                        v = w
                        if a:
                            v = lat.shm[v][dirs[0]]
                        if b:
                            v = lat.shm[v][dirs[1]]
                        if c:
                            v = lat.shm[v][dirs[2]]
                        out.add((v, dirs))
    return out


def plaquettes_near(lat, verts):
    s = set()
    for w in verts:
        s.update(lat.plaqs_at_vertex(w))
    return s


def link_neighbours(lat, r):
    """plaquettes sharing a link with r (other than r)."""
    s = set()
    for e in lat.plaq_edges(r):
        for c in lat.plaqs_at_edge(e):
            if c != r:
                s.add(c)
    return s


def main(D, L):
    lat = Lattice(D, L)
    o = lat.origin()
    pe = {}

    def PE(p):
        r = pe.get(p)
        if r is None:
            r = frozenset(lat.plaq_edges(p))
            pe[p] = r
        return r

    # ---------------- six-edge fronts ----------------
    six = fronts_through_vertex(lat, o, 6)
    kinds = Counter()
    bad = 0
    for S in six:
        V = lat.verts_of_edges(S)
        deg = Counter()
        for e in S:
            a, b = lat.ends(e)
            deg[a] += 1
            deg[b] += 1
        assert len(V) == 6 and set(deg.values()) == {2}, "six-edge front not a cycle"
        near = plaquettes_near(lat, V)
        share = {r: len(PE(r) & S) for r in near}
        # dominoes: pairs (r, c) sharing a link with dr ^ dc = S
        doms = set()
        for r in near:
            if share[r] == 3:
                for c in link_neighbours(lat, r):
                    if (PE(r) ^ PE(c)) == S:
                        doms.add(frozenset((r, c)))
        # F1
        f1_ok = max(share.values()) <= 3
        for r in near:
            has = any(r in dpair for dpair in doms)
            if (share[r] == 3) != has:
                f1_ok = False
        # skew hexagon: V = 6 vertices of a unit cube missing two antipodal ones
        skew_cubes = []
        for (v, dirs) in all_cubes_near(lat, V):
            cv = cube_verts(lat, v, dirs)
            if V <= set(cv):
                miss = [w for w in cv if w not in V]
                # antipodal in the cube: indices 0 and 7 in cube_verts ordering, etc.
                idx = sorted(cv.index(w) for w in miss)
                if idx[0] + idx[1] == 7:
                    skew_cubes.append((v, dirs))
        if len(doms) == 1 and not skew_cubes:
            (r, c) = tuple(next(iter(doms)))
            planar = lat.plaq_dirs(r) == lat.plaq_dirs(c)
            kind = "planar domino" if planar else "bent domino"
        elif not doms and len(skew_cubes) == 1:
            kind = "skew hexagon"
        else:
            kind = f"OTHER doms={len(doms)} skew={len(skew_cubes)}"
            bad += 1
        kinds[kind] += 1
        # F2
        f2_ok = True
        two = sorted(r for r in near if share[r] == 2)
        if kind == "planar domino":
            pair = next(iter(doms))
            f2_ok = all(share[r] <= 1 for r in near if r not in pair)
        elif kind == "bent domino":
            r, c = tuple(next(iter(doms)))
            e = next(iter(PE(r) & PE(c)))
            cubes = [(v, dirs) for (v, dirs, k) in cubes_containing(lat, r) if c in cube_faces(lat, v, dirs)]
            assert len(cubes) == 1
            v, dirs = cubes[0]
            faces = cube_faces(lat, v, dirs)
            i = lat.direction(e)
            perp = [f for f in faces if i not in lat.plaq_dirs(f)]
            ends = set(lat.ends(e))
            cond = (sorted(perp) == two and len(perp) == 2
                    and not (set(lat.plaq_verts(perp[0])) & set(lat.plaq_verts(perp[1])))
                    and all(len(set(lat.plaq_verts(f)) & ends) == 1 for f in perp)
                    and all(share[x] <= 1 for x in near if x not in (r, c) and x not in perp))
            f2_ok = cond
        elif kind == "skew hexagon":
            v, dirs = skew_cubes[0]
            faces = cube_faces(lat, v, dirs)
            cond = sorted(faces) == two
            for f in faces:
                l2 = list(PE(f) & S)
                if len(l2) != 2 or not (set(lat.ends(l2[0])) & set(lat.ends(l2[1]))):
                    cond = False
            cond = cond and all(share[x] <= 1 for x in near if x not in faces)
            f2_ok = cond
        if not (f1_ok and f2_ok):
            bad += 1
            print("FAIL", kind, f1_ok, f2_ok)
    print(f"D={D} L={L}: {len(six)} six-edge fronts through a vertex: {dict(kinds)}; F1/F2/classification failures: {bad}")

    # ---------------- seven-edge fronts ----------------
    seven = fronts_through_vertex(lat, o, 7)
    kinds7 = Counter()
    bad7 = 0
    for T in seven:
        V = lat.verts_of_edges(T)
        near = plaquettes_near(lat, V)
        pairs = set()
        for r in near:
            if PE(r) <= T:
                for c in link_neighbours(lat, r):
                    if (PE(r) | PE(c)) == T:
                        pairs.add(frozenset((r, c)))
        if len(pairs) != 1:
            bad7 += 1
            kinds7["OTHER"] += 1
            continue
        r0, c0 = tuple(next(iter(pairs)))
        e0 = next(iter(PE(r0) & PE(c0)))
        planar = lat.plaq_dirs(r0) == lat.plaq_dirs(c0)
        kinds7["coplanar pair" if planar else "bent pair"] += 1
        share = {x: len(PE(x) & T) for x in near if x not in (r0, c0)}
        two = sorted(x for x, m in share.items() if m == 2)
        ok = max(share.values()) <= 2
        if planar:
            ok = ok and not two
        else:
            cubes = [(v, dirs) for (v, dirs, k) in cubes_containing(lat, r0) if c0 in cube_faces(lat, v, dirs)]
            v, dirs = cubes[0]
            i = lat.direction(e0)
            perp = sorted(f for f in cube_faces(lat, v, dirs) if i not in lat.plaq_dirs(f))
            ok = ok and two == perp
        if not ok:
            bad7 += 1
    print(f"D={D} L={L}: {len(seven)} seven-edge fronts through a vertex: {dict(kinds7)}; failures: {bad7}")


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))

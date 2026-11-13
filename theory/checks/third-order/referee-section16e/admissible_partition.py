# Referee check (Part IV, Lemma mg:lem:uppertails, Remark mg:rem:upperexact): the admissible sets, counted from the definition.
#
# Usage: python3 admissible_partition.py d A|Q nmax [--list]
#
# Definition (Lemma mg:lem:uppertails).  (S, C) = (S_A, {p0}) with p0 the plaquette at the origin in directions (1,2), S_A = E(p0);
# or (S_Q, Q) with Q the unit 3-cube at the origin in directions (0,1,2), C its six faces, S_Q its twelve edges.  A finite set Y of
# plaquettes of Z^d is admissible if (1) Y and C are disjoint, (2) every link of E(Y) \ S lies in at least two plaquettes of Y, and
# (3) every link-connected component of Y contains a plaquette of P_S \ C (P_S = plaquettes with a link in S).
#
# Method (written for this report; it does not use the authors' code).  A component of an admissible set is admissible and connected,
# and an admissible set is the union of link-disjoint connected admissible sets.  The connected ones are found by a partition search:
# the root is the smallest plaquette of P_S \ C in the set (smaller roots are excluded); in a state Y with a deficient link (a link
# of E(Y) \ S in exactly one plaquette of Y) the final set contains a further plaquette through the smallest deficient link, and we
# branch on which candidate is the first one in the final set (include it, exclude the earlier ones); in a state without deficient
# link Y itself is recorded, and a larger connected final set contains a plaquette sharing a link with Y, on which we branch in the
# same way.  Each connected admissible set is reached exactly once (asserted).  Pruning: a superset must cover every deficient link
# by a new plaquette, and a plaquette has four links, so it needs at least ceil(#deficient/4) more plaquettes.
import sys
from collections import Counter, defaultdict

TORUS = None   # None: Z^d; an integer L >= 3: the periodic lattice (Z/LZ)^d (option --torus L)

def geometry(d, kind):
    E = [tuple(1 if i == j else 0 for i in range(d)) for j in range(d)]
    if TORUS is None:
        add = lambda x, y: tuple(a + b for a, b in zip(x, y))
        sub = lambda x, y: tuple(a - b for a, b in zip(x, y))
    else:
        add = lambda x, y: tuple((a + b) % TORUS for a, b in zip(x, y))
        sub = lambda x, y: tuple((a - b) % TORUS for a, b in zip(x, y))
    def links(p):
        x, m, n = p
        return ((x, m), (add(x, E[m]), n), (add(x, E[n]), m), (x, n))
    def plaqs(l):
        x, m = l
        out = []
        for n in range(d):
            if n == m: continue
            a, b = min(m, n), max(m, n)
            out.append((x, a, b)); out.append((sub(x, E[n]), a, b))
        return out
    O = (0,) * d
    if kind == 'A':
        C = {(O, 1, 2)}
    else:
        C = set()
        for a, b in ((0, 1), (0, 2), (1, 2)):
            c = 3 - a - b
            C.add((O, a, b)); C.add((E[c], a, b))
    S = {l for p in C for l in links(p)}
    PS = {p for l in S for p in plaqs(l)}
    return links, plaqs, C, S, PS, E, add

def connected_admissible(d, kind, nmax):
    links, plaqs, C, S, PS, E, add = geometry(d, kind)
    roots = sorted(PS - C)
    found = {}
    cnt = Counter(); Y = []; Yset = set(); excl = set(C); deficient = set()
    def add_p(p):
        Y.append(p); Yset.add(p)
        for l in links(p):
            c = cnt[l] + 1; cnt[l] = c
            if l in S: continue
            if c == 1: deficient.add(l)
            elif c == 2: deficient.discard(l)
    def rem_p(p):
        q = Y.pop(); assert q == p; Yset.discard(p)
        for l in links(p):
            c = cnt[l] - 1
            if c: cnt[l] = c
            else: del cnt[l]
            if l in S: continue
            if c == 1: deficient.add(l)
            elif c == 0: deficient.discard(l)
    nodes = [0]
    def branch(cands):
        added = []
        for c in cands:
            add_p(c); rec(); rem_p(c)
            excl.add(c); added.append(c)
        for c in added: excl.discard(c)
    def rec():
        nodes[0] += 1
        n = len(Y)
        if n + (len(deficient) + 3) // 4 > nmax: return
        if deficient:
            l = min(deficient)
            branch([p for p in plaqs(l) if p not in Yset and p not in excl])
        else:
            key = frozenset(Y)
            assert key not in found, "a set was reached twice"
            found[key] = True
            if n < nmax:
                branch(sorted({q for p in Y for l in links(p) for q in plaqs(l)} - Yset - excl))
    for r in roots:
        add_p(r); rec(); rem_p(r)
        excl.add(r)
    return list(found), nodes[0], (links, plaqs, C, S, PS, E, add)

def all_admissible(d, kind, nmax):
    conn, nodes, geo = connected_admissible(d, kind, nmax)
    links = geo[0]
    E_of = {Y: frozenset(l for p in Y for l in links(p)) for Y in conn}
    allsets = [(Y,) for Y in conn]
    # two components (each has >= 5 plaquettes, so at most two when nmax <= 14)
    assert nmax <= 14
    byn = defaultdict(list)
    for Y in conn: byn[len(Y)].append(Y)
    for i, Y1 in enumerate(conn):
        for Y2 in conn[i + 1:]:
            if len(Y1) + len(Y2) <= nmax and not (E_of[Y1] & E_of[Y2]):
                allsets.append((Y1, Y2))
    return conn, allsets, nodes, geo

def is_cap(Y, geo):
    """Y = the five faces other than F of a unit cube, with F in C"""
    links, plaqs, C, S, PS, E, add = geo
    if len(Y) != 5: return False
    for F in C:
        x, a, b = F
        d = len(x)
        for c in range(d):
            if c in (a, b): continue
            for base in (x, tuple((xi - (1 if i == c else 0)) % (TORUS or 10**9) if TORUS else xi - (1 if i == c else 0) for i, xi in enumerate(x))):
                dirs = sorted((a, b, c))
                faces = set()
                for (m, n) in ((dirs[0], dirs[1]), (dirs[0], dirs[2]), (dirs[1], dirs[2])):
                    k = 3 - dirs.index(m) - dirs.index(n)
                    faces.add((base, m, n)); faces.add((add(base, E[dirs[k]]), m, n))
                if F in faces and faces - {F} == set(Y):
                    return True
    return False

def is_cube(Y, geo):
    links, plaqs, C, S, PS, E, add = geo
    if len(Y) != 6: return False
    pts = [v for p in Y for v in (p[0],)]
    # a unit cube: every link of E(Y) in exactly two plaquettes of Y, and 12 links, 8 vertices
    cnt = Counter(l for p in Y for l in links(p))
    verts = {l[0] for l in cnt} | {add(l[0], E[l[1]]) for l in cnt}
    return len(cnt) == 12 and all(c == 2 for c in cnt.values()) and len(verts) == 8

if __name__ == "__main__":
    import time
    d, kind, nmax = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
    if '--torus' in sys.argv:
        TORUS = int(sys.argv[sys.argv.index('--torus') + 1]); assert TORUS >= 3
        print(f"periodic lattice with all sides {TORUS}")
    t0 = time.time()
    conn, allsets, nodes, geo = all_admissible(d, kind, nmax)
    links, plaqs, C, S, PS, E, add = geo
    print(f"d={d} (S,C)=({'S_A,{p0}' if kind == 'A' else 'S_Q,Q'}): |S|={len(S)} |C|={len(C)} |P_S|={len(PS)} J'={len(PS - C)}; nmax={nmax}; "
          f"search nodes {nodes}, {time.time() - t0:.1f} s")
    cc = Counter(len(Y) for Y in conn)
    ac = Counter(sum(len(Y) for Y in s) for s in allsets)
    two = Counter(sum(len(Y) for Y in s) for s in allsets if len(s) == 2)
    caps = [Y for Y in conn if is_cap(Y, geo)]
    cubes = [Y for Y in conn if is_cube(Y, geo)]
    print("  connected admissible sets by size:", dict(sorted(cc.items())))
    print("  all admissible sets by size      :", dict(sorted(ac.items())), " of which with two components:", dict(sorted(two.items())))
    print(f"  caps: {len(caps)}; five-plaquette sets that are not caps: {cc[5] - len(caps)}; six-plaquette sets that are unit cubes: {len(cubes)} of {cc[6]}")
    nn = {n: ac[n] - (len(caps) if n == 5 else 0) for n in range(5, nmax + 1)}
    print("  N^cap =", len(caps), " N_n (non-caps, n = 5..nmax):", [nn[n] for n in range(5, nmax + 1)])
    if '--list' in sys.argv:
        for s in sorted(allsets, key=lambda s: (sum(len(Y) for Y in s), sorted(p for Y in s for p in Y))):
            print("   ", sum(len(Y) for Y in s), [sorted(Y) for Y in s])

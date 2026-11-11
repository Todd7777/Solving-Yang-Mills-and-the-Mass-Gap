# Section 11 of Part IV: exhaustive checks of the two combinatorial lemmas, on spatial tori with all sides >= 6
# (so that the fronts with four edges are the plaquette boundaries). A compass, not part of the proofs.
#  (A) Lemma mg:lem:threeplaq: no junction without support has exactly two front-pairs, both with four edges, and exactly three
#      plaquettes (Z_{3,4} = 0). Enumeration: entry phi = boundary of q1 at layer s-1; exit phi' = boundary of a plaquette q2 within
#      distance R of q1, at layer s or s-1; plaquette sets Q of size 3 among the plaquettes with at least two edges in phi u phi'
#      (Lemma mg:lem:twofront(i)); admissibility n(e) >= 2 on E(J), vertex-disjoint fronts of one layer, every front touched by a
#      plaquette of Q, connectivity of the junction. Also the counts for |Q| = 1, 2 (cap; shifts and jumps), as a calibration.
#  (B) Lemma mg:lem:oneplaq: a junction without support with exactly one plaquette c and exactly three front-pairs has a front-pair
#      with at least five edges. Enumeration: U0 = union of vertex-disjoint fronts of layer s-1 through vertices of c (sizes in
#      SIZES0), D a subset of the boundary of c, U1 = U0 xor D; the junction is admissible iff boundary(c) is contained in U0 u U1
#      and every component of U1 is a front (minimum degree >= 2) touching c. All admissible junctions with >= 3 front-pairs are
#      classified by the sizes of their fronts.
import itertools, sys, os
from collections import Counter, defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lower-bound'))
from small_fronts import torus_graph
from junctions_small import fronts_of_size

def setup(shape):
    n = len(shape)
    sites, edges, ends, at, nbr = torus_graph(shape)
    def add(x, i):
        y = list(x); y[i] = (y[i] + 1) % shape[i]; return tuple(y)
    plaqs = {(x, i, j): frozenset([(x, i), (add(x, i), j), (add(x, j), i), (x, j)])
             for x in sites for i in range(n) for j in range(i + 1, n)}
    return sites, edges, ends, at, nbr, plaqs

def components(E, ends):
    E = set(E); comps = []
    while E:
        e = E.pop(); comp = {e}; stack = [e]
        while stack:
            f = stack.pop()
            for g in list(E):
                if set(ends[g]) & set(ends[f]): E.discard(g); comp.add(g); stack.append(g)
        comps.append(frozenset(comp))
    return comps

def min_degree_ok(E, ends):
    deg = Counter(v for e in E for v in ends[e])
    return all(c >= 2 for c in deg.values())

def check_A(shape, R=3):
    sites, edges, ends, at, nbr, plaqs = setup(shape)
    V = lambda E: {v for e in E for v in ends[e]}
    o = tuple(0 for _ in shape); q1 = (o, 0, 1); phi = plaqs[q1]
    # plaquettes within vertex-distance R of q1
    near = set(V(phi))
    for _ in range(R):
        near |= {v for e in edges if set(ends[e]) & near for v in ends[e]}
    Q2 = [p for p, E in plaqs.items() if V(E) <= near]
    res = Counter()
    for q2 in Q2:
        phi2 = plaqs[q2]
        for lay in ('s', 's-1'):
            if lay == 's-1' and (q2 == q1 or V(phi2) & V(phi)): continue      # distinct fronts of one layer are vertex-disjoint
            Fm = set(phi) | (set(phi2) if lay == 's-1' else set()); Fp = set(phi2) if lay == 's' else set()
            cand = [p for p, E in plaqs.items() if len(E & (phi | phi2)) >= 2]
            for nQ in (1, 2, 3):
                for Q in itertools.combinations(cand, nQ):
                    EQ = [plaqs[p] for p in Q]; cnt = Counter(e for E in EQ for e in E)
                    if any((e in Fm) + (e in Fp) + cnt[e] < 2 for e in set().union(*EQ) | phi | phi2): continue
                    if not any(V(E) & V(phi) for E in EQ) or not any(V(E) & V(phi2) for E in EQ): continue
                    nodes = [V(E) for E in EQ] + [V(phi), V(phi2)]           # connectivity (common vertices)
                    seen = {0}; st = [0]
                    while st:
                        a = st.pop()
                        for b in range(len(nodes)):
                            if b not in seen and nodes[a] & nodes[b]: seen.add(b); st.append(b)
                    if len(seen) != len(nodes): continue
                    kind = 'cap' if (nQ == 1 and q2 == q1) else ('shares a link' if len(phi & phi2) == 1 else 'shares a vertex only')
                    res[(nQ, lay, kind)] += 1
    return res

def check_B(shape, SIZES0=(4, 6)):
    sites, edges, ends, at, nbr, plaqs = setup(shape)
    V = lambda E: {v for e in E for v in ends[e]}
    o = tuple(0 for _ in shape); c = plaqs[(o, 0, 1)]; Vc = V(c)
    fronts = sorted({f for k in SIZES0 for f in fronts_of_size(shape, k, sites, edges, ends, at, nbr, Vc)}, key=lambda f: (len(f), sorted(f)))
    fronts = [f for f in fronts if V(f) & Vc]
    out = Counter(); bad = []
    # unions of vertex-disjoint fronts through c (at most 4, one per vertex of c at most)
    def unions(start, chosen, usedV):
        yield list(chosen)
        if len(chosen) == 4: return
        for i in range(start, len(fronts)):
            f = fronts[i]
            if V(f) & usedV: continue
            chosen.append(f); yield from unions(i + 1, chosen, usedV | V(f)); chosen.pop()
    cl = sorted(c)
    for U0list in unions(0, [], set()):
        U0 = set().union(*U0list) if U0list else set()
        for r in range(5):
            for D in itertools.combinations(cl, r):
                U1 = U0 ^ set(D)
                if not set(c) <= (U0 | U1): continue
                comps1 = components(U1, ends) if U1 else []
                if not all(min_degree_ok(C, ends) and (V(C) & Vc) for C in comps1): continue
                sizes = tuple(sorted(len(f) for f in U0list)), tuple(sorted(len(C) for C in comps1))
                nfp = len(U0list) + len(comps1)
                if nfp >= 3:
                    out[(nfp, sizes)] += 1
                    if nfp == 3 and all(s == 4 for s in sizes[0] + sizes[1]): bad.append((U0list, D))
    return out, bad, len(fronts)

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else 'AB'
    if 'A' in which:
        for shape in ((7, 7), (6, 6, 6), (6, 6, 6, 6)):
            res = check_A(shape, R=3)
            print(f"(A) spatial torus {shape}: admissible junctions, entry and exit plaquette boundaries, |Q| <= 3:")
            for key in sorted(res): print(f"     |Q|={key[0]}  exit layer {key[1]:4s} exit {key[2]:21s}: {res[key]}")
            print(f"     |Q|=3 with a four-edge exit: {sum(v for k, v in res.items() if k[0] == 3)}")
            sys.stdout.flush()
    if 'B' in which:
        for shape, sizes in (((7, 7), (4, 6)), ((7, 7, 7), (4,)), ((7, 7, 7), (4, 6))):
            out, bad, nf = check_B(shape, sizes)
            print(f"(B) spatial torus {shape}, fronts of layer s-1 with sizes {sizes} through the plaquette ({nf} fronts): junctions with one plaquette and >= 3 front-pairs:")
            for key in sorted(out): print(f"     {key[0]} front-pairs, sizes at layer s-1 {key[1][0]}, at layer s {key[1][1]}: {out[key]}")
            print(f"     with exactly three front-pairs of four edges: {len(bad)}")
            sys.stdout.flush()

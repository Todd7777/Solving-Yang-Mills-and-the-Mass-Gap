# Section 13 of Part IV, Lemma mg:lem:jsums6 (a compass): for every six-cycle phi through the origin of Z^2 and Z^3, count the junctions (local
# conditions as in junction_lemmas_abcd.py) with one plaquette c and the two front-pairs (phi, s-1) and (psi, layer), psi a four-cycle, either layer.
# The lemma says: 2 (the two narrowings, by c and by q) if phi = bd q xor bd c is the boundary of a domino (flat or bent), and 0 otherwise.
# Expected output: n=2: 12 six-cycles, all dominos with 2; n=3: 132 six-cycles (= N^v_6 for d=4), 108 dominos with 2, 24 skew hexagons with 0.
# Usage: python3 six_edge_narrowings.py   (about a minute)
import sys, os, itertools
here = os.path.dirname(os.path.abspath(__file__))
import importlib.util
spec = importlib.util.spec_from_file_location('jl', os.path.join(here, 'junction_lemmas_abcd.py'))
jl = importlib.util.module_from_spec(spec)
src = open(spec.origin).read().split("if __name__")[0]
exec(compile(src, spec.origin, 'exec'), jl.__dict__)
from collections import Counter
for n in (2, 3):
    R = 2; o = tuple([0]*n)
    near = jl.plaqs_near(n, R)
    cyc6 = [E for E in jl.cycles(n, 6, [o], R)]
    plaqE = {p: jl.plaq_edges(p, n) for p in near}
    tally = Counter()
    for phi in cyc6:
        Vphi = jl.verts(phi, n)
        # is phi a domino boundary? exists q ~ c with bd q xor bd c = phi
        dom = sum(1 for q, c in itertools.permutations(near, 2) if len(plaqE[q] & plaqE[c]) == 1 and plaqE[q] ^ plaqE[c] == phi)
        cnt = 0
        for c in near:
            if not (jl.verts(plaqE[c], n) & Vphi): continue
            for p in near:
                psi = plaqE[p]
                for lay in (-1, 1):
                    if jl.is_junction({c}, [(phi, -1), (psi, lay)], n): cnt += 1
        tally[(dom, cnt)] += 1
    print(f"n={n}: {len(cyc6)} six-cycles through the origin; (ordered domino decompositions, narrowing junctions) -> number of cycles: {dict(tally)}")

# Written by the independent referee of Section 13 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs. Code written from the definitions alone (no lemma of the paper is used).
# Independent check of Lemma 13.3 (mg:lem:jsums6), first bullet: for every six-edge front phi through the origin of Z^n (n = 2, 3), the junctions
# (definitions J1-J4, see lattice.py) with one plaquette c and exactly two front-pairs (phi, -1), (psi, lay), psi any front with <= KMAX edges.
# Lemma: with |psi| = 4 there are exactly two (psi = bd q, plaquette c; psi = bd c, plaquette q; other layer) if phi = bd q xor bd c (q ~ c), none
# otherwise.  We also tabulate the exits with more edges (they are counted crudely in Sigma^(6)).
import sys, time
from collections import Counter
from lattice import *

def shift(F, v):
    return frozenset(E(add(a, v), add(b, v)) for (a, b) in F)

def run(n, KMAX):
    t0 = time.time(); o = tuple([0] * n)
    F0 = fronts_through(o, KMAX, n)
    six = [f for f in F0 if len(f) == 6]
    print(f"n={n}: {len(six)} six-edge fronts through the origin; fronts with <= {KMAX} edges per vertex: {len(F0)}")
    def domino_pairs(phi):
        out = set()
        for v in verts(phi):
            for qq in plaqs_at_vertex(v, n):
                for e in qq:
                    for cc in plaqs_containing_edge(e, n):
                        if cc != qq and len(qq & cc) == 1 and (qq ^ cc) == phi: out.add(frozenset([qq, cc]))
        return out
    tal = Counter(); exits = Counter()
    for phi in six:
        dps = domino_pairs(phi)
        Vphi = verts(phi)
        cs = set()
        for v in Vphi: cs |= plaqs_at_vertex(v, n)
        n4 = 0; kinds = []
        for c in cs:
            Vc = verts(c)
            cand = set()
            for v in Vc:
                for f in F0: cand.add(shift(f, v))
            for psi in cand:
                if not (verts(psi) & Vc): continue
                for lay in (-1, 1):
                    if psi == phi and lay == -1: continue
                    if lay == -1 and (verts(psi) & Vphi): continue
                    cnt = Counter()
                    for e in phi: cnt[e] += 1
                    for e in psi: cnt[e] += 1
                    for e in c: cnt[e] += 1
                    if min(cnt.values()) < 2: continue
                    exits[(len(dps) > 0, len(psi), 'same' if lay == -1 else 'other')] += 1
                    if len(psi) == 4:
                        n4 += 1
                        ok = any(c in pr and psi == next(iter(pr - {c})) for pr in dps) and lay == 1
                        kinds.append(ok)
        tal[(len(dps), n4, all(kinds))] += 1
    print(f"  (number of domino decompositions {{q,c}}, narrowings to four edges, all of the expected form) -> number of six-cycles: {dict(tal)}")
    print(f"  all one-plaquette exits (domino?, |psi|, layer of psi) summed over the six-cycles: {dict(sorted(exits.items()))}")
    print(f"  ({time.time()-t0:.0f}s)"); sys.stdout.flush()

if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]))

#!/usr/bin/env python3
"""Lemma mg:lem:threelow and Lemma mg:lem:jsums5: enumeration from the definitions against the lemma's families.

For the torus (Z/L)^D (D = d-1), the plaquette q at the origin in the plane (0,1), and both layers lam of p = (dq, lam):
  (i)   |P| = 3, other front-pairs with 4, 4 edges;  compare with the families 'caps' and 'cube';
  (ii)  |P| = 2, other front-pairs with 4, 6 edges;  compare with (a), (b), (c);
  (iii) |P| = 2, other front-pairs with 4, 7 edges;  compare with (a), (b) for Theta = dr u dc.
Then the pairs (J, p') of these junctions are tabulated by (n, E), E = (|p'|-4) + sum |psi|/2, and by whether p' has the
layer of p (the factor eta of Lemma mg:lem:jsums4).  Also: no junction with |P| = 2 and front-pairs of 4, 4, 4 (or 4, 4, 5)
edges (Lemma mg:lem:betafour(iv)), and none with |P| = 1, three front-pairs, p with four edges and the others with at most 7
edges (Lemma mg:lem:betafour(iii)).
Usage: python3 analyze_threelow.py D L [orientation]
"""
import sys
import time
from collections import Counter
from fractions import Fraction
from lattice16b import Lattice
from fronts16b import FrontCache
from enum_junctions import enumerate_junctions, is_junction
from families16b import family_i, family_ii, near, PV


def table_E(found, p):
    tab = Counter()
    for (P, Phi) in found:
        n = len(P)
        others = [x for x in Phi if x != p]
        for pp in others:
            E = Fraction(len(pp[0]) - 4) + sum(Fraction(len(psi[0]), 2) for psi in others if psi != pp)
            same = (pp[1] == p[1])
            tab[(n, E, 'p\' in layer of p' if same else 'p\' in other layer')] += 1
    return tab


def nu2(lat, q):
    Nl = sorted(near(lat, q))
    c = 0
    for a in range(len(Nl)):
        for b in range(a + 1, len(Nl)):
            if not (PV(lat, Nl[a]) & PV(lat, Nl[b])):
                c += 1
    return c


def triples_first(lat, q):
    """#{(q, b, c): c vertex-disjoint from q, b sharing a vertex with q and c}"""
    cnt = 0
    for b in near(lat, q):
        for c in near(lat, b):
            if c != q and not (PV(lat, c) & PV(lat, q)):
                cnt += 1
    return cnt


def main():
    D, L = int(sys.argv[1]), int(sys.argv[2])
    orient = (0, 1) if len(sys.argv) < 4 else tuple(int(x) for x in sys.argv[3].split(','))
    d = D + 1
    lat = Lattice(D, L)
    fc = FrontCache(lat, [4, 5, 6, 7])
    q = lat.plaq(lat.origin(), *orient)
    Pv = 2 * (d - 1) * (d - 2)
    KS = 4 * (2 * d - 5)
    Jv = 4 * (Pv - 1) - 2 * KS
    n2 = nu2(lat, q)
    N3 = 3 * n2 + 12 * (d - 3)
    N2 = KS * KS // 2 + KS * Jv + 8 * (d - 3)
    N2p = KS * KS // 2 + KS * Jv
    print(f"D={D} (d={d}) L={L} q in plane {orient}: P_v={Pv} K_S={KS} J_v={Jv} nu_2={n2} "
          f"#triples(q,b,c)={triples_first(lat, q)} (2 nu_2 = {2 * n2})")
    print(f"  formulas: N_3 = 3 nu_2 + 12(d-3) = {N3},  N_2 = {N2},  N'_2 = {N2p}")
    allok = True
    for lam in (0, 1):
        p = (frozenset(lat.plaq_edges(q)), lam)
        for label, n, sizes, famfun in (("(i)", 3, (4, 4), lambda: family_i(lat, q, lam)),
                                        ("(ii)", 2, (4, 6), lambda: family_ii(lat, q, lam)),
                                        ("(iii)", 2, (4, 7), lambda: family_ii(lat, q, lam, seven=True))):
            t0 = time.time()
            found = enumerate_junctions(lat, fc, q, lam, n, sizes)
            fset = set(found)
            fam = famfun()
            valid = sum(1 for k in fam if is_junction(lat, k[0], k[1], p=p, nfp=3))
            missing = fset - set(fam)        # found by the search, not in the lemma's list
            extra = set(fam) - fset           # in the lemma's list, not found
            per = Counter(fam[k] for k in fam if k in fset)
            ok = (not missing) and (not extra) and valid == len(fam)
            allok &= ok
            print(f"  lam={lam} {label:5s} |P|={n} sizes={sizes}: search {len(fset):4d}, lemma's families {len(fam):4d} "
                  f"(all valid: {valid == len(fam)}), not in lemma {len(missing)}, lemma-only {len(extra)}  "
                  f"[{time.time() - t0:.1f}s]")
            print(f"        by family: {dict(sorted(per.items()))}")
            for (nn, E, lay), m in sorted(table_E(found, p).items()):
                print(f"        pairs (J,p') with (n,E)=({nn},{E}), {lay}: {m}")
        # negative checks
        for n, sizes in ((2, (4, 4)), (2, (4, 5)), (1, (4, 4)), (1, (4, 5)), (1, (4, 6)), (1, (4, 7)),
                         (1, (5, 5)), (1, (5, 6)), (1, (6, 6)), (1, (6, 7)), (1, (7, 7))):
            m = len(enumerate_junctions(lat, fc, q, lam, n, sizes))
            allok &= (m == 0)
            print(f"  lam={lam} |P|={n} other sizes {sizes}: {m} junctions")
    print("ALL CHECKS PASSED" if allok else "SOME CHECK FAILED")


if __name__ == "__main__":
    main()

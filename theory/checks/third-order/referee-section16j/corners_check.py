# Referee check (tenth pass on Section 16, Corollary mg:cor:sectorscap), item 2: the d=3 set with the empty exponents.
#  (a) Independently of orders/assemble_d3_s9z.py: the exponent set E of Lemma mg:lem:corners, the set S''_9 of Theorem mg:thm:thirdcap(c), the added set
#      Z = {(n,2): 4<=n<=40} u {(n,3),(n,7/2): 5<=n<=40}, and the minimal elements of E \ (S''_9 u Z), computed exactly on a window large enough that every
#      exponent outside it dominates one of the corners found (checked). Compared with orders/d3_s9z.json, whose terms must equal those of orders/d3_s9.json and
#      avoid Z.
#  (b) The terms of d3_s9.json: no pair at an exponent of Z (Lemma mg:lem:planar(ii) predicts none), and E = (|p'| - 4) + sum_psi |psi|/2 of each term (the class
#      written as (count, n, E, |p'|, number of further front-pairs)) consistent with E in the exponent set.
#  (c) A brute-force check of the planar face argument of Lemma mg:lem:planar(ii): for all placements of three fronts of 4,4,4 edges (unit squares), or 4,4,6
#      (two squares and a domino) or 4,4,7 (two squares and a domino with its middle link) in a window of Z^2, every bounded face of the union W of the
#      front edges is a front square or the domino (so a set P of plaquettes that is a union of bounded faces of W has at most 3, resp. 4, plaquettes).
# Usage: python3 corners_check.py
import json, os, itertools
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__)); orders = os.path.join(here, '..', 'orders')

def in_E(n, E):  # Lemma mg:lem:corners: n >= 2, E in {2, 3, 7/2, 4, 9/2, ...}, E >= 3 if n = 2
    return n >= 2 and (E == 2 or (E >= 3 and (2*E).denominator == 1)) and (n > 2 or E >= 3)

def S9(n, E): return n + E <= Fr(17, 2) or (n + E <= Fr(19, 2) and n >= 4)          # S(8) u {n+E <= 19/2, n >= 4}
def Zset(n, E): return (E == 2 and 4 <= n <= 40) or (E in (Fr(3), Fr(7, 2)) and 5 <= n <= 40)

def corners(inS, nmax=60, Emax=Fr(30)):
    """the minimal elements (componentwise order) of the exponents outside S, by brute force over the window n <= nmax, E <= Emax"""
    pts = [(n, Fr(e, 2)) for n in range(2, nmax + 1) for e in range(4, int(2*Emax) + 1) if in_E(n, Fr(e, 2))]
    out = [p for p in pts if not inS(*p)]
    mins = sorted(p for p in out if not any(q != p and q[0] <= p[0] and q[1] <= p[1] for q in out))
    for (n, E) in out:                                        # every exponent outside S (on the window) dominates a corner
        assert any(m[0] <= n and m[1] <= E for m in mins), (n, E)
    return mins

if __name__ == "__main__":
    J9 = json.load(open(os.path.join(orders, 'd3_s9.json')))
    Jz = json.load(open(os.path.join(orders, 'd3_s9z.json')))
    # (a)
    c9 = corners(lambda n, E: S9(n, E))
    cz = corners(lambda n, E: S9(n, E) or Zset(n, E))
    print("(a) corners of S''_9:      ", [(n, float(E)) for n, E in c9])
    print("    corners of S''_9 u Z:  ", [(n, float(E)) for n, E in cz])
    print("    d3_s9.json corners:    ", J9['corners'], ' agree:', sorted((n, Fr(E).limit_denominator(2)) for n, E in J9['corners']) == c9)
    print("    d3_s9z.json corners:   ", Jz['corners'], ' agree:', sorted((n, Fr(E).limit_denominator(2)) for n, E in Jz['corners']) == cz)
    # window sufficiency: every exponent with n > 60 dominates (41,2); every exponent with E > 30 dominates (2,7)
    assert (41, 2) in cz and (2, 7) in cz
    print("    window: n <= 60, E <= 30; beyond it every exponent dominates (41,2) or (2,7): True")
    print("    terms of d3_s9z.json equal those of d3_s9.json:", J9['terms'] == Jz['terms'])
    # (b)
    bad = [t for t in J9['terms'] if Zset(t[1], Fr(t[2]).limit_denominator(2))]
    print(f"(b) terms of S''_9 at an exponent of Z: {bad}  (none expected: Lemma mg:lem:planar(ii))")
    allin = all(S9(t[1], Fr(t[2]).limit_denominator(2)) and in_E(t[1], Fr(t[2]).limit_denominator(2)) for t in J9['terms'])
    print(f"    every term lies in S''_9 and in the exponent set: {allin}")
    E2 = sorted({(t[1], t[2]) for t in J9['terms'] if t[2] == 2.0}); E3 = sorted({(t[1], t[2]) for t in J9['terms'] if t[2] in (3.0, 3.5)})
    print(f"    exponents with E = 2 that carry pairs: {E2};  with E in (3, 7/2): {E3}")
    # (c) brute force of the face argument
    def sq(x, y):            # the four edges of the unit square with lower-left corner (x,y); an edge is ((x1,y1),(x2,y2)) sorted
        return {((x, y), (x + 1, y)), ((x, y + 1), (x + 1, y + 1)), ((x, y), (x, y + 1)), ((x + 1, y), (x + 1, y + 1))}
    def domino(x, y, horiz, middle):
        a, b = (x, y), ((x + 1, y) if horiz else (x, y + 1))
        e = sq(*a) ^ sq(*b)
        return (e | (sq(*a) & sq(*b))) if middle else e, {a, b}
    def bounded_faces(W, lo, hi):
        """the connected components of the unit cells of [lo, hi)^2 (cells adjacent across an edge not in W), minus the one touching the border"""
        cells = {(x, y) for x in range(lo, hi) for y in range(lo, hi)}
        seen, comps = set(), []
        for c in cells:
            if c in seen: continue
            comp, stack, border = set(), [c], False
            seen.add(c)
            while stack:
                x, y = stack.pop(); comp.add((x, y))
                if x in (lo, hi - 1) or y in (lo, hi - 1): border = True
                for (dx, dy, e) in ((1, 0, ((x + 1, y), (x + 1, y + 1))), (-1, 0, ((x, y), (x, y + 1))), (0, 1, ((x, y + 1), (x + 1, y + 1))), (0, -1, ((x, y), (x + 1, y)))):
                    nb = (x + dx, y + dy)
                    if e in W or nb not in cells or nb in seen: continue
                    seen.add(nb); stack.append(nb)
            if not border: comps.append(frozenset(comp))
        return comps
    R = 4                      # fronts placed with their lower-left corners in [0, R)^2; window [-3, R+4)
    pos = [(x, y) for x in range(R) for y in range(R)]
    n444 = n446 = n447 = 0; ok = True; maxP = {444: 0, 446: 0, 447: 0}
    for p1, p2, p3 in itertools.product(pos, repeat=3):
        W = sq(*p1) | sq(*p2) | sq(*p3); fronts = {p1, p2, p3}
        for F in bounded_faces(W, -3, R + 4):
            if not (len(F) == 1 and next(iter(F)) in fronts): ok = False; print('444 counterexample', p1, p2, p3, sorted(F))
        maxP[444] = max(maxP[444], sum(len(F) for F in bounded_faces(W, -3, R + 4))); n444 += 1
    for p1, p2, p3 in itertools.product(pos, repeat=3):
        for horiz in (True, False):
            for middle, key in ((False, 446), (True, 447)):
                De, Dsq = domino(*p3, horiz, middle)
                W = sq(*p1) | sq(*p2) | De; squares = {p1, p2} | Dsq
                faces = bounded_faces(W, -3, R + 4)
                for F in faces:
                    if not ((len(F) == 1 and next(iter(F)) in squares) or F == frozenset(Dsq)):
                        ok = False; print(key, 'counterexample', p1, p2, p3, horiz, sorted(F))
                maxP[key] = max(maxP[key], sum(len(F) for F in faces))
                if key == 446: n446 += 1
                else: n447 += 1
    print(f"(c) placements checked: {n444} (4,4,4), {n446} (4,4,6), {n447} (4,4,7); every bounded face is a front square or the domino: {ok};")
    print(f"    largest union of all bounded faces: {maxP[444]} plaquettes (4,4,4; Lemma: <= 3), {maxP[446]} (4,4,6; <= 4), {maxP[447]} (4,4,7; <= 4)")

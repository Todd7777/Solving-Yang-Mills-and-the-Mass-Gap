# Section 16 (Corollary mg:cor:sectorscap): the set S''_9 in d=3 extended by the exponents that carry no pairs by Lemma mg:lem:planar(ii) (in d=3, in infinite
# volume, every pair with E = 2 has n = 3, and every pair with E in {3, 7/2} has n <= 4).
# Usage: python3 orders/assemble_d3_s9z.py      (writes orders/d3_s9z.json from orders/d3_s9.json)
"""S''_9 u Z with Z = {(n,2): 4 <= n <= 40} u {(n,3), (n,7/2): 5 <= n <= 40}. The counts at the exponents of Z are 0 (the Theorem), so the terms are those of S''_9
(none of which lies in Z: asserted); the corners are recomputed as the minimal exponents of the exponent set of Lemma mg:lem:corners outside S''_9 u Z:
(2,7), (3,6), (5,5), (6,4) and (41,2). The corners (7,3) and (8,2) of S''_9, which set its majorant, disappear."""
import json, os
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__))
NZ = 40

def in_exp(n, E): return n >= 2 and E >= 2 and (n > 2 or E >= 3)
def in_s9(n, E): return in_exp(n, E) and (n + E <= Fr(17, 2) or (n + E <= Fr(19, 2) and n >= 4))
def in_zero(n, E): return (E == 2 and 4 <= n <= NZ) or (E in (3, Fr(7, 2)) and 5 <= n <= NZ)
def inS(n, E): return in_exp(n, E) and (in_s9(n, E) or in_zero(n, E))

if __name__ == "__main__":
    J = json.load(open(os.path.join(here, 'd3_s9.json')))
    assert not [t for t in J['terms'] if in_zero(t[1], Fr(t[2]).limit_denominator(2))], 'a term of S\'\'_9 at an exponent that the Theorem makes empty'
    EXP = [(n, Fr(e, 2)) for n in range(2, NZ + 4) for e in [4] + list(range(6, 40)) if in_exp(n, Fr(e, 2))]
    out = [(n, E) for n, E in EXP if not inS(n, E)]
    corners = sorted((n, E) for n, E in out if not any((n2, E2) != (n, E) and n2 <= n and E2 <= E for n2, E2 in out))
    json.dump({'terms': J['terms'], 'corners': [(n, float(E)) for n, E in corners]}, open(os.path.join(here, 'd3_s9z.json'), 'w'))
    print(f"{len(J['terms'])} terms (those of S''_9); corners {[(n, float(E)) for n, E in corners]}")

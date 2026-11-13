# Referee (Section 16, fourth pass): the deficit at beta_0(3) for smaller sets S than S''_9 (Lemma mg:lem:corners allows any finite S).
# Usage: python3 alt_sets.py            (about a minute per set; the authors' evaluator practical_beta.deficit)
"""S_a = S(8) u {(6,3), (6,7/2), (7,2)}: only the three pairs of order 9 whose counts are zero (K = (4,6), (4,7) with six plaquettes, K = (4,4) with
seven), which the planar face argument of the report proves; S_b = S_a u {(5,4), (5,9/2)}. The exact terms are those of orders/d3_m8.json (the order-9
pairs of S_a, S_b that are not zero are added from orders/d3_s9.json), the corners are computed from the definition (corner_analysis.corners)."""
import os, sys, json
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import corner_analysis as CA
import interp_frozen as I                       # sets up practical_beta / third_order_effective with d3_s9 and the closed polymers
from mpmath import mp, mpf, log
PB = I.PB; T = I.T

def run(name, inS):
    J9 = json.load(open(os.path.join(here, '..', 'orders', 'd3_s9.json')))
    terms = [tuple(t) for t in J9['terms'] if inS(t[1], Fr(t[2]).limit_denominator(2))]
    C = CA.corners(inS)
    T.EXACT = {3: terms}; T.CORNERS = {3: [(n, float(E)) for n, E in C]}
    D, best, g = PB.deficit(3, I.beta0(3))
    b0 = I.beta0(3)
    print(f"{name}: {len(terms)} exact terms, corners {[(n, str(E)) for n, E in C]}: D/beta^4 = {float(D/b0**4):.4f} against 1/(432 beta_0) = {float(1/(432*b0)):.4f};"
          f" with the scaling bound (1+beta_0)^3 D/beta^4 = {float((1+b0)**3*D/b0**4):.4f}; bars {best['bars']}")
    sys.stdout.flush()

if __name__ == "__main__":
    mp.dps = 60
    lo8, hi9 = Fr(17, 2), Fr(19, 2)
    zeros = {(6, Fr(3)), (6, Fr(7, 2)), (7, Fr(2))}
    run("S''_9", lambda n, E: n + E <= lo8 or (n + E <= hi9 and n >= 4))
    run("S_a = S(8) + the three zero pairs", lambda n, E: n + E <= lo8 or (n, E) in zeros)
    run("S_b = S_a + (5,4), (5,9/2)", lambda n, E: n + E <= lo8 or (n, E) in zeros | {(5, Fr(4)), (5, Fr(9, 2))})
    run("S(8)", lambda n, E: n + E <= lo8)

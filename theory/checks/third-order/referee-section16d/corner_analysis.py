# Referee (Section 16, fourth pass): which pairs of S''_9 matter at beta_0(3)?  Usage: python3 corner_analysis.py
"""With the parameters that the authors' search finds at beta_0(3) (interp_frozen.capture), print
  * the corner factors (b/bbar)^n (a/abar)^E of the corners of S''_9, and the tree-sum factor S3;
  * the contribution of every exact term of orders/d3_s9.json to delta' (relative to beta^4);
  * the majorant if S''_9 is replaced by smaller sets (all minimal elements recomputed), to show which of the order-9 counts the result depends on.
The corners of any finite S are computed here from the definition (minimal elements of the complement in the exponent set E of Lemma
mg:lem:corners), independently of orders/assemble_d3_s9.py."""
import os, sys, json
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import interp_frozen as I
from mpmath import mp, mpf, exp, log, fsum

def exponent_set(nmax=14, emax=16):
    """E of Lemma mg:lem:corners: n >= 2, E in {2, 3, 7/2, 4, 9/2, ...}, E >= 3 if n = 2"""
    Es = [Fr(2)] + [Fr(k, 2) for k in range(6, 2*emax + 1)]
    return [(n, E) for n in range(2, nmax + 1) for E in Es if not (n == 2 and E < 3)]

def corners(inS):
    out = [(n, E) for (n, E) in exponent_set() if not inS(n, E)]
    return sorted((n, E) for (n, E) in out if not any((n2, E2) != (n, E) and n2 <= n and E2 <= E for (n2, E2) in out))

if __name__ == "__main__":
    mp.dps = 60
    b0 = I.beta0(I.d); g = dict(I.T.GROUPS)['SU(3)'](b0)
    P, lo = I.capture(b0, g, mpf(5)/2)
    q = I.basic(b0, g); a, b = q['a'], q['b']; bb, ab, kind = P['bars']
    out, ex = I.evaluate(b0, g, P)
    beta4 = b0**4
    print(f"beta_0(3) = {float(b0):.8g}: b = {float(b):.5e}, a = {float(a):.5e}, bbar = {float(bb)}, abar = {float(ab)}, b/bbar = {float(b/bb):.4f}, a/abar = {float(a/ab):.4f}")
    print(f"  delta'/beta^4 = {float(P['dl']/beta4):.4f}; parts/beta^4: " + ', '.join(f"{k} {float(v/beta4):.4g}" for k, v in ex['parts'].items()))
    # S3 from the corner term
    J = json.load(open(os.path.join(here, '..', 'orders', 'd3_s9.json')))
    cn = [tuple(c) for c in J['corners']]
    vals = {c: (b/bb)**c[0]*(a/ab)**c[1] for c in cn}
    cmax = max(vals.values())
    exact = fsum(c*b**n*a**E*exp(m*P['cT'])*1 for (c, n, E, kp, m) in J['terms'] if kp == 4)
    S3 = (ex['parts']['multi'] - fsum(c*b**n*a**E*exp(m*P['cT'])*((I.T.EXACT[3] and 1) if kp == 4 else 0) for (c, n, E, kp, m) in J['terms'])) / cmax
    print("  corners of S''_9 (from the JSON) and (b/bbar)^n (a/abar)^E: " + ', '.join(f"{c}: {float(v):.3e}" for c, v in sorted(vals.items())))
    # the exact terms' contributions (weights as in lower3e: f(kp)/a^(kp-4) for kp > 4)
    xF, F, F6 = P['xF'], P['F'], P['F6']
    def kap(kp):
        if kp == 4: return mpf(1)
        if kp == 6: return 36*F6*exp(6*xF)
        return F*kp**2*exp(xF*kp)
    terms = [(c*b**n*a**mpf(E)*exp(m*P['cT'])*kap(kp), (c, n, E, kp, m)) for (c, n, E, kp, m) in J['terms']]
    tot_exact = fsum(t for t, _ in terms)
    print(f"  exact terms in all: {float(tot_exact/beta4):.4e} beta^4;  those with n + E >= 9 (order 9): {float(fsum(t for t, k in terms if k[1] + k[2] >= 9)/beta4):.3e} beta^4")
    for t, k in sorted(terms, reverse=True)[:6]: print(f"     {k}: {float(t/beta4):.4e} beta^4")
    corner_term = ex['parts']['multi'] - tot_exact
    S3 = corner_term/cmax
    print(f"  corner term {float(corner_term/beta4):.4f} beta^4 = max corner factor {float(cmax):.4e} x S3 = {float(S3):.4e}")
    # alternatives: which sets S give which majorant (the exact sums change only negligibly; the corner term is max_C(...) x S3)
    lo8 = Fr(17, 2); hi9 = Fr(19, 2)
    alts = {
        "S(8) (as in (vi))": lambda n, E: n + E <= lo8,
        "S''_9 (as in (vii))": lambda n, E: n + E <= lo8 or (n + E <= hi9 and n >= 4),
        "S(8) + the zero pairs (6,3),(6,7/2),(7,2)": lambda n, E: n + E <= lo8 or (n, E) in {(6, Fr(3)), (6, Fr(7, 2)), (7, Fr(2))},
        "S(8) + zeros + (5,4),(5,9/2)": lambda n, E: n + E <= lo8 or (n, E) in {(6, Fr(3)), (6, Fr(7, 2)), (7, Fr(2)), (5, Fr(4)), (5, Fr(9, 2))},
        "S(8) + zeros + (4,5),(4,11/2)": lambda n, E: n + E <= lo8 or (n, E) in {(6, Fr(3)), (6, Fr(7, 2)), (7, Fr(2)), (4, Fr(5)), (4, Fr(11, 2))},
    }
    for name, inS in alts.items():
        C = corners(inS)
        cm = max((b/bb)**n*(a/ab)**(mpf(E.numerator)/E.denominator) for (n, E) in C)
        arg = max(C, key=lambda c: (b/bb)**c[0]*(a/ab)**(mpf(c[1].numerator)/c[1].denominator))
        dl = P['dl'] - corner_term/q['Ubar'] + cm*S3/q['Ubar']        # delta' with the other corner term (first order: delta' = Sigma/Ubar)
        D = 8*(q['K'] + 1)*q['th'] + log(q['Ubar']/q['Ups']) + log(1 + dl)
        print(f"  {name}: corners {[(n, str(E)) for n, E in C]}; largest at {(arg[0], str(arg[1]))}: corner term {float(cm*S3/beta4):.3f} beta^4,"
              f" D/beta^4 ~ {float(D/beta4):.3f}  (beta^3/432 = {float(1/(432*b0)):.3f} beta^4)")
    print("  corners of S''_9 recomputed from the definition:", [(n, str(E)) for n, E in corners(alts["S''_9 (as in (vii))"])], " JSON:", cn)

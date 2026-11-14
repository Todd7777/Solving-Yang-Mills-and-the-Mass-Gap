# Referee (Section 16, seventh pass): the parts of delta' in the setting of Remark mg:rem:nextorder(x), term by term, with the authors' evaluator.
# Usage: python3 parts_x.py [beta]          (default beta_0(4); about 30 s)
"""Prints, in units of beta^4: the exits (one plaquette: K_S (a f(6) + b f(7)) + (4(P_v-1) - 2K_S) b f(8); Z_(2,k), Z_(3,k), Z_(4,k) exact; n >= 5 crude,
by n), the junctions with two four-edge front-pairs (exact Z_(n,4) to 12 and crude beyond, by n), and, at the chosen (bbar, abar), the factor
(b/bbar)^n (a/abar)^E of every corner and the exact terms of M_S grouped by n + E.  The formulas are those of third_order_effective.lower3e, recomputed
from the returned parameters (F, F_6, xi_F, a, b, c_J)."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..'))
from mpmath import mpf, exp, binomial, fsum, log
import practical_beta as PB
import third_order_effective as T
from junction_constants import counts

beta = mpf(sys.argv[1]) if len(sys.argv) > 1 else mpf('1.53766448e-4')
PB.load_terms(4, os.path.join(here, '..', 'orders', 'd4_s10x.json'), {10: 60, 11: 66, 12: 320, 13: 0, 14: 1316, 15: 2520, 16: 12552}, 17)
D, best, g = PB.deficit(4, beta)
d = 4; b4 = beta**4
a, b, F, F6, xF, cT = best['a'], best['b'], best['F'], best['F6'], best['xF'], best['cT']
Pv, KS, Tv, N = counts(d, 'infinite')
for k, v in T.NOVR[d].items(): N[k] = mpf(v)
for k in N:
    if k >= 6 and k not in T.NOVR[d]: N[k] = min(N[k], T.decision_bound(d - 1, k))
def f(k):
    if k == 4: return mpf(1)
    if k == 6: return 36*F6*exp(6*xF)*a**2
    return F*k**2*exp(xF*k)*a**(k - 4)
def X(n, k): return 0 if k == 5 else min(binomial(4*n + 4, k), (4*n + 4)*N[k])
print(f"beta = {float(beta):.6g}: D/beta^4 = {float(D/b4):.4f}; a = {float(a):.4e} ({float(a/beta):.5f} beta), b = {float(b):.4e}, F = {float(F):.4g},"
      f" F_6 = {float(F6):.4g}, xi_F = {float(xF):.3g}, c_J = {float(cT):.6g}, varpi = b^{float(log(best['varpi'])/log(b)):.3g}")
print("parts of delta' (beta^4):", {k: round(float(v/b4), 4) for k, v in best['parts'].items()})
one = KS*(a*f(6) + b*f(7)) + (4*(Pv - 1) - 2*KS)*b*f(8)
print(f"exits: one plaquette {float(one/b4):.4f} (K_S a f(6) {float(KS*a*f(6)/b4):.4f}, K_S b f(7) {float(KS*b*f(7)/b4):.4f},"
      f" b f(8) {float((4*(Pv-1)-2*KS)*b*f(8)/b4):.4f})")
z2 = b**2*fsum(c*f(k) for k, c in T.ZEXIT[d].items()); print(f"       Z_(2,k) exact {float(z2/b4):.4f}")
for n in (3, 4):
    zn = b**n*fsum(c*f(k) for k, c in T.ZEXITN[d][n].items()); print(f"       Z_({n},k) exact {float(zn/b4):.4f}")
tot = 0
for n in range(5, 40):
    cr = 2*binomial(4*Pv + Tv*n, n)*b**n*fsum(X(n, k)*f(k) for k in range(5, min(4*n + 5, 199)))
    tot += cr
    if n <= 8: print(f"       n = {n} crude {float(cr/b4):.4g}")
print(f"       n >= 5 crude, total {float(tot/b4):.4g}")
four_exact = fsum(c*b**n for n, c in T.ZSMALL[d].items()); nmx = max(T.ZSMALL[d])
four_crude = fsum(2*binomial(4*Pv + Tv*n, n)*(n + 1)*(2*d - 3)*b**n for n in range(nmx + 1, 40))
print(f"two four-edge front-pairs: exact n <= {nmx}: {float(four_exact/b4):.4g}; crude n > {nmx}: {float(four_crude/b4):.4g}; cubes 10a^4 (d-3): {float(10*a**4/b4):.4g}")
bars = best['bars']; print(f"majorant at (bbar, abar) = (2^{float(log(bars[0], 2)):.0f}, 2^{float(log(bars[1], 2)):.0f}):")
for n, E in T.CORNERS[d]:
    print(f"   corner ({n}, {E}): (b/bbar)^n (a/abar)^E = {float((b/bars[0])**n*(a/bars[1])**E):.3e}")
groups = {}
for (c, n, E, kp, m) in T.EXACT[d]:
    w = c*b**n*a**E*exp(m*cT)*(f(kp)/a**(kp - 4) if kp > 4 else 1)
    groups[n + E] = groups.get(n + E, 0) + w
print("exact terms of M_S by n + E (beta^4):", {k: f"{float(v/b4):.2e}" for k, v in sorted(groups.items())})

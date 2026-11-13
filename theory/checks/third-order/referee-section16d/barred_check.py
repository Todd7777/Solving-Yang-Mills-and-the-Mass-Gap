# Referee (Section 16, fourth pass): the barred constants of Lemma mg:lem:jsums4 at beta_0(3), checked against their definition.
# Usage: python3 barred_check.py
"""At (bbar, abar) the authors' tree_constants and marked return xi_F, xi_C, Ybar_j, Abar_j and gbar_F, ..., kappabar_C. Lemma mg:lem:jsums4 needs:
P_v Ybar_0 + Abar_0 <= xibar_F, 4 P_v Ybar_0 + 8 Abar_0 <= xibar_C (the first two conditions of Lemma mg:lem:jsums), and the six inequalities of
Definition mg:def:jconst with bbar in Y, abar in the front weights of A, and fbar(4) = 1, fbar(k) = f(k) (abar/a)^(k-4) in B; c_J that of the bound.
Here every quantity is recomputed from the definitions (with KMAX = 200, NMAX = 200 terms, to see the effect of the truncations at 60, 40) and the
conditions are evaluated with the authors' values; also the ratios of the series to their radii of convergence."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import interp_frozen as I
from mpmath import mp, mpf, exp, binomial, fsum

if __name__ == "__main__":
    mp.dps = 60
    b0 = I.beta0(3); g = dict(I.T.GROUPS)['SU(3)'](b0)
    P, lo = I.capture(b0, g, mpf(5)/2)
    bb, ab, kind = P['bars']; cT, xF, F, F6 = P['cT'], P['xF'], P['F'], P['F6']
    Pv, KS, Tv, N = I.counts(3, 'infinite')
    for k, v in I.T.NOVR[3].items(): N[k] = mpf(v)
    for k in range(200, 400): N[k] = 2*binomial(6*k, k - 1)
    xFb, xCb, Yb, Ab = I.tree_constants(bb, ab, cT, Pv, KS, N, 60, 40)
    def fb(k):
        if k == 4: return mpf(1)
        if k == 6: return 36*F6*exp(6*xF)*ab**2
        return F*k**2*exp(xF*k)*ab**(k - 4)
    gF, gC, g2F, g2C, kF, kC = I.marked(xFb, Yb, Ab, fb, Pv, N, 60)
    # recompute the series with many more terms
    Y = [fsum(m**j*binomial(KS*m, m - 1)*(bb*exp(xCb))**m for m in range(1, 200)) for j in range(3)]
    A = [exp(cT)*fsum(k**j*N[k]*ab**(mpf(k)/2)*exp(xFb*k) for k in range(4, 200)) for j in range(3)]
    B0 = fsum(N[k]*fb(k)*exp(xFb*k) for k in range(4, 200)); B1 = fsum(k*N[k]*fb(k)*exp(xFb*k) for k in range(4, 200))
    conds = {
        'xi_F':  (Pv*Y[0] + A[0], xFb), 'xi_C': (4*Pv*Y[0] + 8*A[0], xCb),
        'six_1': (Pv*gC*Y[1] + B0 + gF*A[1], gF), 'six_2': (4*Pv*gC*Y[1] + 8*B0 + 8*gF*A[1], gC),
        'six_3': (Pv*g2C*Y[1] + A[0] + g2F*A[1], g2F), 'six_4': (4*Pv*g2C*Y[1] + 8*A[0] + 8*g2F*A[1], g2C),
        'six_5': (gF*g2F + Pv*kC*Y[2] + g2F*B1 + gF*A[1] + kF*A[2], kF),
        'six_6': (gC*g2C + 4*Pv*kC*Y[2] + 8*(g2F*B1 + gF*A[1] + kF*A[2]), kC)}
    print(f"beta_0(3): (bbar, abar) = ({float(bb)}, {float(ab)}), c_J = {float(cT):.8f}; xibar_F = {float(xFb):.6g}, xibar_C = {float(xCb):.6g}")
    for k, (l, r) in conds.items():
        print(f"  {k}: lhs {float(l):.10e}  rhs {float(r):.10e}  relative margin {float((r - l)/r):+.3e}")
    rY = 4**4/3**3*bb*exp(xCb)                     # binom(4m, m-1)^(1/m) -> 4^4/3^3 (K_S = 4): ratio to the radius of convergence of Y
    rA = 6**6/5**5*ab**mpf(0.5)*exp(xFb)           # crude N^v_k ~ (6^6/5^5)^k for large k: ratio for the A-series (the exact N^v_k are far smaller)
    print(f"  ratio to the radius of convergence: Ybar {float(rY):.3f}, Abar (with the crude N^v_k of the tail) {float(rA):.3f}")

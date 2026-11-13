# Referee (Section 16, fourth pass): the monotonicity facts quoted in the scaling paragraph of Remark mg:rem:nextorder (vii), for SU(3), d = 3:
#   a/beta, b/beta, eps/beta, b_2/beta^2, sigma_6/beta^6, theta/beta^6 and log(Ubar/Ups)/beta^4 increase with beta, and 1 <= Ubar <= 1 + kappa (kappa = beta_0(3)).
# Usage: python3 monotone_check.py [NPOINTS]
"""The remark says these were checked at 302 points from 1e-24 to beta_0(3); no script in the repository does this (whole_range_d3.py checks only that
8(K+1)theta and log(Ubar/Ups) increase, on [1e-12, beta_0(3)]).  Here: NPOINTS log-spaced points of [1e-24, beta_0(3)], 220 digits.
Ubar - Ups is computed from the formulas of Lemma mg:lem:dressed directly (n_B^2/(Ups - c_2) + 2 e^{2 beta} sigma_6/(upsilon_1 - u^4 (1+gamma))), so that
log(Ubar/Ups) = log1p((Ubar - Ups)/Ups) has no cancellation; u and Ups come from the authors' su3 (Weyl quadrature), which is the definition used."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); up = os.path.join(here, '..')
for sub in ('second-order', 'upper-bound', 'lower-bound'): sys.path.insert(0, os.path.join(up, '..', sub))
from mpmath import mp, mpf, e, exp, expm1, log1p, sqrt, log
from second_order_constants import su3          # (this module sets mp.dps = 40 on import: the precision is set after the imports)
from constants_upper import beta0
mp.dps = 220

d = 3; K = 4*(2*d - 3); CLOSED = {10: 10, 11: 11, 12: 24, 13: 0, 14: 70}; CFROM = 15

def quantities(beta):
    g = su3(beta); N, iota, u, Ups = 3, 0, g['u'], g['Ups']
    eps = expm1(2*beta); u2 = expm1(beta) - beta; gam = expm1(2*beta) + beta**2*exp(beta)
    s6 = (2 - iota)*N**2*u**6 + u2**5*expm1(beta)
    nB = (u2/u)**2*eps + sqrt(2 - iota)*(N/u)*sqrt(u2**5*expm1(beta))
    c2 = (u2/u)**4*exp(beta)
    d1 = nB**2/(Ups - c2); ups1 = Ups + d1
    diff = d1 + 2*exp(2*beta)*s6/(ups1 - u**4*(1 + gam))          # Ubar - Ups
    x = e**2*K*eps
    th = 2*(d - 2)*e**6*s6 + x**CFROM/(1 - x) + sum(c*(e*eps)**n for n, c in CLOSED.items())
    lam = exp(2*(K + 1)*th); a = lam*u; b = lam*gam
    b2 = lam**2*(2*s6 + s6**2)/u**4
    return dict(a_b=a/beta, b_b=b/beta, eps_b=eps/beta, b2_b2=b2/beta**2, s6_b6=s6/beta**6, th_b6=th/beta**6,
                lu_b4=log1p(diff/Ups)/beta**4, Ubar=Ups + diff, u_b=u/beta)

if __name__ == "__main__":
    npts = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    B0 = beta0(d); kappa = B0
    lo, hi = log(mpf('1e-24')), log(B0)
    grid = [exp(lo + (hi - lo)*i/(npts - 1)) for i in range(npts)]
    Q = [quantities(b) for b in grid]
    names = ['a_b', 'b_b', 'eps_b', 'b2_b2', 's6_b6', 'th_b6', 'lu_b4', 'u_b']
    for nme in names:
        inc = all(p[nme] < q[nme] for p, q in zip(Q, Q[1:]))
        print(f"{nme:6s}: strictly increasing on {npts} points of [1e-24, beta_0(3)]: {inc};  from {float(Q[0][nme]):.6e} to {float(Q[-1][nme]):.6e}")
    ub = [q['Ubar'] for q in Q]
    print(f"1 <= Ubar <= 1 + kappa: {all(1 <= x <= 1 + kappa for x in ub)};  Ubar(beta_0) - 1 = {float(ub[-1] - 1):.6e} = {float((ub[-1]-1)/kappa):.4f} kappa;"
          f"  Ubar increasing: {all(p < q for p, q in zip(ub, ub[1:]))}")

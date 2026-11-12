#!/usr/bin/env python3
"""Size of the exact part e^{c_J}(2N_3 b^3 a^2 + N_2 b^2 a^3 + N'_2 b^2 a^{7/2}) of M_5 (Lemma mg:lem:jsums5), divided by beta^4, for SU(3) at the
thresholds of Remark mg:rem:practical(iii), with theta_7; c_J is bounded by 1.1 here (c_J -> 1 as beta -> 0), i.e. e^{c_J} <= 3.1."""
import sys, os
sys.path.insert(0, '..')
import toe_rmax as T
from mpmath import mpf, expm1, exp, e
for d, be, N3, N2, N2p in ((4, '1.7e-6', 876, 320, 312), (3, '1.8e-5', 48, 24, 24)):
    beta = mpf(be); g = dict(T.GROUPS)['SU(3)'](beta)
    K = 4*(2*d-3); eps = expm1(2*beta); x = e**2*K*eps
    th7 = 2*(d-2)*e**6*T.dressed(beta, g)[1] + x**8/(1-x)
    lam = exp(2*(K+1)*th7); a = lam*g['u']; b = lam*(expm1(2*beta) + beta**2*exp(beta))
    ex = exp(mpf('1.1'))*(2*N3*b**3*a**2 + N2*b**2*a**3 + N2p*b**2*a**mpf(3.5))
    print(f"SU(3) d={d} beta={be}: b={float(b):.3g} a={float(a):.3g}; exact part of M_5 / beta^4 <= {float(ex/beta**4):.2e} "
          f"(of which 2N_3 b^3a^2: {float(exp(mpf('1.1'))*2*N3*b**3*a**2/beta**4):.2e}); third-order term beta^3/432 = {float(1/(432*beta)):.4g} beta^4")

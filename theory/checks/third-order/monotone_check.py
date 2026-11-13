# Section 16 of Part IV (Remark mg:rem:nextorder (vii), the scaling; Remark mg:rem:upperexact (ii)): the monotonicity facts they use, SU(3), d = 3.
# Usage: python3 monotone_check.py [--d4] [POINTS]   (default d=3, 400 points of [1e-24, beta_0(d)], 220 digits; about 2 minutes)
"""Checks on a logarithmic grid that a/beta, b/beta, eps/beta, gamma/beta, b_2/beta^2, sigma_5/beta^5, sigma_6/beta^6, theta/beta^6 and
log(Ubar/Ups)/beta^4 increase with beta, that Ubar - a increases and is >= 1 (so that Ubar(1 + delta') - a >= 1), that Ubar(beta_0) - 1 <= beta_0/5 (it is beta_0/6 to 9 digits),
and that Ups >= 1.  Here theta is that of Remark mg:rem:nextorder (vi) (closed polymers of 10 to 14 plaquettes counted exactly), lambda =
e^{2(K+1)theta}, a = lambda u, b = lambda gamma.  Precision: second_order_constants sets 40 digits on import, at which u/beta is not resolved below
1e-23; and log(Ubar/Ups) is computed as log1p((Ubar - Ups)/Ups) with Ubar - Ups from the formula of Lemma mg:lem:dressed, since log Ubar - log Ups loses
all digits for beta below 1e-22 even at 130 digits (both points made by the fourth referee)."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
for sub in ('second-order', 'upper-bound'): sys.path.insert(0, os.path.join(here, '..', sub))
sys.path.insert(0, here)
from mpmath import mp, mpf, e, exp, expm1, log1p, sqrt, fsum
import third_order_effective as T
from constants_upper import beta0
mp.dps = 220
su3 = dict(T.GROUPS)['SU(3)']
CLOSED_D = {3: ({10: 10, 11: 11, 12: 24, 13: 0, 14: 70}, 15),      # the theta of Remark mg:rem:nextorder (vi): closed polymers counted exactly, crude from
             4: ({10: 60, 11: 66}, 12)}                              # 15 on (d=3); the theta of (ix) in d=4: 10, 11 exactly, crude from 12 on

def quantities(beta, d=3):
    CLOSED, CF = CLOSED_D[d]
    g = su3(beta); N, iota, u, Ups = g['N'], g['iota'], g['u'], g['Ups']
    eps = expm1(2*beta); u2 = expm1(beta) - beta; gam = expm1(2*beta) + beta**2*exp(beta)
    s6 = (2 - iota)*N**2*u**6 + u2**5*expm1(beta); s5 = (2 - iota)*N**2*u**5 + u2**4*expm1(beta)
    nB = (u2/u)**2*eps + sqrt(2 - iota)*(N/u)*sqrt(u2**5*expm1(beta)); c2 = (u2/u)**4*exp(beta)
    up1 = nB**2/(Ups - c2)                                   # Lemma mg:lem:dressed: Ubar = Ups + up1 + 2 e^{2 beta} s6/(Ups + up1 - u^4(1 + gamma))
    dU = up1 + 2*exp(2*beta)*s6/(Ups + up1 - u**4*(1 + gam)); Ubar = Ups + dU
    K = 4*(2*d - 3); x = e**2*K*eps
    th = 2*(d - 2)*e**6*s6 + x**CF/(1 - x) + fsum(c*(e*eps)**n for n, c in CLOSED.items())
    lam = exp(2*(K + 1)*th); a = lam*u; b = lam*gam; b2 = lam**2*(2*s6 + s6**2)/u**4
    return dict(a=a/beta, b=b/beta, eps=eps/beta, gam=gam/beta, b2=b2/beta**2, s5=s5/beta**5, s6=s6/beta**6, th=th/beta**6,
                ups=log1p(dU/Ups)/beta**4, Ubar_minus_a=Ubar - a), Ubar, Ups

if __name__ == "__main__":
    args = sys.argv[1:]; d = 3
    if '--d4' in args: d = 4; args.remove('--d4')
    P = int(args[0]) if args else 400
    B0 = beta0(d)
    span = -24 - mp.log10(B0)
    pts = sorted([mpf(10)**(-24 - j*span/(P - 1)) for j in range(P)])
    pts[-1] = B0
    rows = [quantities(p, d) for p in pts]
    ok = True
    for key in rows[0][0]:
        mono = all(r1[0][key] < r2[0][key] for r1, r2 in zip(rows, rows[1:]))
        ok = ok and mono
        print(f"d={d}: {key:13s} increasing on {P} points of [1e-24, beta_0(d)]: {mono};  from {float(rows[0][0][key]):.10g} to {float(rows[-1][0][key]):.10g}")
    low = all(r[0]['Ubar_minus_a'] >= 1 for r in rows); ups1 = all(r[2] >= 1 for r in rows)
    Ub0 = rows[-1][1]
    top = Ub0 - 1 <= B0/5
    print(f"Ubar - a >= 1 everywhere: {low};  Ups >= 1 everywhere: {ups1};  Ubar(beta_0) - 1 = {float(Ub0 - 1):.6g} = {float((Ub0 - 1)/B0):.10f} beta_0"
          f" (<= kappa/5 with kappa = beta_0: {top})")
    print("all facts hold:", ok and low and ups1 and top)

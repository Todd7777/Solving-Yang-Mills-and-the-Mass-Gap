# Referee check of Remark mg:rem:nextorder (viii): the error D_2 + log Ups of Theorem mg:thm:upper2 with the crude tails T_Q, T_p of Lemma mg:lem:twopoint
# (theta = x/(1-x), rooted count from n = 5), SU(3); thresholds where it reaches beta^3/432.  Own implementation of the formulas; SU(3) constants from upper_eval.py.
from mpmath import mp, mpf, e, exp, log, sqrt, expm1
from upper_eval import su3_consts, beta0
mp.dps = 50
def crude(beta, d):
    beta = mpf(beta); K = 4*(2*d - 3); eps = expm1(2*beta); x = e**2*K*eps; th = x/(1 - x)
    qs = e*K*exp((K + 1)*th)*eps
    TQ = (1 + eps)**6*exp((24*d - 42)*th)*(1 + mpf(1)/K)**(24*(d - 2))*qs**5/(1 - qs)
    Tp = (1 + eps)*exp((8*d - 11)*th)*(1 + mpf(1)/K)**K*qs**5/(1 - qs)
    g = su3_consts(beta); u, Ups = g['u'], g['Ups']; N = 3; tau0 = mpf(1)/2
    psi1 = sqrt(2)*3*u; a2 = (1 + Tp)*psi1 + Tp
    V2 = (1 + Tp)*(Ups + psi1**2) + Tp - max((1 - Tp)*psi1 - Tp, 0)**2
    L2 = (1 - TQ)/(1 + eps**6)*max(Ups - Tp*(1 + psi1)*sqrt(2)*N*u, 0)**2 - (1 + a2**2)*exp(16*beta)*(mpf(N)**4/tau0**2)*TQ
    return log(V2/L2) + log(Ups), TQ, Tp
for d in (3, 4):
    lo, hi = mpf('1e-10'), mpf('1e-6')
    for _ in range(60):
        mid = (lo + hi)/2
        if crude(mid, d)[0] < mid**3/432: lo = mid
        else: hi = mid
    err = crude(lo, d)[0]
    print(f"d={d}: crude-tail error reaches beta^3/432 at beta = {float(lo):.3e}; error there = {float(err/lo**5):.3e} beta^5")
b0 = beta0(3); err, TQ, Tp = crude(b0, 3)
print(f"d=3 at beta_0(3): T_Q = {float(TQ):.3e} = {float(TQ/b0**4):.3e} beta^4, T_p = {float(Tp/b0**4):.3e} beta^4")

# Referee check: which ingredients of Lemma mg:lem:uppertails matter for Remark mg:rem:upperexact (ii)?  E(beta_0(3))/beta_0^4 for SU(3), d = 3, n_0 = 10,
# (a) as in the paper; (b) the caps counted with w^5 instead of sigma_5 e^{5(K+1)theta_5}; (c) theta = x/(1-x) instead of theta_5 (with the caps by characters).
from mpmath import mp, mpf, e, exp, log, sqrt, expm1, fsum
import upper_eval as U
mp.dps = 50
b = U.beta0(3); g = U.su3_consts(b); d = 3; K = 12
def E_with(theta_kind, caps_by_chars):
    eps = expm1(2*b); x = e**2*K*eps
    th = x**5/(1 - x) if theta_kind == 'theta5' else x/(1 - x)
    w = eps*exp((K + 1)*th); q = e*K*w; u = g['u']
    s5 = 18*u**5 + (expm1(b) - b)**4*expm1(b)
    T = {}
    for S, nC, PS, Jp in (('Q', 6, 24*d - 42, 24*(d - 2)), ('A', 1, 8*d - 11, K)):
        ncap, Nn = U.COUNTS[(3, 10)][S]
        cap = ncap*s5*exp(5*(K + 1)*th) if caps_by_chars else ncap*w**5
        T[S] = (1 + eps)**nC*exp(th*PS)*(cap + fsum(c*w**n for n, c in Nn.items()) + (1 + mpf(1)/K)**Jp*q**11/(1 - q))
    a = sqrt(2)*3*u; ap = (1 + T['A'])*a + T['A']
    B = T['Q'] + eps**6 + 2*a*(1 + a)*T['A'] + (1 + ap**2)*exp(16*b)*324*T['Q']
    return T['A']*(2 + 3*a + 3*a**2) + B/(1 - B) if B < 1 else mpf('inf')
for lab, tk, cc in (('paper (theta_5, caps by characters)', 'theta5', True), ('caps counted with w^5', 'theta5', False), ('theta instead of theta_5', 'theta', True)):
    print(f"{lab:40s}: E(beta_0)/beta_0^4 = {float(E_with(tk, cc)/b**4):.4g}   (beta^3/432 = {float(1/(432*b)):.4f} beta^4)")
# (d) suggestion: the six-plaquette sets are unit cubes sharing exactly one link l with S; integrating their other eleven links gives a constant (gauge
# invariance at the ends of l) equal to the cube activity sigma_6, so they may be weighted sigma_6 e^{6(K+1)theta_5} instead of w^6.  n_0 = 10 and 12.
def E_cubes(n0key):
    eps = expm1(2*b); x = e**2*K*eps; th = x**5/(1 - x)
    w = eps*exp((K + 1)*th); q = e*K*w; u = g['u']; u2 = expm1(b) - b
    s5 = 18*u**5 + u2**4*expm1(b); s6 = 18*u**6 + u2**5*expm1(b)
    T = {}
    for S, nC, PS, Jp in (('Q', 6, 24*d - 42, 24*(d - 2)), ('A', 1, 8*d - 11, K)):
        ncap, Nn = U.COUNTS[n0key][S]; n0 = max(Nn)
        body = ncap*s5*exp(5*(K + 1)*th) + Nn[6]*s6*exp(6*(K + 1)*th) + fsum(c*w**n for n, c in Nn.items() if n != 6) + (1 + mpf(1)/K)**Jp*q**(n0 + 1)/(1 - q)
        T[S] = (1 + eps)**nC*exp(th*PS)*body
    a = sqrt(2)*3*u; ap = (1 + T['A'])*a + T['A']
    B = T['Q'] + eps**6 + 2*a*(1 + a)*T['A'] + (1 + ap**2)*exp(16*b)*324*T['Q']
    return T['A']*(2 + 3*a + 3*a**2) + B/(1 - B)
for key in ((3, 10), (3, 12)):
    print(f"cubes by sigma_6, n_0 = {key[1]:2d}                : E(beta_0)/beta_0^4 = {float(E_cubes(key)/b**4):.4g}")

# Constants of Theorem mg:thm:fronts (revised 10 November 2026, "near beta_0"), Lemma mg:lem:frontsrange and Remark mg:rem:fronts.
#
#   theta_5 = x^5/(1-x)                     (Lemma mg:lem:closed: only closed polymers, >= 5 plaquettes, have non-zero activity)
#   gamma_g = sup|g| <= gamma_2 = e^{2b}-1+b^2 e^b   (Lemma mg:lem:dup(d); rho without trivial summand, e.g. SU(2) fundamental)
#   a = e^{2(K+1)theta_5} u_*,  b = e^{2(K+1)theta_5} gamma_g
#   P_v = 2(d-1)(d-2), K_S = 4(2d-5); N_k^v from Lemma mg:lem:smallfronts:
#       infinite volume (spatial sides >= 6): N_4^v = P_v, N_5^v = 0;  finite volume (sides >= 4): N_4^v = (d-1)(2d-3), N_5^v = d-1;
#       N_k^v <= (d-1) C((4d-6)k, k-1) for k >= 6
#   eta_C = 2 P_v sum_n C(K_S n, n-1) y^n,  kappa_C = 8 P_v sum_n n C(K_S n, n-1) y^n   (y = b e^{lambda_C})
#   eta_T = 8 sum_k N_k^v e^{lambda_T k} a^k/(1-a^k),  nu = 2 sum_k k N_k^v e^{lambda_T k} a^{k-4}
#   conditions eta_C <= lambda_T, eta_T <= lambda_C, a^4(1+delta) < 1;  delta = kappa_C nu
#   lower bound  m >= 4 log(1/u_*) - 8(K+1) theta_5 - log(1+delta)
# lambda_C, lambda_T: the minimal fixed point, by iteration from lambda_C = 0 (the referee's remark E1 on the first version).
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'upper-bound'))
from mpmath import mp, mpf, e, exp, log, besseli, expm1, binomial
from constants_upper import consts, beta0
mp.dps = 30

def u_su2(b): return besseli(2, b) / besseli(1, b)

def Nkv(k, d, volume):
    n = d - 1
    if k == 4: return 2*n*(n-1) + (n if volume == 'finite' else 0)
    if k == 5: return n if volume == 'finite' else 0
    return n * binomial((4*d-6)*k, k-1)

def series(d, a, b, volume, KMAX=80, NMAX=80):
    K = 4*(2*d-3); KS = 4*(2*d-5); Pv = 2*(d-1)*(d-2)
    clus = lambda lC, wn: sum(wn(n)*binomial(KS*n, n-1)*(b*exp(lC))**n for n in range(1, NMAX))
    etaC = lambda lC: 2*Pv*clus(lC, lambda n: 1)
    kappaC = lambda lC: 8*Pv*clus(lC, lambda n: n)
    etaT = lambda lT: 8*sum(Nkv(k, d, volume)*exp(lT*k)*a**k/(1-a**k) for k in range(4, KMAX))
    nu = lambda lT: 2*sum(k*Nkv(k, d, volume)*exp(lT*k)*a**(k-4) for k in range(4, KMAX))
    return etaC, kappaC, etaT, nu

def theorem_fronts(beta, d, u, volume='infinite', gamma='gamma2'):
    K = 4*(2*d-3)
    eps = expm1(2*beta); x = e**2*K*eps; th5 = x**5/(1-x); lam = exp(2*(K+1)*th5)
    gam = expm1(2*beta) + beta**2*exp(beta) if gamma == 'gamma2' else expm1(4*beta)
    a = lam*u; b = lam*gam
    etaC, kappaC, etaT, nu = series(d, a, b, volume)
    lC = mpf(0)
    for _ in range(500):
        lT = etaC(lC); new = etaT(lT)
        if abs(new - lC) <= mpf(10)**-28*(1 + lC): lC = new; break
        lC = new
    lC *= 1 + mpf(10)**-15; lT = etaC(lC)*(1 + mpf(10)**-15)
    assert etaC(lC) <= lT and etaT(lT) <= lC
    delta = kappaC(lC)*nu(lT)
    assert a**4*(1 + delta) < 1
    return dict(th5=th5, a=a, b=b, lT=lT, lC=lC, delta=delta, loss=8*(K+1)*th5 + log(1 + delta),
                loss_sharp=4*log(e*K) + 8*(K+1)*th5)

def range_lemma(d):
    """Lemma mg:lem:frontsrange: at beta_0(d), with the worst cases u_* = eps, gamma_g = e^{4 beta}-1 and the finite-volume counts,
    the hypotheses hold with lambda_C = lambda_T = 1/10."""
    K = 4*(2*d-3); b0 = beta0(d)
    eps = expm1(2*b0); x = e**2*K*eps; th5 = x**5/(1-x); lam = exp(2*(K+1)*th5)
    a = lam*eps; b = lam*expm1(4*b0)
    etaC, kappaC, etaT, nu = series(d, a, b, 'finite')
    l = mpf(1)/10
    delta = kappaC(l)*nu(l)
    return dict(etaC=etaC(l), etaT=etaT(l), delta=delta, a4delta=a**4*(1+delta), lam=lam, x=x)

if __name__ == "__main__":
    print("Lemma mg:lem:frontsrange (beta = beta_0(d), worst case u_* = eps, gamma_g = e^{4beta}-1, finite-volume counts, lambda = 1/10):")
    for d in range(3, 13):
        r = range_lemma(d)
        ok = r['etaC'] <= mpf(1)/10 and r['etaT'] <= mpf(1)/10 and r['a4delta'] < 1
        print(f"  d={d:2d}: e^(2(K+1)theta_5)-1={float(r['lam']-1):.2e}, eta_C={float(r['etaC']):.4f}, eta_T={float(r['etaT']):.2e}, "
              f"delta={float(r['delta']):.2f}, a^4(1+delta)={float(r['a4delta']):.1e}  {'ok' if ok else 'FAILS'}")
    for d in (3, 4):
        b0 = beta0(d); K = 4*(2*d-3); Pv = 2*(d-1)*(d-2)
        print(f"\nSU(2), d={d}: beta_0={float(b0):.5e}; Theorem mg:thm:sharp loses 4log(eK)+8(K+1)theta_5 = {float(4*log(e*K)):.3f}+...")
        for volume in ('infinite', 'finite'):
            N4 = Nkv(4, d, volume)
            print(f"  {volume} volume (N_4^v={N4}, N_5^v={Nkv(5, d, volume)}); delta/beta -> 128 P_v N_4^v = {128*Pv*N4}:")
            for beta in (mpf('1e-7'), mpf('1e-6'), mpf('1e-5'), b0/2, b0*(1 - mpf(10)**-12)):
                r = theorem_fronts(beta, d, u_su2(beta), volume)
                D = consts(beta, d, 2, 1)['D']
                print(f"    beta={float(beta):.4e}: delta={float(r['delta']):.4e} (delta/beta={float(r['delta']/beta):.0f}), lambda_T={float(r['lT']):.2e}, "
                      f"lambda_C={float(r['lC']):.1e}; loss {float(r['loss']):.4f} (Thm sharp {float(r['loss_sharp']):.3f}); "
                      f"window with the upper bound {float(r['loss'] + D):.4f}; 4log(1/u)-loss = {float(4*log(1/u_su2(beta)) - r['loss']):.3f}")

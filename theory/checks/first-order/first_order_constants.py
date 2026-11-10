# Section 9 of Part IV ("The cap by characters"): the constants of Lemma mg:lem:capchar, Theorem mg:thm:caplower and
# Theorem mg:thm:upperchar, for SU(2) and SU(3) (fundamental representation), d = 3, 4, infinite and finite volume.
#   kappa(beta)   = tau0^{-1} int (Re chi_rho)^2 f_beta dU        (SU(2): 3 I_3/I_1; SU(3): Weyl integration)
#   kbar_cap      = kappa + (u2/u)^2 (1+(u2/u)^2) eps + sigma6/u^4 + sqrt2 e^beta sigma6^{1/2} + u^4 gamma_g
#                   with the proved bounds u2 <= e^beta-1-beta, sigma6 <= (2-iota)N^2 u^6 + u2^5 (e^beta-1)
#   lower bound   : Theorem mg:thm:junctions with b -> b_cap in the cap term of Sigma_4, and F(1+delta'-a) in place of
#                   F delta' (exits with >= 5 edges); delta_chi = smallest admissible delta'
#   upper bound   : D_chi = log(V_chi/L_chi); whole range: D_0 = log((1+2T_p)/L_0) >= D_chi, increasing in beta
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', 'upper-bound')); sys.path.insert(0, os.path.join(here, '..', 'lower-bound'))
from mpmath import mp, mpf, e, exp, log, besseli, expm1, binomial, sqrt, pi, cos, sin, expj, re, conj, fsum
from constants_upper import consts, beta0
from junction_constants import counts
mp.dps = 30

def su2(beta):
    u = besseli(2, beta)/besseli(1, beta)
    return dict(name='SU(2)', u=u, kappa=3*besseli(3, beta)/besseli(1, beta), N=2, iota=1)

def su3(beta, M=48):
    # trapezoidal rule on the maximal torus (exact up to exponentially small terms for these entire integrands)
    c0 = c3 = k2 = mpf(0)
    vals = []
    for i in range(M):
        for j in range(M):
            t1, t2 = 2*pi*i/M, 2*pi*j/M
            z = [expj(t1), expj(t2), expj(-t1-t2)]
            X = z[0]+z[1]+z[2]
            vd = abs((z[0]-z[1])*(z[0]-z[2])*(z[1]-z[2]))**2/6
            w = exp(beta/3*re(X))
            vals.append((vd, w, X))
    c0 = fsum(vd*w for vd, w, X in vals)/M**2
    c3 = fsum(vd*w*re(conj(X)) for vd, w, X in vals)/M**2
    kap = fsum(vd*re(X)**2*(w/c0-1) for vd, w, X in vals)/M**2/mpf('0.5')
    return dict(name='SU(3)', u=c3/(3*c0), kappa=kap, N=3, iota=0)

def kbar_cap(beta, g):
    u, kap, N, iota = g['u'], g['kappa'], g['N'], g['iota']
    eps = expm1(2*beta); u2 = expm1(beta)-beta; gam = expm1(2*beta)+beta**2*exp(beta)
    s6 = (2-iota)*N**2*u**6 + u2**5*expm1(beta)
    return kap + (u2/u)**2*(1+(u2/u)**2)*eps + s6/u**4 + sqrt(2)*exp(beta)*sqrt(s6) + u**4*gam

def upper_chi(beta, d, g):
    N, iota = g['N'], g['iota']; tau0 = mpf(1+iota)/2; c = consts(beta, d, N, iota)
    ubar = ((1+iota)*beta+beta**2*exp(beta))/(2*N**2)
    astar = N*ubar*(1+c['Tp'])+sqrt(tau0)*c['Tp']; eta2 = astar*N*ubar/tau0; kap = g['kappa']; eps = c['eps']
    tail = (1+astar**2/tau0)*exp(16*beta)*(mpf(N)**4/tau0**2)*c['TQ']
    L = (1+kap-eta2)**2*(1-c['TQ'])/(1+eps**6) - tail
    V = (1+kap)*(1+c['Tp'])+c['Tp']
    L0 = (1-eta2)**2*(1-c['TQ'])/(1+eps**6) - tail
    return dict(D=c['D'], Dchi=log(V/L), D0=log((1+2*c['Tp'])/L0), eta2=eta2)

def lower_chi(beta, d, g, volume='infinite', KMAX=60, NMAX=40):
    """smallest admissible delta' of Theorem mg:thm:caplower (rho without trivial summand: gamma_g <= gamma_2)."""
    K = 4*(2*d-3); eps = expm1(2*beta); x = e**2*K*eps; th5 = x**5/(1-x); lam = exp(2*(K+1)*th5)
    u = g['u']; a = lam*u; gam = expm1(2*beta) + beta**2*exp(beta); b = lam*gam
    bcap = lam*min(kbar_cap(beta, g), gam)
    Pv, KS, Tv, N = counts(d, volume)
    def Yf(z, p=0): return fsum(m**p*binomial(KS*m, m-1)*(b*exp(z))**m for m in range(1, NMAX))
    cT = mpf(1); xF = mpf(0); xC = mpf(0)
    for _ in range(200):
        A = exp(cT)*fsum(N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX))
        xF_n = Pv*Yf(xC) + A; xC_n = 4*Pv*Yf(xC) + 8*A
        hJ = max(a**(mpf(k)/2)*k*Pv*Yf(xC)*exp(xF*k) for k in range(4, 12))
        cT_n = exp(hJ)/(1-a**4)
        done = abs(xF_n-xF) + abs(xC_n-xC) + abs(cT_n-cT) < mpf(10)**-25
        xF, xC, cT = xF_n, xC_n, cT_n
        if done: break
    A0 = exp(cT)*fsum(N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX))
    A1 = exp(cT)*fsum(k*N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX))
    A2 = exp(cT)*fsum(k**2*N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX))
    Y0, Y1, Y2 = Yf(xC), Yf(xC, 1), Yf(xC, 2)
    def solve(F):
        f = lambda k: mpf(1) if k == 4 else F*k**2*exp(xF*k)*a**(k-4)
        B0 = fsum(N[k]*f(k)*exp(xF*k) for k in range(4, KMAX)); B1 = fsum(k*N[k]*f(k)*exp(xF*k) for k in range(4, KMAX))
        det = (1-A1)*(1-4*Pv*Y1) - (Pv*Y1)*(8*A1)
        gF = (B0*(1-4*Pv*Y1) + Pv*Y1*8*B0)/det; gC = ((1-A1)*8*B0 + 8*A1*B0)/det
        g2F = (A0*(1-4*Pv*Y1) + Pv*Y1*8*A0)/det; g2C = ((1-A1)*8*A0 + 8*A1*A0)/det
        r1 = gC*g2C + 8*(g2F*B1 + gF*A1); r2 = gF*g2F + g2F*B1 + gF*A1
        dt = (1-4*Pv*Y2)*(1-A2) - (8*A2)*(Pv*Y2)
        kC = (r1*(1-A2) + 8*A2*r2)/dt; kF = ((1-4*Pv*Y2)*r2 + Pv*Y2*r1)/dt
        def X(n, k):
            if k == 5: return mpf(4*n+4)/5 if (N[5] > 0 and n >= 4) else 0
            return min(binomial(4*n+4, k), (4*n+4)*N[k])
        parts = dict(cap=bcap, two=(4*(Pv-1)-KS)*b**2, three=48*binomial(18, 2)*(d-2)*(2*d-3)*b**3,
                     four_plus=fsum(2*binomial(4*Pv+Tv*n, n)*(n+1)*(2*d-3)*b**n for n in range(4, NMAX)),
                     big_exit=fsum(2*binomial(4*Pv+Tv*n, n)*b**n*fsum(X(n, k)*f(k) for k in range(5, min(4*n+5, 199))) for n in range(1, NMAX)),
                     three_fronts=4*Pv*exp(4*xF)*(kC*Y2 + 4*(gC*g2F + g2C*gF)*Y1 + 16*kF*Y0))
        return fsum(parts.values()), parts, Pv*(gC*Y1/5 + gF*Y0)
    F = mpf(1)
    for _ in range(200):
        S4, parts, need = solve(F); dl = S4*(1+mpf(10)**-12)
        Fn = need/(1+dl-a)*(1+mpf(10)**-9)
        done = abs(Fn-F) < mpf(10)**-20*F
        F = Fn
        if done: break
    S4, parts, need = solve(F)
    assert need < F*(1+S4-a) and a**4*(1+S4) <= 1
    return dict(delta=S4, parts=parts, b=b, bcap=bcap, a=a, F=F, loss=8*(K+1)*th5 + log(1+S4))

def main():
    for gname, gf in (('SU(2)', su2), ('SU(3)', su3)):
        for d in (3, 4):
            b0 = beta0(d)
            print(f"{gname}, d={d}:")
            for beta in (mpf('1e-7'), mpf('1e-6'), mpf('1e-5'), b0/2, b0*(1-mpf(10)**-12)):
                g = gf(beta); up = upper_chi(beta, d, g)
                out = []
                for vol in ('infinite', 'finite'):
                    r = lower_chi(beta, d, g, vol)
                    out.append(r)
                ri = out[0]
                pc = ', '.join(f'{k}={float(v/ri["b"]):.2e}' for k, v in ri['parts'].items())
                print(f"  beta={float(beta):.4e}: u={float(g['u']):.5e} kappa={float(g['kappa']):.4e} (kappa/beta={float(g['kappa']/beta):.4f}) "
                      f"kbar_cap={float(kbar_cap(beta, g)):.4e}")
                print(f"     lower: delta_chi={float(ri['delta']):.4e} (delta_chi/beta^2={float(ri['delta']/beta**2):.2f}, delta_chi/b={float(ri['delta']/ri['b']):.4f}); "
                      f"loss {float(ri['loss']):.3e} (finite volume, sides>=4: {float(out[1]['loss']):.3e}); F={float(ri['F']):.3g}")
                print(f"     parts/b: {pc}")
                print(f"     upper: D={float(up['D']):.3e}  D_chi={float(up['Dchi']):.4e} (D_chi/beta^2={float(up['Dchi']/beta**2):.3f})  D_0={float(up['D0']):.4e}")
                sys.stdout.flush()

def monotone_check(n=60):
    """b_cap(beta) = e^{2(K+1)theta_5} min(kbar_cap, gamma_2) increases on a grid of (0, beta_0(d)] (used for the whole-range claim)."""
    for gname, gf in (('SU(2)', su2), ('SU(3)', lambda bb: su3(bb, M=32))):
        for d in (3, 4):
            K = 4*(2*d-3); b0 = beta0(d); prev = mpf(0); ok = True
            for i in range(1, n+1):
                beta = b0*i/n*(1-mpf(10)**-12); g = gf(beta)
                eps = expm1(2*beta); x = e**2*K*eps; th5 = x**5/(1-x)
                bc = exp(2*(K+1)*th5)*min(kbar_cap(beta, g), expm1(2*beta)+beta**2*exp(beta))
                ok = ok and bc > prev; prev = bc
            print(f"  {gname} d={d}: b_cap increasing on {n} grid points of (0, beta_0]: {ok}")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == 'monotone': monotone_check()
    else: main()

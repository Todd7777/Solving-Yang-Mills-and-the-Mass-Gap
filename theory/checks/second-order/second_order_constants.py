# Section 10 of Part IV ("The second-order term"): the constants of Lemma mg:lem:dressed, Theorem mg:thm:upper2 and
# Theorem mg:thm:lower2, for SU(2) and SU(3) (fundamental representation) and U(1), d = 3, 4.
#   Upsilon(beta) = largest variance of a real unit plaquette function in span(chi, chi-bar) under w/c0
#                   SU(2): 1 + 3 I_3/I_1 - 4 (I_2/I_1)^2;  U(1): C-odd 1 - I_2/I_0, C-even 1 + I_2/I_0 - 2 (I_1/I_0)^2;
#                   SU(3): Weyl integration (trapezoidal rule on the maximal torus)
#   upper bound   : D_2 = log(V_2/L_2) (Theorem mg:thm:upper2); D_2 + log Upsilon = O(beta^5)
#   lower bound   : Upsilon-bar of Lemma mg:lem:dressed (with the proved bounds u2 <= e^beta-1-beta, sigma6 <= (2-iota)N^2u^6+u2^5(e^beta-1)),
#                   varpi = b, C = 1+max(1/varpi, (1+gamma_g)/Upsilon-bar), omega = C a^4 (1+gamma_g); the junction constants of
#                   Definition mg:def:jconst with c_J = e^{h_J}/(1-omega) and h_J >= C k P_v Y_0 a^{k/2} e^{xi_F k}; Sigma_4^(2) = Sigma_4 without
#                   the cap and with varpi*min(b^2, b_2) for the shifts and jumps; delta' = smallest admissible value.
#                   Loss: m - 4 log(1/u) >= -(8(K+1)theta_5 + log Upsilon-bar + log(1+delta')).
# Usage: python3 second_order_constants.py            (the table)
#        python3 second_order_constants.py asym       (delta'/beta^3 and D_2 + log Upsilon as beta -> 0)
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', 'upper-bound')); sys.path.insert(0, os.path.join(here, '..', 'lower-bound'))
from mpmath import mp, mpf, e, exp, log, besseli, expm1, binomial, sqrt, pi, expj, re, im, fsum
from constants_upper import consts, beta0
from junction_constants import counts
mp.dps = 40

def su2(beta):
    I1, I2, I3 = besseli(1, beta), besseli(2, beta), besseli(3, beta)
    u = I2/I1
    return dict(name='SU(2)', u=u, Ups=1 + 3*I3/I1 - 4*u**2, psi1=2*u, N=2, iota=1)

def u1(beta):
    I0, I1, I2 = besseli(0, beta), besseli(1, beta), besseli(2, beta)
    Um, Up = 1 - I2/I0, 1 + I2/I0 - 2*(I1/I0)**2
    # C-odd sqrt2 sin: <f phi> = 0; C-even sqrt2 cos: <f phi> = sqrt2 u
    return dict(name='U(1)', u=I1/I0, Ups=max(Um, Up), psi1=(0 if Um >= Up else sqrt(2)*I1/I0), N=1, iota=0, Upm=(Up, Um))

def su3(beta, M=48):
    vals = []
    for i in range(M):
        for j in range(M):
            t1, t2 = 2*pi*i/M, 2*pi*j/M
            z = [expj(t1), expj(t2), expj(-t1-t2)]
            X = z[0]+z[1]+z[2]
            vd = abs((z[0]-z[1])*(z[0]-z[2])*(z[1]-z[2]))**2/6
            vals.append((vd*exp(beta/3*re(X)), X))
    c0 = fsum(p for p, X in vals)
    mR = fsum(p*re(X) for p, X in vals)/c0; mI = fsum(p*im(X) for p, X in vals)/c0
    Up = 2*(fsum(p*re(X)**2 for p, X in vals)/c0 - mR**2)
    Um = 2*(fsum(p*im(X)**2 for p, X in vals)/c0 - mI**2)
    u = mR/3                                  # <f chi-bar> = <(1+f) Re chi> = N u
    return dict(name='SU(3)', u=u, Ups=max(Up, Um), psi1=(sqrt(2)*3*u if Up >= Um else 0), N=3, iota=0, Upm=(Up, Um))

def upper2(beta, d, g):
    """D_2 = log(V_2/L_2) of Theorem mg:thm:upper2."""
    N, iota, u, Ups, psi1 = g['N'], g['iota'], g['u'], g['Ups'], g['psi1']
    tau0 = mpf(1+iota)/2; c = consts(beta, d, N, iota); Tp, TQ, eps = c['Tp'], c['TQ'], c['eps']
    a2 = (1+Tp)*psi1 + Tp
    V2 = (1+Tp)*(Ups+psi1**2) + Tp - max((1-Tp)*psi1 - Tp, 0)**2
    L2 = (1-TQ)/(1+eps**6)*max(Ups - Tp*(1+psi1)*sqrt(2-iota)*N*u, 0)**2 - (1+a2**2)*exp(16*beta)*(mpf(N)**4/tau0**2)*TQ
    return dict(V2=V2, L2=L2, D2=log(V2/L2) if L2 > 0 else None)

def dressed(beta, g):
    """Upsilon-bar of Lemma mg:lem:dressed."""
    N, iota, u, Ups = g['N'], g['iota'], g['u'], g['Ups']
    eps = expm1(2*beta); u2 = g.get('u2', expm1(beta)-beta); gam = expm1(2*beta) + beta**2*exp(beta)   # u2: a sharper bound if the group data give one
    s6 = (2-iota)*N**2*u**6 + u2**5*expm1(beta)
    nB = (u2/u)**2*eps + sqrt(2-iota)*(N/u)*sqrt(u2**5*expm1(beta))
    c2 = (u2/u)**4*exp(beta); assert Ups > c2
    ups1 = Ups + nB**2/(Ups-c2)
    assert ups1 > u**4*(1+gam)
    return ups1 + 2*exp(2*beta)*s6/(ups1 - u**4*(1+gam)), s6

def lower2(beta, d, g, volume='infinite', KMAX=60, NMAX=40):
    """smallest admissible delta' of Theorem mg:thm:lower2 (rho without trivial summand: gamma_g <= gamma_2), with varpi = b."""
    K = 4*(2*d-3); eps = expm1(2*beta); x = e**2*K*eps; th5 = x**5/(1-x); lam = exp(2*(K+1)*th5)
    u = g['u']; a = lam*u; gam = expm1(2*beta) + beta**2*exp(beta); b = lam*gam
    Ubar, s6 = dressed(beta, g)
    varpi = b; C = 1 + max(1/varpi, (1+gam)/Ubar); omega = C*a**4*(1+gam)
    b2 = lam**2*(2*s6 + s6**2)/u**4
    Pv, KS, Tv, N = counts(d, volume)
    def Yf(z, p=0): return fsum(m**p*binomial(KS*m, m-1)*(b*exp(z))**m for m in range(1, NMAX))
    cT = mpf(1); xF = mpf(0); xC = mpf(0)
    for _ in range(200):
        A = exp(cT)*fsum(N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX))
        xF_n = Pv*Yf(xC) + A; xC_n = 4*Pv*Yf(xC) + 8*A
        hJ = C*max(a**(mpf(k)/2)*k*Pv*Yf(xC)*exp(xF*k) for k in range(4, 12))
        cT_n = exp(hJ)/(1-omega)
        done = abs(xF_n-xF) + abs(xC_n-xC) + abs(cT_n-cT) < mpf(10)**-30
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
        rest = dict(three=48*binomial(18, 2)*(d-2)*(2*d-3)*b**3,
                    four_plus=fsum(2*binomial(4*Pv+Tv*n, n)*(n+1)*(2*d-3)*b**n for n in range(4, NMAX)),
                    big_exit=fsum(2*binomial(4*Pv+Tv*n, n)*b**n*fsum(X(n, k)*f(k) for k in range(5, min(4*n+5, 199))) for n in range(1, NMAX)),
                    three_fronts=4*Pv*exp(4*xF)*(kC*Y2 + 4*(gC*g2F + g2C*gF)*Y1 + 16*kF*Y0))
        parts = dict(two=(4*(Pv-1)-KS)*varpi*min(b**2, b2), **rest)
        full = b + (4*(Pv-1)-KS)*b**2 + fsum(rest.values())          # Sigma_4 with the same constants (old junctions)
        return fsum(parts.values()), parts, Pv*(gC*Y1/5 + gF*Y0), full
    F = mpf(1)
    for _ in range(300):
        S, parts, need, full = solve(F); dl = S/Ubar*(1+mpf(10)**-12)
        Fn = need/(Ubar*(1+dl)-a)*(1+mpf(10)**-9)
        done = abs(Fn-F) < mpf(10)**-25*F
        F = Fn
        if done: break
    S, parts, need, full = solve(F); dl = S/Ubar*(1+mpf(10)**-12)
    # the hypotheses of Theorem mg:thm:lower2
    hyp = dict(sigma=S < Ubar*dl, exits=need < F*(Ubar*(1+dl)-a), omega=omega < 1,
               old=omega + C*a**4*(2*Ubar*(1+dl) + full)/(1-omega) <= 1,
               beyond=a**4*Ubar*(1+dl) <= varpi, growth=a**4*Ubar*(1+dl) <= 1)
    assert all(hyp.values()), hyp
    return dict(delta=dl, parts=parts, b=b, b2=b2, a=a, F=F, Ubar=Ubar, C=C, omega=omega, full=full,
                loss=8*(K+1)*th5 + log(Ubar) + log(1+dl))

GROUPS = (('SU(2)', su2), ('SU(3)', su3), ('U(1)', u1))

def main():
    for gname, gf in GROUPS:
        for d in (3, 4):
            b0 = beta0(d)
            print(f"{gname}, d={d}:")
            for beta in (mpf('1e-7'), mpf('1e-6'), mpf('1e-5'), b0/2, b0*(1-mpf(10)**-12)):
                g = gf(beta); up = upper2(beta, d, g)
                lo = lower2(beta, d, g); lof = lower2(beta, d, g, 'finite')
                mlU = -log(g['Ups'])
                pc = ', '.join(f'{k}={float(v):.2e}' for k, v in lo['parts'].items())
                print(f"  beta={float(beta):.4e}: u={float(g['u']):.6e}  -log Ups={float(mlU):+.6e} ((-log Ups)/beta^2={float(mlU/beta**2):+.5f})"
                      f"  log(Ubar/Ups)={float(log(lo['Ubar']/g['Ups'])):.2e}")
                print(f"     upper: D_2={float(up['D2']):+.6e}  D_2+log Ups={float(up['D2']+log(g['Ups'])):.3e}")
                print(f"     lower: delta'={float(lo['delta']):.4e} (delta'/beta^3={float(lo['delta']/beta**3):.4g}); "
                      f"m-4log(1/u) >= {float(-lo['loss']):+.4e} (finite volume, sides>=5: {float(-lof['loss']):+.4e}); "
                      f"C a^4={float(lo['C']*lo['a']**4):.2e}")
                print(f"     Sigma_4^(2) parts: {pc}")
                sys.stdout.flush()

def asym():
    for gname, gf in GROUPS:
        for d in (3, 4):
            print(f"{gname}, d={d}:")
            for beta in (mpf('1e-6'), mpf('1e-7'), mpf('1e-8'), mpf('1e-9')):
                g = gf(beta); up = upper2(beta, d, g); lo = lower2(beta, d, g)
                print(f"  beta={float(beta):.0e}: (D_2+log Ups)/beta^5={float((up['D2']+log(g['Ups']))/beta**5):.4g}  "
                      f"delta'/beta^3={float(lo['delta']/beta**3):.6g}  (log Ubar - log Ups)/beta^6={float(log(lo['Ubar']/g['Ups'])/beta**6):.4g}")
                sys.stdout.flush()

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == 'asym': asym()
    else: main()

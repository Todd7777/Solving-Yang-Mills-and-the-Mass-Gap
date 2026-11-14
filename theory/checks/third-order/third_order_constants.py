# Section 11 of Part IV ("The third-order term"): the constants of Theorem mg:thm:lower3 for SU(2), SU(3), U(1), d = 3, 4.
# As Theorem mg:thm:lower2 (theory/checks/second-order/second_order_constants.py), with:
#   varpi = b^2 (so the shifts and jumps give varpi*min(b^2, b_2) = O(beta^4));
#   no cap term and no term of junctions with three plaquettes and a four-edge exit (Lemma mg:lem:threeplaq: Z_{3,4} = 0);
#   junctions with >= 3 front-pairs: one plaquette (Lemma mg:lem:oneplaq) 4 P_v b [F1' A + F1 A' + F1 (e^A - 1 - A)],
#     F1 = 8 sum_k N_k f(k), F1' = 8 sum_{k>=5} N_k f(k), A = 8 e^{c_J} sum_k N_k a^{k/2}, A' = 8 e^{c_J} sum_{k>=5} N_k a^{k/2};
#     two or more plaquettes: (b/bbar)^2 times the doubly marked tree sum with bbar in place of b in the cluster weights.
# Spatial sides >= 6 (N_5 = 0).  Usage: python3 third_order_constants.py [asym]
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
for sub in ('second-order', 'upper-bound', 'lower-bound'): sys.path.insert(0, os.path.join(here, '..', sub))
from mpmath import mp, mpf, e, exp, log, expm1, binomial, fsum, sqrt
from second_order_constants import su2, su3, u1, dressed, upper2
from constants_upper import beta0
from junction_constants import counts
mp.dps = 40

INFLATE = 1 + mpf(10)**-15   # constants defined by fixed points or by truncated series are raised by this factor, which makes their defining inequalities
                             # strict (the referee of Remark mg:rem:nextorder (vi)-(vii) found them to hold only up to 1e-28 to 1e-32 otherwise)
def tree_constants(bc, a, cT, Pv, KS, N, KMAX, NMAX):
    """xi_F, xi_C, Y_j(bc), A_j for cluster weight bc (decoration factor e^{c_J} fixed)."""
    def Yf(z, p=0): return fsum(m**p*binomial(KS*m, m-1)*(bc*exp(z))**m for m in range(1, NMAX))
    xF = mpf(0); xC = mpf(0)
    for _ in range(500):
        A0 = exp(cT)*fsum(N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX))
        xF_n = Pv*Yf(xC) + A0; xC_n = 4*Pv*Yf(xC) + 8*A0
        if xC_n > 50: return None
        done = abs(xF_n-xF) + abs(xC_n-xC) < mpf(10)**-30
        xF, xC = xF_n, xC_n
        if done: break
    else:
        return None
    xF, xC = xF*INFLATE, xC*INFLATE      # the fixed point is approached from below and the series are truncated: raise by a relative 1e-15
    A = [exp(cT)*fsum(k**j*N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX)) for j in range(3)]
    Y = [Yf(xC, j) for j in range(3)]
    return xF, xC, Y, A

def marked(xF, Y, A, f, Pv, N, KMAX):
    """g_F, g_C, g'_F, g'_C, kappa_F, kappa_C solving the six inequalities with equality (None if not contractive)."""
    B0 = fsum(N[k]*f(k)*exp(xF*k) for k in range(4, KMAX)); B1 = fsum(k*N[k]*f(k)*exp(xF*k) for k in range(4, KMAX))
    Y0, Y1, Y2 = Y; A0, A1, A2 = A
    det = (1-A1)*(1-4*Pv*Y1) - (Pv*Y1)*(8*A1)
    dt = (1-4*Pv*Y2)*(1-A2) - (8*A2)*(Pv*Y2)
    if det <= 0 or dt <= 0 or A1 >= 1 or 4*Pv*Y1 >= 1 or A2 >= 1 or 4*Pv*Y2 >= 1: return None
    gF = (B0*(1-4*Pv*Y1) + Pv*Y1*8*B0)/det; gC = ((1-A1)*8*B0 + 8*A1*B0)/det
    g2F = (A0*(1-4*Pv*Y1) + Pv*Y1*8*A0)/det; g2C = ((1-A1)*8*A0 + 8*A1*A0)/det
    gF, gC, g2F, g2C = (INFLATE*v for v in (gF, gC, g2F, g2C))      # raised first: the last two inequalities are quadratic in these
    r1 = gC*g2C + 8*(g2F*B1 + gF*A1); r2 = gF*g2F + g2F*B1 + gF*A1     # sources of the kappa system, from the raised g's
    kC = (r1*(1-A2) + 8*A2*r2)/dt; kF = ((1-4*Pv*Y2)*r2 + Pv*Y2*r1)/dt
    return gF, gC, g2F, g2C, INFLATE*kF, INFLATE*kC     # all six inequalities strict, with relative room ~1e-15 (beyond the truncated tails of B_j, A_j, Y_j)

def lower3(beta, d, g, volume='infinite', KMAX=60, NMAX=40, varpi_pow=2, bbar_grid=None):
    K = 4*(2*d-3); eps = expm1(2*beta); x = e**2*K*eps; th5 = x**5/(1-x); lam = exp(2*(K+1)*th5)
    u = g['u']; a = lam*u; gam = expm1(2*beta) + beta**2*exp(beta); b = lam*gam
    Ubar, s6 = dressed(beta, g)
    varpi = b**varpi_pow; C = 1 + max(1/varpi, (1+gam)/Ubar); omega = C*a**4*(1+gam)
    b2 = lam**2*(2*s6 + s6**2)/u**4
    Pv, KS, Tv, N = counts(d, volume)
    if volume == 'finite': N[4] = Pv; N[5] = 0          # spatial sides >= 6
    # decoration constants (actual b): xi_F, xi_C, c_J with h_J >= C k P_v Y_0 a^{k/2} e^{xi_F k}
    def Yf(z, p=0): return fsum(m**p*binomial(KS*m, m-1)*(b*exp(z))**m for m in range(1, NMAX))
    cT = mpf(1); xF = mpf(0); xC = mpf(0)
    for _ in range(300):
        A0 = exp(cT)*fsum(N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX))
        xF_n = Pv*Yf(xC) + A0; xC_n = 4*Pv*Yf(xC) + 8*A0
        hJ = C*max(a**(mpf(k)/2)*k*Pv*Yf(xC)*exp(xF*k) for k in range(4, 12))
        cT_n = exp(hJ)/(1-omega)
        done = abs(xF_n-xF) + abs(xC_n-xC) + abs(cT_n-cT) < mpf(10)**-30
        xF, xC, cT = xF_n, xC_n, cT_n
        if done: break
    A = [exp(cT)*fsum(k**j*N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX)) for j in range(3)]
    Y = [Yf(xC, j) for j in range(3)]
    def X(n, k):
        if k == 5: return mpf(4*n+4)/5 if (N[5] > 0 and n >= 4) else 0
        return min(binomial(4*n+4, k), (4*n+4)*N[k])
    Aone = 8*exp(cT)*fsum(N[k]*a**(mpf(k)/2) for k in range(4, KMAX)); Aone5 = 8*exp(cT)*fsum(N[k]*a**(mpf(k)/2) for k in range(5, KMAX))
    if bbar_grid is None: bbar_grid = [mpf(2)**(-j) for j in range(3, 16)]
    def solve(F):
        f = lambda k: mpf(1) if k == 4 else F*k**2*exp(xF*k)*a**(k-4)
        mk = marked(xF, Y, A, f, Pv, N, KMAX)
        gF, gC, g2F, g2C, kF, kC = mk
        F1 = 8*fsum(N[k]*f(k) for k in range(4, KMAX)); F15 = 8*fsum(N[k]*f(k) for k in range(5, KMAX))
        best = None
        for bb in bbar_grid:                                  # the majorant: (b/bbar)^2 x (tree sum at bbar)
            if bb < b: continue
            tc = tree_constants(bb, a, cT, Pv, KS, N, KMAX, NMAX)
            if tc is None: continue
            xFb, xCb, Yb, Ab = tc
            mkb = marked(xFb, Yb, Ab, f, Pv, N, KMAX)
            if mkb is None: continue
            gFb, gCb, g2Fb, g2Cb, kFb, kCb = mkb
            S3 = 4*Pv*exp(4*xFb)*(kCb*Yb[2] + 4*(gCb*g2Fb + g2Cb*gFb)*Yb[1] + 16*kFb*Yb[0])
            val = (b/bb)**2*S3
            if best is None or val < best[0]: best = (val, bb)
        parts = dict(two=(4*(Pv-1)-KS)*varpi*min(b**2, b2),
                     four_plus=fsum(2*binomial(4*Pv+Tv*n, n)*(n+1)*(2*d-3)*b**n for n in range(4, NMAX)),
                     big_exit=fsum(2*binomial(4*Pv+Tv*n, n)*b**n*fsum(X(n, k)*f(k) for k in range(5, min(4*n+5, 199))) for n in range(1, NMAX)),
                     one_plaq=4*Pv*b*(F15*Aone + F1*Aone5 + F1*(exp(Aone)-1-Aone)),
                     majorant=best[0])
        full = b + (4*(Pv-1)-KS)*b**2 + 48*binomial(18, 2)*(d-2)*(2*d-3)*b**3 + parts['four_plus'] + parts['big_exit'] \
               + 4*Pv*exp(4*xF)*(kC*Y[2] + 4*(gC*g2F + g2C*gF)*Y[1] + 16*kF*Y[0])
        return fsum(parts.values()), parts, Pv*(gC*Y[1]/5 + gF*Y[0]), full, best[1]
    F = mpf(1)
    for _ in range(300):
        S, parts, need, full, bb = solve(F); dl = S/Ubar*(1+mpf(10)**-12)
        Fn = need/(Ubar*(1+dl)-a)*(1+mpf(10)**-9)
        done = abs(Fn-F) < mpf(10)**-25*F
        F = Fn
        if done: break
    S, parts, need, full, bb = solve(F); dl = S/Ubar*(1+mpf(10)**-12)
    hyp = dict(sigma=S < Ubar*dl, exits=need < F*(Ubar*(1+dl)-a), omega=omega < 1,
               old=omega + C*a**4*(2*Ubar*(1+dl) + full)/(1-omega) <= 1, beyond=a**4*Ubar*(1+dl) <= varpi, growth=a**4*Ubar*(1+dl) <= 1)
    assert all(hyp.values()), hyp
    return dict(delta=dl, parts=parts, b=b, a=a, F=F, Ubar=Ubar, bbar=bb, loss=8*(K+1)*th5 + log(Ubar) + log(1+dl))

GROUPS = (('SU(2)', su2), ('SU(3)', su3), ('U(1)', u1))

def main():
    for gname, gf in GROUPS:
        for d in (3, 4):
            print(f"{gname}, d={d}:")
            for beta in (mpf('1e-9'), mpf('1e-8'), mpf('1e-7'), mpf('1e-6'), mpf('1e-5')):
                g = gf(beta); lo = lower3(beta, d, g); up = upper2(beta, d, g)
                pc = ', '.join(f'{k}={float(v):.2e}' for k, v in lo['parts'].items())
                print(f"  beta={float(beta):.0e}: -log Ups={float(-log(g['Ups'])):+.6e}  upper D_2={float(up['D2']):+.6e}  lower {float(-lo['loss']):+.6e}"
                      f"  delta'/beta^4={float(lo['delta']/beta**4):.4g}  (bbar={float(lo['bbar']):.3g})")
                print(f"     Sigma_4^(3) parts: {pc}")
                sys.stdout.flush()

if __name__ == "__main__":
    main()

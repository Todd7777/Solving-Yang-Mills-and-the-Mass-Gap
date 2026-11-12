# Referee check for Part IV, Section 16, Remark mg:rem:practical (i) and (iii).
# The functions lower_chi (Section 9, first_order_constants.py) and lower3e (Sections 13-14, third_order_effective.py) hard-code
# theta_5 = x^5/(1-x) and the general counts. They are copied here with the ratio constant theta and the counts as parameters; nothing
# else is changed. theta_7 = 2(d-2) e^6 sigma_6 + x^8/(1-x) (Lemma mg:lem:theta7).
# Usage: python3 ref_numbers.py A|B|C|D        (A: losses 8(K+1)theta; B: SU(2) whole range; C: SU(3), d=4, third order; D: threshold)
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
base = os.path.join(here, '..', '..')
for sub in ('first-order', 'second-order', 'upper-bound', 'lower-bound', 'third-order'):
    sys.path.insert(0, os.path.join(base, sub))
from mpmath import mp, mpf, e, exp, log, expm1, binomial, fsum, besseli, sqrt
import first_order_constants as FO
import third_order_effective as TE
from second_order_constants import su2 as su2_2, su3 as su3_2, dressed
from third_order_constants import tree_constants, marked
from junction_constants import counts
from constants_upper import beta0
mp.dps = 60

def x_of(beta, d):
    K = 4*(2*d-3); return K, e**2*K*expm1(2*beta)

def sigma6_su2_exact(beta, jmax=40):
    I1 = besseli(1, beta)
    return fsum((2*j+1)**2*(besseli(2*j+1, beta)/I1)**6 for j in [mpf(k)/2 for k in range(1, 2*jmax)])

def sigma6_bound(beta, N, iota, u):
    u2 = expm1(beta) - beta
    return (2-iota)*N**2*u**6 + u2**5*expm1(beta)          # Lemma mg:lem:capelement(iv)

def theta5(beta, d, **kw):
    K, x = x_of(beta, d); return x**5/(1-x)

def theta7(beta, d, s6):
    K, x = x_of(beta, d); return 2*(d-2)*e**6*s6 + x**8/(1-x)

# ---------------------------------------------------------------------------------------------------------------------------- (A)
def part_A():
    print("(A) the loss 8(K+1)theta (SU(2) for sigma_6; the x^8 term dominates theta_7 at these beta)")
    for d, beta in ((4, beta0(4)*(1-mpf(10)**-12)), (3, beta0(3)*(1-mpf(10)**-12)), (4, mpf('1e-5')), (4, mpf('1e-7'))):
        K, x = x_of(beta, d); u = besseli(2, beta)/besseli(1, beta)
        s6e = sigma6_su2_exact(beta); s6b = sigma6_bound(beta, 2, 1, u)
        t5 = theta5(beta, d); t7e = theta7(beta, d, s6e); t7b = theta7(beta, d, s6b)
        print(f"  d={d} beta={float(beta):.4e}: 8(K+1)theta_5 = {float(8*(K+1)*t5):.3e};  8(K+1)theta_7 = {float(8*(K+1)*t7e):.3e} (sigma_6 exact), "
              f"{float(8*(K+1)*t7b):.3e} (sigma_6 bound);  cube part {float(8*(K+1)*2*(d-2)*e**6*s6e):.2e};  theta_7 <= theta_5: {t7b <= t5}")
    # SU(3), d = 4, beta = 1e-7: theta_5 loss relative to beta^4 (Section 16, first paragraph)
    beta = mpf('1e-7'); K, x = x_of(beta, 4)
    print(f"  d=4, beta=1e-7: 8(K+1)theta_5/beta^4 = {float(8*(K+1)*theta5(beta, 4)/beta**4):.3e}; (beta^3/432)/beta^4 = {float(1/(432*beta)):.3e}")
    # theta_7 increasing in beta on a grid (SU(2), sigma_6 by the bound used in the evaluations, and exact)
    for d in (3, 4):
        b0 = beta0(d); prev = (mpf(0), mpf(0)); ok = True
        for i in range(1, 61):
            beta = b0*i/60*(1-mpf(10)**-12); u = besseli(2, beta)/besseli(1, beta)
            cur = (theta7(beta, d, sigma6_bound(beta, 2, 1, u)), theta7(beta, d, sigma6_su2_exact(beta)))
            ok = ok and cur[0] > prev[0] and cur[1] > prev[1]; prev = cur
        print(f"  d={d}: theta_7 (SU(2), both sigma_6) increasing on 60 grid points of (0, beta_0]: {ok}")

# ---------------------------------------------------------------------------------------------------------------------------- (B)
def lower_chi_theta(beta, d, g, th, volume='infinite', KMAX=60, NMAX=40):
    """copy of first_order_constants.lower_chi with th5 replaced by the argument th"""
    K = 4*(2*d-3); lam = exp(2*(K+1)*th)
    u = g['u']; a = lam*u; gam = expm1(2*beta) + beta**2*exp(beta); b = lam*gam
    bcap = lam*min(FO.kbar_cap(beta, g), gam)
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
    return dict(delta=S4, loss=8*(K+1)*th + log(1+S4), bcap=bcap, a=a, b=b)

def part_B():
    print("(B) SU(2): lower bound of Theorem mg:thm:caplower, infinite volume, loss = 8(K+1)theta + log(1+delta_chi)")
    for d in (4, 3):
        b0 = beta0(d)
        for beta in (b0*(1-mpf(10)**-12), b0/2, mpf('1e-5')):
            g = FO.su2(beta); u = g['u']
            t5 = theta5(beta, d); t7 = theta7(beta, d, sigma6_bound(beta, 2, 1, u)); t7e = theta7(beta, d, sigma6_su2_exact(beta))
            r5 = lower_chi_theta(beta, d, g, t5); r7 = lower_chi_theta(beta, d, g, t7); r7e = lower_chi_theta(beta, d, g, t7e)
            print(f"  d={d} beta={float(beta):.4e}: loss with theta_5 {float(r5['loss']):.4e}; with theta_7 {float(r7['loss']):.4e} "
                  f"(sigma_6 exact: {float(r7e['loss']):.4e}); delta_chi/beta^2 = {float(r7['delta']/beta**2):.1f}")
        # monotonicity used for "the whole range": loss and b_cap increase with beta on a grid
        ok_loss = ok_bcap = True; prev_loss = prev_bcap = mpf(0); worst = mpf(0)
        for i in range(1, 21):
            beta = b0*i/20*(1-mpf(10)**-12); g = FO.su2(beta); t7 = theta7(beta, d, sigma6_bound(beta, 2, 1, g['u']))
            r = lower_chi_theta(beta, d, g, t7)
            ok_loss = ok_loss and r['loss'] > prev_loss; ok_bcap = ok_bcap and r['bcap'] > prev_bcap
            prev_loss, prev_bcap = r['loss'], r['bcap']; worst = max(worst, r['loss'])
        print(f"  d={d}: with theta_7, loss and b_cap increasing on 20 grid points of (0, beta_0]: {ok_loss}, {ok_bcap}; max loss {float(worst):.4e}")

# ---------------------------------------------------------------------------------------------------------------------------- (C)
def lower3e_mod(beta, d, g, th, Nover=None, Zexact=None, Zfull=False, KMAX=60, NMAX=40, varpi_pow=2, chars=True, volume='infinite', M5=None):
    """copy of third_order_effective.lower3e (sharp6 = uturn = seven = True) with theta, N^v_k (Nover) and Z_{n,4} (Zexact: {n: count},
    used for 3 <= n <= max key in the four-edge exit terms; Zfull: also in the full sum Sigma_4 that enters eta) as parameters"""
    K = 4*(2*d-3); lam = exp(2*(K+1)*th)
    u = g['u']; a = lam*u; gam = expm1(2*beta) + beta**2*exp(beta); b = lam*gam
    Ubar, s6 = dressed(beta, g)
    varpi = b**varpi_pow; C = 1 + max(1/varpi, (1+gam)/Ubar); omega = C*a**4*(1+gam)
    b2 = lam**2*(2*s6 + s6**2)/u**4
    Pv, KS, Tv, N = counts(d, volume)
    if Nover:
        for k, v in Nover.items(): N[k] = mpf(v)
    nZ = max(Zexact) if Zexact else 0
    def crude4(n): return 2*binomial(4*Pv+Tv*n, n)*(n+1)*(2*d-3)
    def four_count(n): return mpf(Zexact[n]) if (Zexact and n <= nZ) else crude4(n)
    def Yf(z, p=0): return fsum(m**p*binomial(KS*m, m-1)*(b*exp(z))**m for m in range(1, NMAX))
    cT = mpf(1); xF = mpf(0); xC = mpf(0)
    if omega >= 1: raise AssertionError('omega >= 1')
    for _ in range(300):
        if cT > 50 or xF > 50 or xC > 50: raise AssertionError('decoration constants diverge')
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
    Ak = lambda kmin: 8*exp(cT)*fsum(N[k]*a**(mpf(k)/2) for k in range(kmin, KMAX))
    bbar_grid = [mpf(2)**(-j) for j in range(3, 18)]; abar_grid = [mpf(2)**(-j) for j in range(2, 26)]
    tcache = {}
    for bb in bbar_grid:
        if bb < b: continue
        for ab in abar_grid:
            if ab < a: continue
            tc = tree_constants(bb, ab, cT, Pv, KS, N, KMAX, NMAX)
            if tc is not None: tcache[(bb, ab)] = tc
    tcache11 = {}
    for bb in bbar_grid:
        if bb < b: continue
        tc = tree_constants(bb, a, cT, Pv, KS, N, KMAX, NMAX)
        if tc is not None: tcache11[bb] = tc
    def solve(F, F6, dl):
        def f(k):
            if k == 4: return mpf(1)
            if k == 6: return 36*F6*exp(6*xF)*a**2
            return F*k**2*exp(xF*k)*a**(k-4)
        gF, gC, g2F, g2C, kF, kC = marked(xF, Y, A, f, Pv, N, KMAX)
        Fk = lambda kmin: 8*fsum(N[k]*f(k) for k in range(kmin, KMAX))
        best = None
        for (bb, ab), tc in tcache.items():
            xFb, xCb, Yb, Ab = tc
            def fb(k, ab=ab):
                if k == 4: return mpf(1)
                if k == 6: return 36*F6*exp(6*xF)*ab**2
                return F*k**2*exp(xF*k)*ab**(k-4)
            mkb = marked(xFb, Yb, Ab, fb, Pv, N, KMAX)
            if mkb is None: continue
            gFb, gCb, g2Fb, g2Cb, kFb, kCb = mkb
            S3 = 4*Pv*exp(4*xFb)*(kCb*Yb[2] + 4*(gCb*g2Fb + g2Cb*gFb)*Yb[1] + 16*kFb*Yb[0])
            val = ((b/bb)**3*(a/ab)**2 + (b/bb)**2*(a/ab)**3)*S3
            if best is None or val < best[0]: best = (val, (bb, ab))
            if M5 is not None:                  # Lemma mg:lem:jsums5: exact lowest orders + R x tree sum
                N3, N2, N2p = M5
                R = (b/bb)**4*(a/ab)**2 + (b/bb)**3*(a/ab)**3 + (b/bb)**2*(a/ab)**4
                val5 = exp(cT)*(2*N3*b**3*a**2 + N2*b**2*a**3 + N2p*b**2*a**mpf(3.5)) + R*S3
                if val5 < best[0]: best = (val5, (bb, ab, 'M5'))
        for bb, tc in tcache11.items():
            xFb, xCb, Yb, Ab = tc
            mkb = marked(xFb, Yb, Ab, f, Pv, N, KMAX)
            if mkb is None: continue
            gFb, gCb, g2Fb, g2Cb, kFb, kCb = mkb
            S3 = 4*Pv*exp(4*xFb)*(kCb*Yb[2] + 4*(gCb*g2Fb + g2Cb*gFb)*Yb[1] + 16*kFb*Yb[0])
            val = (b/bb)**2*S3
            if best is None or val < best[0]: best = (val, (bb, 'a'))
        if Zfull:
            full = b + (4*(Pv-1)-KS)*b**2 + four_count(3)*b**3 + fsum(four_count(n)*b**n for n in range(4, NMAX))
        else:
            full = b + (4*(Pv-1)-KS)*b**2 + 48*binomial(18, 2)*(d-2)*(2*d-3)*b**3 + fsum(crude4(n)*b**n for n in range(4, NMAX))
        full += fsum(2*binomial(4*Pv+Tv*n, n)*b**n*fsum(X(n, k)*f(k) for k in range(5, min(4*n+5, 199))) for n in range(1, NMAX)) \
                + 4*Pv*exp(4*xF)*(kC*Y[2] + 4*(gC*g2F + g2C*gF)*Y[1] + 16*kF*Y[0])
        eta = C*a**4*(2*Ubar*(1+dl) + full)/(1-omega)
        parts = dict(two=(4*(Pv-1)-KS)*varpi*min(b**2, b2),
                     four=(10*(a**4 if chars else b**4) + 2*eta*b**4)*(d-3) + fsum(four_count(n)*b**n for n in range(5, NMAX)),
                     exits=KS*((a if chars else b)*f(6) + b*f(7)) + (4*(Pv-1)-2*KS)*b*f(8)
                           + fsum(2*binomial(4*Pv+Tv*n, n)*b**n*fsum(X(n, k)*f(k) for k in range(5, min(4*n+5, 199))) for n in range(2, NMAX)),
                     one_plaq=4*Pv*b*(Fk(8)*Ak(4) + Fk(4)*Ak(8) + Fk(4)*TE.em1x(Ak(4))),
                     multi=best[0])
        X6 = lambda n, k: 0 if k == 5 else min(binomial(4*n+6, k), (4*n+6)*N[k])
        sig6 = 2*(a if chars else b) + fsum(2*(2*d-3)*(mpf(2*n+3)/2)*binomial(6*Pv+Tv*n, n)*b**n for n in range(2, NMAX)) \
               + fsum(2*binomial(6*Pv+Tv*n, n)*b**n*fsum(X6(n, k)*f(k) for k in range(5, min(4*n+7, 199))) for n in range(1, NMAX)) \
               + 6*Pv*exp(6*xF)*(kC*Y[2] + 6*(gC*g2F + g2C*gF)*Y[1] + 36*kF*Y[0])
        return fsum(parts.values()), parts, Pv*(gC*Y[1]/7 + gF*Y[0]), full, best[1], sig6
    F = mpf(1); F6 = mpf(1); dl = mpf(0)
    for _ in range(100):
        S, parts, need, full, bb, sig6 = solve(F, F6, dl); dln = S/Ubar*(1+mpf(10)**-12)
        Fn = need/(Ubar*(1+dln)-a)*(1+mpf(10)**-9)
        F6n = sig6/(36*exp(6*xF)*(Ubar*(1+dln)-a))*(1+mpf(10)**-9)
        done = abs(Fn-F) < mpf(10)**-12*F and abs(F6n-F6) < mpf(10)**-12*F6 and abs(dln-dl) < mpf(10)**-12*dln
        F, F6, dl = Fn, F6n, dln
        if done: break
    dl = dl*(1+mpf(10)**-9)
    S, parts, need, full, bb, sig6 = solve(F, F6, dl)
    hyp = dict(sigma=S < Ubar*dl, exits=need < F*(Ubar*(1+dl)-a), omega=omega < 1, F6=F6 > 0,
               six=sig6 < 36*F6*exp(6*xF)*(Ubar*(1+dl)-a),
               old=omega + C*a**4*(2*Ubar*(1+dl) + full)/(1-omega) <= 1, beyond=a**4*Ubar*(1+dl) <= varpi, growth=a**4*Ubar*(1+dl) <= 1)
    assert all(hyp.values()), hyp
    return dict(delta=dl, parts=parts, Ubar=Ubar, loss=8*(K+1)*th + log(Ubar) + log(1+dl))

def best3(beta, d, g, th, **kw):
    out = None
    for vp in (2, mpf(5)/2, 3):
        try:
            lo = lower3e_mod(beta, d, g, th, varpi_pow=vp, **kw); lo['vp'] = vp
        except AssertionError:
            continue
        if out is None or lo['delta'] < out['delta']: out = lo
    return out

NEX4 = {6: 132, 7: 108}                                      # Remark mg:rem:practical(ii), d = 4 (confirmed by ref_fronts.py)
NEX4_MORE = {6: 132, 7: 108, 8: 2034, 9: 4344}              # with N^v_8, N^v_9 of Z^3 from ref_fronts.py (not used in the paper)
ZEX4 = {3: 0, 4: 12, 5: 26, 6: 302, 7: 288}                 # Remark mg:rem:practical(ii), d = 4 (n = 3..7)

def part_C(betas=('1e-5', '1e-6', '1e-7')):
    print("(C) SU(3), d=4: smallest admissible delta' of Theorem mg:thm:lower5 (chars), best of varpi = b^2, b^(5/2), b^3")
    for bs in betas:
        beta = mpf(bs); g = su3_2(beta); u = g['u']
        t5 = theta5(beta, 4); t7 = theta7(beta, 4, sigma6_bound(beta, 3, 0, u))
        old = best3(beta, 4, g, t5)
        new = best3(beta, 4, g, t7, Nover=NEX4, Zexact=ZEX4)
        newf = best3(beta, 4, g, t7, Nover=NEX4, Zexact=ZEX4, Zfull=True)
        more = best3(beta, 4, g, t7, Nover=NEX4_MORE, Zexact=ZEX4)
        thz = best3(beta, 4, g, t7)                           # theta_7 only, general counts
        pc = ', '.join(f"{k}={float(v/beta**4):.3g}" for k, v in new['parts'].items())
        print(f"  beta={bs}: delta'/beta^4: theta_5 + general counts {float(old['delta']/beta**4):.3g};  theta_7 only {float(thz['delta']/beta**4):.3g};  "
              f"theta_7 + exact counts {float(new['delta']/beta**4):.3g} (varpi=b^{float(new['vp']):g}; with Z also in Sigma_4: {float(newf['delta']/beta**4):.3g}; "
              f"with N^v_8, N^v_9 exact too: {float(more['delta']/beta**4):.3g})")
        print(f"     parts/beta^4 (theta_7 + exact counts): {pc}")
        print(f"     loss before: 8(K+1)theta_5 = {float(8*21*t5/beta**4):.3g} beta^4, total {float(old['loss'] - log(g['Ups']))/float(beta**4):.3g} beta^4;"
              f"  after: total {float(new['loss'] - log(g['Ups']))/float(beta**4):.3g} beta^4;  beta^3/432 = {float(1/(432*beta)):.3g} beta^4")
        sys.stdout.flush()

def part_D():
    print("(D) threshold: total loss (8(K+1)theta + log Ubar - log Ups + log(1+delta')) against beta^3/432, SU(3), d = 4")
    for bs in ('2e-8', '3e-8', '4e-8', '5e-8', '6e-8'):
        beta = mpf(bs); g = su3_2(beta); t7 = theta7(beta, 4, sigma6_bound(beta, 3, 0, g['u']))
        new = best3(beta, 4, g, t7, Nover=NEX4, Zexact=ZEX4)
        L = new['loss'] - log(g['Ups'])
        print(f"  beta={bs}: loss = {float(L):.3e}, beta^3/432 = {float(beta**3/432):.3e}, ratio {float(L/(beta**3/432)):.3f}")
        sys.stdout.flush()
    for bs in ('1e-9', '2e-9', '3e-9'):
        beta = mpf(bs); g = su3_2(beta); t5 = theta5(beta, 4)
        old = best3(beta, 4, g, t5)
        L = old['loss'] - log(g['Ups'])
        print(f"  before (theta_5, general counts), beta={bs}: loss = {float(L):.3e}, beta^3/432 = {float(beta**3/432):.3e}, ratio {float(L/(beta**3/432)):.3f}")
        sys.stdout.flush()

NEX3 = {6: 12, 7: 12, 8: 70, 9: 160}                        # Remark mg:rem:practical(ii), d = 3 (confirmed by ref_fronts.py)
ZEX3 = {3: 0, 4: 0, 5: 0, 6: 0, 7: 0}

def part_E():
    print("(E) SU(3), d=3: delta' of Theorem mg:thm:lower5 (chars) with theta_7 and the exact counts, against beta^3/432")
    for bs in ('1e-5', '3e-6', '1e-6'):
        beta = mpf(bs); g = su3_2(beta); t7 = theta7(beta, 3, sigma6_bound(beta, 3, 0, g['u']))
        new = best3(beta, 3, g, t7, Nover=NEX3, Zexact=ZEX3)
        old = best3(beta, 3, g, theta5(beta, 3))
        L = new['loss'] - log(g['Ups'])
        pc = ', '.join(f"{k}={float(v/beta**4):.3g}" for k, v in new['parts'].items())
        print(f"  beta={bs}: delta'/beta^4 = {float(new['delta']/beta**4):.3g} (varpi=b^{float(new['vp']):g}; before: {float(old['delta']/beta**4):.3g}, "
              f"total loss before {float((old['loss'] - log(g['Ups']))/beta**4):.3g} beta^4); total loss/(beta^3/432) = {float(L/(beta**3/432)):.3f}")
        print(f"     parts/beta^4: {pc}")
        sys.stdout.flush()

NEX4_8 = {6: 132, 7: 108, 8: 2034}                           # the counts listed in the current Remark mg:rem:practical(ii) for d = 4

def part_F():
    print("(F) SU(3), d=4: delta' with theta_7, Z_{n,4} and N^v_k for k <= 8 (the counts now listed in (ii)), and threshold with N^v_8, N^v_9")
    for bs in ('1e-5', '1e-7'):
        beta = mpf(bs); g = su3_2(beta); t7 = theta7(beta, 4, sigma6_bound(beta, 3, 0, g['u']))
        r = best3(beta, 4, g, t7, Nover=NEX4_8, Zexact=ZEX4)
        print(f"  beta={bs}: delta'/beta^4 = {float(r['delta']/beta**4):.3g} (varpi=b^{float(r['vp']):g})")
        sys.stdout.flush()
    for bs in ('8e-8', '9e-8', '1e-7'):
        beta = mpf(bs); g = su3_2(beta); t7 = theta7(beta, 4, sigma6_bound(beta, 3, 0, g['u']))
        r = best3(beta, 4, g, t7, Nover=NEX4_MORE, Zexact=ZEX4)
        L = r['loss'] - log(g['Ups'])
        print(f"  with N^v_8, N^v_9 exact, beta={bs}: delta'/beta^4 = {float(r['delta']/beta**4):.3g}; total loss/(beta^3/432) = {float(L/(beta**3/432)):.3f}")
        sys.stdout.flush()

M5_4 = (876, 320, 312)                                         # Lemma mg:lem:threelow, d = 4 (confirmed by ref_three.py)
M5_3 = (48, 24, 24)

def part_G():
    print("(G) SU(3): delta' with theta_7, exact counts (N^v_k, k <= 7 resp. 9; Z_{n,4}) and M_5 of Lemma mg:lem:jsums5")
    for d, bs in ((4, '1e-5'), (4, '1e-6'), (4, '1e-7')):
        beta = mpf(bs); g = su3_2(beta); t7 = theta7(beta, d, sigma6_bound(beta, 3, 0, g['u']))
        r = best3(beta, d, g, t7, Nover=NEX4, Zexact=ZEX4, M5=M5_4)
        r8 = best3(beta, d, g, t7, Nover=NEX4_8, Zexact=ZEX4, M5=M5_4)
        L = r['loss'] - log(g['Ups']); L8 = r8['loss'] - log(g['Ups'])
        pc = ', '.join(f"{k}={float(v/beta**4):.3g}" for k, v in r['parts'].items())
        print(f"  d={d} beta={bs}: delta'/beta^4 = {float(r['delta']/beta**4):.3g} (with N^v_8: {float(r8['delta']/beta**4):.3g}); "
              f"loss/(beta^3/432) = {float(L/(beta**3/432)):.3g} (with N^v_8: {float(L8/(beta**3/432)):.3g})")
        print(f"     parts/beta^4: {pc}")
        sys.stdout.flush()

def part_H():
    print("(H) SU(3): Remark 16.5(iii) as revised (theta_7, counts of (ii), M_5): threshold points")
    for d, bs in ((4, '1.4e-6'), (3, '1.7e-5'), (3, '1e-5'), (3, '1e-4')):
        beta = mpf(bs); g = su3_2(beta); t7 = theta7(beta, d, sigma6_bound(beta, 3, 0, g['u']))
        if d == 4: r = best3(beta, d, g, t7, Nover=NEX4_8, Zexact=ZEX4, M5=M5_4)
        else: r = best3(beta, d, g, t7, Nover=NEX3, Zexact=ZEX3, M5=M5_3)
        L = r['loss'] - log(g['Ups'])
        pc = ', '.join(f"{k}={float(v/beta**4):.3g}" for k, v in r['parts'].items())
        print(f"  d={d} beta={bs}: delta'/beta^4 = {float(r['delta']/beta**4):.3g}; (beta^3/432)/beta^4 = {float(1/(432*beta)):.3g}; "
              f"total loss/(beta^3/432) = {float(L/(beta**3/432)):.3g}")
        print(f"     parts/beta^4: {pc}")
        sys.stdout.flush()

def part_I():
    print("(I) SU(3), d=4: as (H) with also N^v_9 = 4344 (ref_fronts.py), at beta = 1e-5 and 1.4e-6")
    for bs in ('1e-5', '1.4e-6'):
        beta = mpf(bs); g = su3_2(beta); t7 = theta7(beta, 4, sigma6_bound(beta, 3, 0, g['u']))
        r = best3(beta, 4, g, t7, Nover=NEX4_MORE, Zexact=ZEX4, M5=M5_4)
        L = r['loss'] - log(g['Ups'])
        pc = ', '.join(f"{k}={float(v/beta**4):.3g}" for k, v in r['parts'].items())
        print(f"  beta={bs}: delta'/beta^4 = {float(r['delta']/beta**4):.3g}; total loss/(beta^3/432) = {float(L/(beta**3/432)):.3g}")
        print(f"     parts/beta^4: {pc}")
        sys.stdout.flush()

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == 'I': part_I(); sys.exit()
    if len(sys.argv) > 1 and sys.argv[1] == 'G': part_G(); sys.exit()
    if len(sys.argv) > 1 and sys.argv[1] == 'H': part_H(); sys.exit()
    what = sys.argv[1] if len(sys.argv) > 1 else 'ABCDEF'
    if 'E' in what: part_E()
    if 'F' in what: part_F()
    if 'A' in what: part_A()
    if 'B' in what: part_B()
    if 'C' in what: part_C()
    if 'D' in what: part_D()

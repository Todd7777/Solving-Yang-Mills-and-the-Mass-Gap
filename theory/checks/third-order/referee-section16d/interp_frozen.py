# Referee (Section 16, fourth pass): the interpolation argument of Remark mg:rem:nextorder (vii).
# Usage: python3 interp_frozen.py BETA_HI [BETA ...]       (BETA_HI may be 'b0' for beta_0(3))
"""Freeze the parameters that lower3e finds at beta_hi = beta_{i+1} (varpi as a number, C_varpi, the majorant (bbar, abar, kind), F, F_6, delta',
the decoration constants xi_F, xi_C, h_J, c_J, and the marked constants g_F, ..., kappa_C), and evaluate at each beta < beta_hi, from the
formulas of Theorems mg:thm:lower4/lower5, Lemmas mg:lem:jsums, mg:lem:jsums4, mg:lem:jsums6, mg:lem:corners and Definition mg:def:jconst:

  * the conditions of Lemma mg:lem:jsums with the frozen constants: P_v Y_0 + A_0 <= xi_F, 4 P_v Y_0 + 8 A_0 <= xi_C,
    C_varpi k P_v Y_0 a^{k/2} e^{xi_F k} <= h_J for k = 4..KMAX-1, and e^{h_J}/(1 - omega) <= c_J;
  * the six inequalities of Definition mg:def:jconst with the FROZEN g, kappa (and, separately, the minimal solution at beta);
  * the hypotheses of Theorem mg:thm:lower4/lower5 as lower3e checks them ('sigma', 'exits', 'six', 'old', 'beyond', 'growth', 'omega');
  * the deficit with the frozen delta': D_frozen(beta) = 8(K+1)theta(beta) + log(Ubar/Ups)(beta) + log(1 + delta'_hi), against D(beta_hi).

The beta-dependent quantities (theta, lambda, a, b, b_2, Ubar, Ups, gamma) are recomputed at each beta exactly as in the authors' lower3e, whose
code is mirrored below (with the same globals: practical_beta.load_terms with orders/d3_s9.json, the closed polymers 10..14, crude from 15).
Printed: the margins (right side minus left side, relative) of every condition."""
import os, sys, time
here = os.path.dirname(os.path.abspath(__file__)); up = os.path.join(here, '..')
sys.path.insert(0, up)
for sub in ('second-order', 'upper-bound', 'lower-bound'): sys.path.insert(0, os.path.join(up, '..', sub))
from mpmath import mp, mpf, e, exp, log, expm1, binomial, fsum
import practical_beta as PB
import third_order_effective as T
from third_order_constants import tree_constants, marked
from second_order_constants import dressed
from junction_constants import counts
from constants_upper import beta0

d = 3
PB.load_terms(d, os.path.join(up, 'orders', 'd3_s9.json'), {10: 10, 11: 11, 12: 24, 13: 0, 14: 70}, 15)
KMAX, NMAX = 60, 40

def basic(beta, g):
    """the beta-dependent primitives, as in lower3e (THETA7 with CLOSED_EXACT, CLOSED_FROM)"""
    K = 4*(2*d-3); eps = expm1(2*beta); x = e**2*K*eps
    th = 2*(d-2)*e**6*dressed(beta, g)[1] + x**T.CLOSED_FROM/(1-x) + fsum(c*(e*eps)**n for n, c in T.CLOSED_EXACT[d].items())
    lam = exp(2*(K+1)*th)
    u = g['u']; a = lam*u; gam = expm1(2*beta) + beta**2*exp(beta); b = lam*gam
    Ubar, s6 = dressed(beta, g)
    b2 = lam**2*(2*s6 + s6**2)/u**4
    return dict(K=K, th=th, lam=lam, u=u, a=a, b=b, gam=gam, Ubar=Ubar, s6=s6, b2=b2, Ups=g['Ups'])

def capture(beta, g, varpi_pow):
    """lower3e(beta, d, g, varpi_pow, chars=True) with every parameter returned (code mirrored from third_order_effective.lower3e)"""
    q = basic(beta, g); a, b, gam, Ubar, s6, b2 = q['a'], q['b'], q['gam'], q['Ubar'], q['s6'], q['b2']
    varpi = b**varpi_pow; C = 1 + max(1/varpi, (1+gam)/Ubar); omega = C*a**4*(1+gam)
    Pv, KS, Tv, N = counts(d, 'infinite')
    for k, v in T.NOVR[d].items(): N[k] = mpf(v)
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
    P = dict(varpi=varpi, C=C, xF=xF, xC=xC, cT=cT, hJ=hJ, beta_hi=beta, varpi_pow=varpi_pow)
    # the search for delta', F, F6 and the majorant: delegate to the authors' lower3e, then re-derive the majorant choice
    lo = T.lower3e(beta, d, g, varpi_pow=varpi_pow, chars=True)
    P.update(F=lo['F'], F6=lo['F6'], dl=lo['delta'], bars=lo['bars'])
    assert lo['a'] == a and lo['b'] == b, 'mirror of lower3e differs'
    return P, lo

def evaluate(beta, g, P, freeze_g=None):
    """all conditions at beta with the parameters P frozen; returns dict name -> (lhs, rhs)"""
    q = basic(beta, g); a, b, gam, Ubar, s6, b2 = q['a'], q['b'], q['gam'], q['Ubar'], q['s6'], q['b2']
    varpi, xF, xC, cT, hJ, F, F6, dl = (P[k] for k in ('varpi', 'xF', 'xC', 'cT', 'hJ', 'F', 'F6', 'dl'))
    C = 1 + max(1/varpi, (1+gam)/Ubar); omega = C*a**4*(1+gam)
    if cT is None: cT = exp(hJ)/(1-omega)              # c_J recomputed as in Theorem mg:thm:lower2 (the scaling of Remark (vii))
    Pv, KS, Tv, N = counts(d, 'infinite')
    for k, v in T.NOVR[d].items(): N[k] = mpf(v)
    def Yf(z, p=0): return fsum(m**p*binomial(KS*m, m-1)*(b*exp(z))**m for m in range(1, NMAX))
    Y = [Yf(xC, j) for j in range(3)]
    A = [exp(cT)*fsum(k**j*N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX)) for j in range(3)]
    out = {}
    out['C_frozen'] = (C, P.get('C'))                  # C_varpi recomputed vs frozen (the max is 1/varpi in both)
    out['dec_xiF'] = (Pv*Y[0] + A[0], xF)
    out['dec_xiC'] = (4*Pv*Y[0] + 8*A[0], xC)
    out['dec_hJ'] = (max(C*k*Pv*Y[0]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX)), hJ)
    out['dec_cJ'] = (exp(hJ)/(1-omega), cT)
    def f(k):
        if k == 4: return mpf(1)
        if k == 6: return 36*F6*exp(6*xF)*a**2
        return F*k**2*exp(xF*k)*a**(k-4)
    mk = marked(xF, Y, A, f, Pv, N, KMAX)
    gF, gC, g2F, g2C, kF, kC = mk if freeze_g is None else freeze_g
    # the six inequalities of Definition mg:def:jconst with the constants used (frozen or minimal)
    B0 = fsum(N[k]*f(k)*exp(xF*k) for k in range(4, KMAX)); B1 = fsum(k*N[k]*f(k)*exp(xF*k) for k in range(4, KMAX))
    out['six_1'] = (Pv*gC*Y[1] + B0 + gF*A[1], gF)
    out['six_2'] = (4*Pv*gC*Y[1] + 8*B0 + 8*gF*A[1], gC)
    out['six_3'] = (Pv*g2C*Y[1] + A[0] + g2F*A[1], g2F)
    out['six_4'] = (4*Pv*g2C*Y[1] + 8*A[0] + 8*g2F*A[1], g2C)
    out['six_5'] = (gF*g2F + Pv*kC*Y[2] + g2F*B1 + gF*A[1] + kF*A[2], kF)
    out['six_6'] = (gC*g2C + 4*Pv*kC*Y[2] + 8*(g2F*B1 + gF*A[1] + kF*A[2]), kC)
    def X(n, k):
        if k == 5: return 0
        return min(binomial(4*n+4, k), (4*n+4)*N[k])
    Ak = lambda kmin: 8*exp(cT)*fsum(N[k]*a**(mpf(k)/2) for k in range(kmin, KMAX))
    Fk = lambda kmin: 8*fsum(N[k]*f(k) for k in range(kmin, KMAX))
    # the majorant with the frozen (bbar, abar) and kind; the barred constants do not depend on beta (cT, F, F6, xF frozen)
    bb, ab, kind = P['bars']
    assert kind == 'corners', P['bars']
    assert b <= bb and a <= ab
    xFb, xCb, Yb, Ab = tree_constants(bb, ab, cT, Pv, KS, N, KMAX, NMAX)
    def fb(k):
        if k == 4: return mpf(1)
        if k == 6: return 36*F6*exp(6*xF)*ab**2
        return F*k**2*exp(xF*k)*ab**(k-4)
    gFb, gCb, g2Fb, g2Cb, kFb, kCb = marked(xFb, Yb, Ab, fb, Pv, N, KMAX)
    S3 = 4*Pv*exp(4*xFb)*(kCb*Yb[2] + 4*(gCb*g2Fb + g2Cb*gFb)*Yb[1] + 16*kFb*Yb[0])
    exact = fsum(c*b**n*a**E*exp(m*cT)*(f(kp)/a**(kp-4) if kp > 4 else 1) for (c, n, E, kp, m) in T.EXACT[d])
    corner = max((b/bb)**cn*(a/ab)**ce for (cn, ce) in T.CORNERS[d])*S3
    multi = exact + corner
    full = b + (4*(Pv-1)-KS)*b**2 + 48*binomial(18, 2)*(d-2)*(2*d-3)*b**3 + fsum(2*binomial(4*Pv+Tv*n, n)*(n+1)*(2*d-3)*b**n for n in range(4, NMAX)) \
           + fsum(2*binomial(4*Pv+Tv*n, n)*b**n*fsum(X(n, k)*f(k) for k in range(5, min(4*n+5, 199))) for n in range(1, NMAX)) \
           + 4*Pv*exp(4*xF)*(kC*Y[2] + 4*(gC*g2F + g2C*gF)*Y[1] + 16*kF*Y[0])
    eta = C*a**4*(2*Ubar*(1+dl) + full)/(1-omega)
    ZS = T.ZSMALL[d]
    parts = dict(two=(4*(Pv-1)-KS)*varpi*min(b**2, b2),
                 four=(10*a**4 + 2*eta*b**4)*(d-3) + fsum(ZS[n]*b**n for n in ZS)
                      + fsum(2*binomial(4*Pv+Tv*n, n)*(n+1)*(2*d-3)*b**n for n in range(max(ZS) + 1, NMAX)),
                 exits=KS*(a*f(6) + b*f(7)) + (4*(Pv-1)-2*KS)*b*f(8) + b**2*fsum(T.ZEXIT[d][k]*f(k) for k in T.ZEXIT[d])
                       + fsum((b**n*fsum(c*f(k) for k, c in T.ZEXITN[d][n].items()) if n in T.ZEXITN[d]
                               else 2*binomial(4*Pv+Tv*n, n)*b**n*fsum(X(n, k)*f(k) for k in range(5, min(4*n+5, 199)))) for n in range(3, NMAX)),
                 one_plaq=4*Pv*b*(Fk(8)*Ak(4) + Fk(4)*Ak(8) + Fk(4)*T.em1x(Ak(4))),
                 multi=multi)
    S = fsum(parts.values())
    X6 = lambda n, k: 0 if k == 5 else min(binomial(4*n+6, k), (4*n+6)*N[k])
    sig6 = 2*a + fsum(2*(2*d-3)*(mpf(2*n+3)/2)*binomial(6*Pv+Tv*n, n)*b**n for n in range(2, NMAX)) \
           + fsum(2*binomial(6*Pv+Tv*n, n)*b**n*fsum(X6(n, k)*f(k) for k in range(5, min(4*n+7, 199))) for n in range(1, NMAX)) \
           + 6*Pv*exp(6*xF)*(kC*Y[2] + 6*(gC*g2F + g2C*gF)*Y[1] + 36*kF*Y[0])
    need = Pv*(gC*Y[1]/7 + gF*Y[0])
    out['sigma'] = (S, Ubar*dl)
    out['exits'] = (need, F*(Ubar*(1+dl)-a))
    out['six'] = (sig6, 36*F6*exp(6*xF)*(Ubar*(1+dl)-a))
    out['old'] = (omega + C*a**4*(2*Ubar*(1+dl) + full)/(1-omega), mpf(1))
    out['eta<=1'] = (eta, mpf(1))
    out['beyond'] = (a**4*Ubar*(1+dl), varpi)
    out['growth'] = (a**4*Ubar*(1+dl), mpf(1))
    out['omega'] = (omega, mpf(1))
    out['barred_b'] = (b, bb); out['barred_a'] = (a, ab)
    Dfr = 8*(q['K']+1)*q['th'] + log(Ubar/q['Ups']) + log(1+dl)
    return out, dict(parts=parts, S=S, D=Dfr, q=q, marked=(gF, gC, g2F, g2C, kF, kC), minimal=mk)

def report(beta, out, extra, Dhi):
    worst = min(((r - l)/abs(r) if r != 0 else mpf(0), k) for k, (l, r) in out.items() if k != 'C_frozen')
    bad = [k for k, (l, r) in out.items() if k != 'C_frozen' and not l <= r]
    print(f"beta = {float(beta):.7g}: all conditions hold: {not bad} {bad if bad else ''}; smallest relative margin {float(worst[0]):.3e} ({worst[1]});"
          f"  D_frozen/beta^4 = {float(extra['D']/beta**4):.5g}, D_frozen <= D(beta_hi): {extra['D'] <= Dhi};  S/beta^4 = {float(extra['S']/beta**4):.5g}")
    print('    margins: ' + ', '.join(f"{k} {float((r-l)/abs(r)) if r else 0:.3g}" for k, (l, r) in out.items() if k != 'C_frozen'))
    sys.stdout.flush()

if __name__ == "__main__":
    mp.dps = 60
    args = sys.argv[1:]
    bhi = beta0(d) if args[0] == 'b0' else mpf(args[0])
    g = dict(T.GROUPS)['SU(3)'](bhi)
    t0 = time.time()
    best = None
    for vp in (2, mpf(5)/2):                         # as practical_beta.deficit: the smaller delta' of varpi = b^2, b^(5/2)
        try: P, lo = capture(bhi, g, vp)
        except AssertionError as ex: print('  varpi_pow', vp, 'fails:', ex); continue
        print(f"  varpi = b^{float(vp):g}: delta'/beta^4 = {float(lo['delta']/bhi**4):.6g}, bars {lo['bars']}")
        if best is None or lo['delta'] < best[1]['delta']: best = (P, lo)
    P, lo = best
    Dhi = lo['loss'] - log(g['Ups'])
    print(f"beta_hi = {float(bhi):.10g}: D/beta^4 = {float(Dhi/bhi**4):.6g} (authors' lower3e), varpi = b^{float(P['varpi_pow']):g}, (bbar, abar) = {P['bars']},"
          f" F = {float(P['F']):.6g}, F6 = {float(P['F6']):.6g}, xi_F = {float(P['xF']):.6g}, xi_C = {float(P['xC']):.6g}, h_J = {float(P['hJ']):.6g},"
          f" c_J = {float(P['cT']):.8g}  ({time.time()-t0:.0f} s)")
    # at beta_hi itself: the mirror must reproduce lower3e's hypotheses, and its constants are then frozen
    out, ex = evaluate(bhi, g, P)
    print(f"  mirror at beta_hi: S = {float(ex['S']):.12g} vs lower3e {float(fsum(lo['parts'].values())):.12g}; D {float(ex['D']):.12g} vs {float(Dhi):.12g}")
    report(bhi, out, ex, Dhi)
    frozen_g = ex['marked']
    for s in args[1:]:
        beta = mpf(s); gb = dict(T.GROUPS)['SU(3)'](beta)
        o1, e1 = evaluate(beta, gb, P)                       # minimal marked constants at beta
        report(beta, o1, e1, Dhi)
        o2, e2 = evaluate(beta, gb, P, freeze_g=frozen_g)    # the marked constants of beta_hi, frozen
        print('    with the marked constants of beta_hi frozen: all conditions hold:', all(l <= r for k, (l, r) in o2.items() if k != 'C_frozen'),
              [k for k, (l, r) in o2.items() if k != 'C_frozen' and not l <= r])
        print(f"    ratios to beta_hi: S {float(e1['S']/ex['S']):.5f}, beta/beta_hi {float(beta/bhi):.5f}, Ubar {float(e1['q']['Ubar']/ex['q']['Ubar']):.8f},"
              f" a/beta {float((e1['q']['a']/beta)/(ex['q']['a']/bhi)):.8f}, b/beta {float((e1['q']['b']/beta)/(ex['q']['b']/bhi)):.8f},"
              f" b2/beta^2 {float((e1['q']['b2']/beta**2)/(ex['q']['b2']/bhi**2)):.8f}")
        sys.stdout.flush()

# Referee check (Part IV, Lemma mg:lem:uppertails, Remark mg:rem:upperexact (ii)-(iv)): the tails T'_S, the error of Theorem mg:thm:upper2 with them,
# the bound E(beta) of Remark (ii), its monotonicity, and the thresholds in d = 3, 4.  Written for this report from the formulas of the paper.
#
# SU(3) constants by Weyl integration on the maximal torus (own implementation; compared with the authors' second_order_constants.su3):
#   c_r = int conj(chi_r) e^{(beta/3) Re chi} dU,  u = c_3/(3 c_0),  Upsilon_+ = 2 Var(Re chi), Upsilon_- = 2 Var(Im chi) under w/c_0,
#   kappa = 2 int (Re chi)^2 f_beta.  sigma_5 <= 18 u^5 + u_2^4 (e^beta - 1), u_2 <= e^beta - 1 - beta (Lemma mg:lem:uppertails, mg:lem:capelement(iv)).
# Usage: python3 upper_eval.py
import sys, os
from mpmath import mp, mpf, e, exp, log, sqrt, expm1, pi, cos, sin, fsum
mp.dps = 50

def su3_consts(beta, M=40):
    beta = mpf(beta)
    acc = dict(c0=0, cr=0, re2=0, im2=0, re=0, im=0)
    S = []
    for i in range(M):
        for j in range(M):
            t1, t2 = 2*pi*i/M, 2*pi*j/M
            t3 = -t1 - t2
            X_re = cos(t1) + cos(t2) + cos(t3); X_im = sin(t1) + sin(t2) + sin(t3)
            # |Delta|^2 = prod_{j<k} |e^{i t_j} - e^{i t_k}|^2 = prod 4 sin^2((t_j - t_k)/2)
            vd = 64*sin((t1 - t2)/2)**2*sin((t1 - t3)/2)**2*sin((t2 - t3)/2)**2/6
            w = vd*exp(beta/3*X_re)
            S.append((vd, w, X_re, X_im))
    Z = fsum(s[0] for s in S)                       # normalisation of the Haar measure (should be M^2)
    c0 = fsum(s[1] for s in S)/Z
    m_re = fsum(s[1]*s[2] for s in S)/Z/c0          # <(1+f) Re chi> = N u
    m_im = fsum(s[1]*s[3] for s in S)/Z/c0
    re2 = fsum(s[1]*s[2]**2 for s in S)/Z/c0
    im2 = fsum(s[1]*s[3]**2 for s in S)/Z/c0
    re2_haar = fsum(s[0]*s[2]**2 for s in S)/Z       # = 1/2
    u = m_re/3
    Up = 2*(re2 - m_re**2); Um = 2*(im2 - m_im**2)
    kappa = 2*(re2 - re2_haar)                       # (1/tau0) int (Re chi)^2 f
    return dict(u=u, Up=Up, Um=Um, Ups=max(Up, Um), kappa=kappa, c0=c0, Z=Z, M=M)

def beta0(d):
    K = 4*(2*d - 3)
    return log(1 + 1/(e**2*K*(K + 2)))/2

# counts: (N^cap, {n: N_n non-cap}) ; d=3 to n=10 (Remark (i), reproduced) and to 12 (this check); d=4 to 8 (Remark (iv)) and to 10 (this check)
COUNTS = {
    (3, 10): {'Q': (6, {5: 0, 6: 12, 7: 0, 8: 0, 9: 30, 10: 93}), 'A': (2, {5: 0, 6: 8, 7: 0, 8: 0, 9: 10, 10: 47})},
    (3, 12): {'Q': (6, {5: 0, 6: 12, 7: 0, 8: 0, 9: 30, 10: 93, 11: 144, 12: 162}), 'A': (2, {5: 0, 6: 8, 7: 0, 8: 0, 9: 10, 10: 47, 11: 60, 12: 96})},
    (4, 8): {'Q': (18, {5: 0, 6: 60, 7: 0, 8: 24}), 'A': (4, {5: 0, 6: 32, 7: 0, 8: 0})},
    (4, 9): {'Q': (18, {5: 0, 6: 60, 7: 0, 8: 24, 9: 262}), 'A': (4, {5: 0, 6: 32, 7: 0, 8: 0, 9: 60})},
    (4, 10): {'Q': (18, {5: 0, 6: 60, 7: 0, 8: 24, 9: 262, 10: 1143}), 'A': (4, {5: 0, 6: 32, 7: 0, 8: 0, 9: 60, 10: 506})},
}

def tails(beta, d, g, counts, theta_kind='theta5'):
    beta = mpf(beta)
    K = 4*(2*d - 3); eps = expm1(2*beta); x = e**2*K*eps
    th = x**5/(1 - x) if theta_kind == 'theta5' else x/(1 - x)
    w = eps*exp((K + 1)*th); q = e*K*w
    u = g['u']; u2 = expm1(beta) - beta
    s5 = 18*u**5 + u2**4*expm1(beta)
    out = {}
    for S, nC, PS, Jp in (('Q', 6, 24*d - 42, 24*(d - 2)), ('A', 1, 8*d - 11, K)):
        ncap, Nn = counts[S]; n0 = max(Nn)
        parts = dict(cap=ncap*s5*exp(5*(K + 1)*th), sets=fsum(c*w**n for n, c in Nn.items()), tail=(1 + mpf(1)/K)**Jp*q**(n0 + 1)/(1 - q))
        pre = (1 + eps)**nC*exp(th*PS)
        out[S] = (pre*sum(parts.values()), {k: pre*v for k, v in parts.items()})
    return out, dict(eps=eps, x=x, th=th, w=w, q=q, s5=s5, K=K)

def errors(beta, d, g, counts):
    """the exact error D'_2 + log Ups of Theorem mg:thm:upper2 with T'_Q, T'_p, and the bound E of Remark (ii)"""
    beta = mpf(beta); N, iota = 3, 0; tau0 = mpf(1)/2
    out, c = tails(beta, d, g, counts)
    TQ, Tp = out['Q'][0], out['A'][0]; eps = c['eps']
    u, Ups = g['u'], g['Ups']
    psi1 = sqrt(2)*3*u if g['Up'] >= g['Um'] else mpf(0)      # a_1 = |<f psi>|
    a2 = (1 + Tp)*psi1 + Tp
    V2 = (1 + Tp)*(Ups + psi1**2) + Tp - max((1 - Tp)*psi1 - Tp, 0)**2
    L2 = (1 - TQ)/(1 + eps**6)*max(Ups - Tp*(1 + psi1)*sqrt(2 - iota)*N*u, 0)**2 - (1 + a2**2)*exp(16*beta)*(mpf(N)**4/tau0**2)*TQ
    exact = log(V2/L2) + log(Ups)
    a = sqrt(2 - iota)*N*u; ap = (1 + Tp)*a + Tp
    B = TQ + eps**6 + 2*a*(1 + a)*Tp + (1 + ap**2)*exp(16*beta)*(mpf(N)**4/tau0**2)*TQ
    E = Tp*(2 + 3*a + 3*a**2) + B/(1 - B)
    return dict(exact=exact, E=E, TQ=TQ, Tp=Tp, B=B, parts=out, c=c)

def threshold(d, counts, lo, hi, which='exact'):
    """largest beta in [lo, hi] (bisection, assuming one crossing) with error < beta^3/432"""
    lo, hi = mpf(lo), mpf(hi)
    f = lambda b: errors(b, d, su3_consts(b), counts)[which] - b**3/432
    if f(hi) < 0: return hi, True
    for _ in range(40):
        mid = (lo + hi)/2
        if f(mid) < 0: lo = mid
        else: hi = mid
    return lo, False

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(here, '..', '..', 'second-order')); sys.path.insert(0, os.path.join(here, '..', '..', 'upper-bound'))
    b03 = beta0(3); b04 = beta0(4)
    print(f"beta_0(3) = {float(b03):.7e}, beta_0(4) = {float(b04):.7e}; beta^3/432 at beta_0(3) = {float(1/(432*b03)):.5f} beta^4")
    # 1. constants against the authors' su3
    try:
        from second_order_constants import su3 as su3_auth
        for b in (b03, mpf('1e-4'), mpf('1e-6')):
            mine = su3_consts(b); th = su3_auth(b)
            print(f"  beta={float(b):.3e}: u/beta = {float(mine['u']/b):.12f} (authors {float(th['u']/b):.12f}); Ups+ - 1 = {float(mine['Up'] - 1):.6e} "
                  f"(authors {float(th['Upm'][0] - 1):.6e}); Ups- - 1 = {float(mine['Um'] - 1):.6e}; kappa = {float(mine['kappa']):.6e}; kappa - 18u^2 = {float(mine['kappa'] - 18*mine['u']**2):.4e}")
    except Exception as ex:
        print("  (authors' su3 not available:", ex, ")")
    # 2. d = 3 at beta_0(3), n0 = 10 and 12
    for key in ((3, 10), (3, 12)):
        g = su3_consts(b03); r = errors(b03, 3, g, COUNTS[key])
        pq = {k: float(v/b03**4) for k, v in r['parts']['Q'][1].items()}; pa = {k: float(v/b03**4) for k, v in r['parts']['A'][1].items()}
        print(f"d=3, n0={key[1]}, beta_0(3): T'_Q = {float(r['TQ']/b03**4):.4e} beta^4 {pq}; T'_A = {float(r['Tp']/b03**4):.4e} beta^4 {pa}")
        print(f"   q = {float(r['c']['q']):.5f}, sigma_5/u^5 = {float(r['c']['s5']/g['u']**5):.6f}, (eps/u) = {float(r['c']['eps']/g['u']):.4f}; "
              f"exact error = {float(r['exact']/b03**4):.5f} beta^4; E = {float(r['E']/b03**4):.5f} beta^4; B = {float(r['B']/b03**4):.5f} beta^4")
    # 3. monotonicity of E/beta^5, u/beta, sigma_5/beta^5 on [1e-12, beta_0(3)]
    pts = sorted([b03*mpf(10)**(-mpf(k)/10) for k in range(0, 91)])
    Ev, uv, sv, okU = [], [], [], True
    for b in pts:
        g = su3_consts(b); r = errors(b, 3, g, COUNTS[(3, 10)])
        Ev.append(r['E']/b**5); uv.append(g['u']/b); sv.append(r['c']['s5']/b**5)
        okU = okU and g['Ups'] >= 1 and g['kappa'] - 18*g['u']**2 > 0
        assert r['exact'] <= r['E'] + mpf(10)**-40
    inc = lambda v: all(x < y for x, y in zip(v, v[1:]))
    print(f"monotone on {len(pts)} points of [{float(pts[0]):.1e}, beta_0(3)]: E/beta^5 {inc(Ev)}, u/beta {inc(uv)}, sigma5/beta^5 {inc(sv)}; Ups >= 1 and kappa > 18u^2 everywhere: {okU};"
          f" exact error <= E everywhere: True")
    print(f"   E(beta)/beta^4 at beta_0(3) = {float(Ev[-1]*b03):.5f};  sup over the grid of error/beta^4 = {float(max(v*b for v, b in zip(Ev, pts))):.5f}")
    # 4. thresholds
    for key in ((4, 8), (4, 9), (4, 10)):
        t, whole = threshold(4, COUNTS[key], mpf('1e-6'), b04)
        g = su3_consts(mpf('1.2e-5')); r = errors(mpf('1.2e-5'), 4, g, COUNTS[key])
        g0 = su3_consts(b04); r0 = errors(b04, 4, g0, COUNTS[key])
        print(f"d=4, n0={key[1]}: error < beta^3/432 for beta <= {float(t):.4e}{' (the whole range)' if whole else ''} = beta_0(4)/{float(b04/t):.3f}; "
              f"at 1.2e-5: error = {float(r['exact']/mpf('1.2e-5')**4):.3e} beta^4; at beta_0(4): error = {float(r0['exact']/b04**4):.4g} beta^4 "
              f"against beta^3/432 = {float(1/(432*b04)):.4g} beta^4 (E = {float(r0['E']/b04**4):.4g})")
    # 5. d = 4 with the counts to 9 and 10 plaquettes of this check: E/beta^5 monotone, and E(beta_0(4)) against beta^3/432
    for key in ((4, 9), (4, 10)):
        pts4 = sorted([b04*mpf(10)**(-mpf(k)/10) for k in range(0, 91)])
        vals = []
        for b in pts4:
            r = errors(b, 4, su3_consts(b), COUNTS[key]); vals.append(r['E']/b**5)
            assert r['exact'] <= r['E'] + mpf(10)**-40
        E0 = vals[-1]*b04**5
        print(f"d=4, n0={key[1]}: E/beta^5 increasing on {len(pts4)} points of [{float(pts4[0]):.1e}, beta_0(4)]: {inc(vals)}; E(beta_0(4)) = {float(E0/b04**4):.4g} beta_0^4, "
              f"so the error is <= {float(E0/b04**4):.4g} beta^4 < beta^3/432 (>= {float(1/(432*b04)):.4g} beta^4) on the whole range: {E0/b04**4 < 1/(432*b04)}")

# Section 13 of Part IV: the constants of Theorem mg:thm:lower4, i.e. Theorem mg:thm:lower3 with the junction sum of Lemma mg:lem:jsums4, built on the
# junction lemmas of Lemma mg:lem:betafour (A = 13.1(i), D = 13.1(ii), B = 13.1(iii), C = 13.1(iv)):
#   (A) Z_{4,4} = 12(d-3): the junctions with four plaquettes and two four-edge front-pairs are the turns, jumps and U-turns through unit cubes;
#   (B) one plaquette, exactly three front-pairs, one of four edges => another has >= 8 edges;
#   (C) two plaquettes, exactly three front-pairs, two of four edges => the third has >= 6 edges;
#   (D) one plaquette, two front-pairs, |p| = 4, |p'| >= 5 => p' is bd(q) xor bd(c) (6 edges) or bd(q) u bd(c) (7 edges), c ~ q, in the other layer
#       (K_S choices each), or has 8 edges and c shares only a vertex with q (4(P_v-1) - 2K_S choices).
# The junction sum Sigma_4^(3e) of the revised Lemma mg:lem:jsums3:
#   shifts and jumps (unchanged); 12(d-3) b^4 + crude counts for n >= 5 plaquettes; exits K_S b (f(6)+f(7)) + (4(P_v-1)-2K_S) b f(8) + crude for n >= 2;
#   one plaquette 4 P_v b (F_8 A + F A_8 + F (e^A - 1 - A));
#   two or more plaquettes and >= 3 front-pairs: [(b/bbar)^3 (a/abar)^2 + (b/bbar)^2 (a/abar)^3] x (doubly marked tree sum at (bbar, abar)), since by (C)
#   such a junction has b^n a^E with n >= 2, E >= 2 and n + E >= 5 (E = total power of a in its weight); or, whichever is smaller, the majorant of
#   Section 11, (b/bbar)^2 x (doubly marked tree sum at bbar, with the actual a). The first is O(beta^5), the second O(beta^4) with a smaller constant.
# Usage: python3 third_order_effective.py
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
for sub in ('second-order', 'upper-bound', 'lower-bound'): sys.path.insert(0, os.path.join(here, '..', sub))
sys.path.insert(0, here)
from mpmath import mp, mpf, e, exp, log, log1p, expm1, binomial, fsum, factorial
from second_order_constants import su2, su3, u1, dressed, upper2
from junction_constants import counts
from third_order_constants import tree_constants, marked, lower3
mp.dps = 60                                   # e^A - 1 - A and the small parts need more than 40 digits below beta ~ 1e-12

def em1x(A):
    """e^A - 1 - A without cancellation"""
    return fsum(A**m/factorial(m) for m in range(2, 60)) if A < mpf(1)/2 else exp(A) - 1 - A

# Section 16 options (off by default, so that the output above is that of Sections 13-14); set by practical_beta.py:
#   THETA7: the ratio bound with theta_7 (Lemma mg:lem:theta7) in place of theta_5;  NOVR[d]: exact N^v_k;  ZSMALL[d][n]: exact Z_{n,4}, n = 5, 6, 7;
#   LOW[d] = (N_3, N_2, N'_2): the exact junctions with three front-pairs of lowest order and the majorant M_5 (Lemmas mg:lem:threelow, mg:lem:jsums5).
THETA7 = False; ZSMALL = None; NOVR = None; LOW = None
CLOSED_EXACT = None   # {d: {n: number of closed polymers of n plaquettes through a plaquette}}, for 10 <= n < CLOSED_FROM
CLOSED_FROM = 8   # the smallest size of a closed polymer other than a cube used in theta (8: Lemma mg:lem:seven; 10: no closed sets of 8, 9 plaquettes)
ZEXITN = None  # {d: {n: {k: Z_(n,k)}}}: exact exits with n >= 3 plaquettes (zexit_count_n.py), in place of the crude counts
ZEXIT = None   # {d: {k: Z_(2,k)}}: exact two-plaquette exits with k >= 5 edges (zexit_count.py), in place of the crude counts
EXACT = None; CORNERS = None   # the general corner majorant (Section 16): see the M_5 block below
RMAX = False   # Lemma mg:lem:jsums5: the corner factors of M_4 and R as maxima instead of sums (each remaining pair dominates one corner)
NV_DECISION = False   # N^v_k beyond the exact values by the decision-tree bound (Lemma mg:lem:frontsdecision) instead of (d-1) binom((4d-6)k, k-1)
_DEC = {}
def decision_bound(D, k, grid=400):
    """min over s of ((1+s)^{2D} - 1 - 2Ds) (y/s)^k, y = max((1+s)^{2D-1} - 1, (1+s)^{D-1}): an upper bound for the number of fronts of k edges through
    a vertex of Z^D (the decision search with the potential y^{min(deg, 2)} per open vertex and weight s/y per edge)"""
    if (D, k) not in _DEC:
        best = None
        for i in range(1, grid):
            s_ = mpf(i)/grid
            y = max((1 + s_)**(2*D - 1) - 1, (1 + s_)**(D - 1))
            v = ((1 + s_)**(2*D) - 1 - 2*D*s_)*(y/s_)**k
            if best is None or v < best: best = v
        _DEC[(D, k)] = best
    return _DEC[(D, k)]
def lower3e(beta, d, g, volume='infinite', KMAX=60, NMAX=40, varpi_pow=2, bbar_grid=None, abar_grid=None, sharp6=True, uturn=True, seven=True, chars=False,
            fixed=None):
    """the smallest admissible delta' and the deficit; with fixed = dict(varpi, xF, xC, hJ, F, F6, dl), these parameters are used as given (no fixed-point
    iterations; c_J = e^{h_J}/(1-omega) as in Theorem mg:thm:lower2) and all hypotheses, including the three conditions of Lemma mg:lem:jsums on the
    decoration constants, are checked (Remark mg:rem:nextorder (vii): the scaling below a computed point)"""
    K = 4*(2*d-3); eps = expm1(2*beta); x = e**2*K*eps; th5 = x**5/(1-x)
    if THETA7:                                    # Lemma mg:lem:theta7: theta_7 = 2(d-2) e^6 sigma_6 + x^8/(1-x); with CLOSED_FROM = 10 (no closed sets of 8 or 9
        th5 = 2*(d-2)*e**6*(g['s6'] if g.get('bounds') else dressed(beta, g)[1]) + x**CLOSED_FROM/(1-x)   # plaquettes, by exhaustive search: Remark mg:rem:practical) theta_10
        if CLOSED_EXACT and d in CLOSED_EXACT:    # closed polymers of n plaquettes through a plaquette counted exactly (activity at most eps^n)
            th5 += fsum(c*(e*eps)**n for n, c in CLOSED_EXACT[d].items())
    lam = exp(2*(K+1)*th5)
    u = g['u']; a = lam*u; gam = expm1(2*beta) + beta**2*exp(beta); b = lam*gam
    if g.get('bounds'):          # bound mode (rigorous_bounds.bound_group): Ubar_lo <= Ubar <= Ubar_hi, u = ubar >= u, u_low <= u (Theorem mg:thm:thirdcap)
        Ubar_lo, Ubar_hi, s6, ulow = g['Ubar_low'], g['Ubar_up'], g['s6'], g['u_low']
    else:
        Ubar, s6 = dressed(beta, g); Ubar_lo = Ubar_hi = Ubar; ulow = u
    varpi = b**varpi_pow if fixed is None else fixed['varpi']
    C = 1 + max(1/varpi, (1+gam)/Ubar_lo); omega = C*a**4*(1+gam)
    b2 = lam**2*(2*s6 + s6**2)/ulow**4
    Pv, KS, Tv, N = counts(d, volume)
    if NOVR and d in NOVR:
        for k, v in NOVR[d].items(): N[k] = mpf(v)
    if NV_DECISION:
        for k in N:
            if k >= 6 and not (NOVR and d in NOVR and k in NOVR[d]): N[k] = min(N[k], decision_bound(d - 1, k))
    if volume == 'finite': N[4] = Pv; N[5] = 0
    def Yf(z, p=0): return fsum(m**p*binomial(KS*m, m-1)*(b*exp(z))**m for m in range(1, NMAX))
    cT = mpf(1); xF = mpf(0); xC = mpf(0)
    if omega >= 1: raise AssertionError('omega >= 1')
    deco_ok = True
    hJ = mpf(0)
    if fixed is not None:                          # the decoration constants as given, c_J = e^{h_J}/(1-omega): check xi_F >= P_v Y_0(xi_C) + A_0,
        xF, xC, hJ = fixed['xF'], fixed['xC'], fixed['hJ']           # xi_C >= 4 P_v Y_0 + 8 A_0 and h_J >= C_varpi k P_v Y_0 a^{k/2} e^{xi_F k} (k >= 4)
        cT = exp(hJ)/(1-omega)
        A0 = exp(cT)*fsum(N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX))
        deco_ok = (xF >= Pv*Yf(xC) + A0 and xC >= 4*Pv*Yf(xC) + 8*A0
                   and hJ >= C*max(a**(mpf(k)/2)*k*Pv*Yf(xC)*exp(xF*k) for k in range(4, 12)))
    for _ in range(300 if fixed is None else 0):
        if cT > 50 or xF > 50 or xC > 50:          # the decoration constants diverge (e.g. varpi too small for this beta): no admissible delta'
            raise AssertionError('decoration constants diverge')
        A0 = exp(cT)*fsum(N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX))
        xF_n = Pv*Yf(xC) + A0; xC_n = 4*Pv*Yf(xC) + 8*A0
        hJ = C*max(a**(mpf(k)/2)*k*Pv*Yf(xC)*exp(xF*k) for k in range(4, 12))
        cT_n = exp(hJ)/(1-omega)
        done = abs(xF_n-xF) + abs(xC_n-xC) + abs(cT_n-cT) < mpf(10)**-30
        xF, xC, cT = xF_n, xC_n, cT_n
        if done: break
    if fixed is None:                          # raise the decoration constants by a relative 1e-15: their conditions then hold strictly (the fixed point
        from third_order_constants import INFLATE   # is approached from below, and the series are truncated)
        xF, xC, hJ = INFLATE*xF, INFLATE*xC, INFLATE*hJ
        cT = exp(hJ)/(1-omega)
    A = [exp(cT)*fsum(k**j*N[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, KMAX)) for j in range(3)]
    Y = [Yf(xC, j) for j in range(3)]
    def X(n, k):
        if k == 5: return mpf(4*n+4)/5 if (N[5] > 0 and n >= 4) else 0
        return min(binomial(4*n+4, k), (4*n+4)*N[k])
    Ak = lambda kmin: 8*exp(cT)*fsum(N[k]*a**(mpf(k)/2) for k in range(kmin, KMAX))
    Ac = lambda kmin, kmax=KMAX: 16*exp(cT)*fsum(N[k]*(1 + mpf(k)/8)*a**(mpf(k)/2) for k in range(kmin, kmax))
    if bbar_grid is None: bbar_grid = [mpf(2)**(-j) for j in range(3, 18)]       # the grid of the referee of Section 13 (finer than the first one)
    if abar_grid is None: abar_grid = [mpf(2)**(-j) for j in range(2, 26)]
    tcache = {}
    for bb in bbar_grid:                      # the tree constants at (bbar, abar) do not depend on F: compute them once
        if bb < b: continue
        for ab in abar_grid:
            if ab < a: continue
            tc = tree_constants(bb, ab, cT, Pv, KS, N, KMAX, NMAX)
            if tc is not None: tcache[(bb, ab)] = tc
    tcache11 = {}
    for bb in [mpf(2)**(-j) for j in range(3, 18)]:   # Section 11's majorant: bbar for the plaquettes only, the actual a
        if bb < b: continue
        tc = tree_constants(bb, a, cT, Pv, KS, N, KMAX, NMAX)
        if tc is not None: tcache11[bb] = tc
    def solve(F, F6, dl):
        def f(k):                              # f(6) has its own constant F6 when sharp6 (the sum at a six-edge front, Sigma^(6))
            if k == 4: return mpf(1)
            if k == 6 and sharp6: return 36*F6*exp(6*xF)*a**2
            return F*k**2*exp(xF*k)*a**(k-4)
        gF, gC, g2F, g2C, kF, kC = marked(xF, Y, A, f, Pv, N, KMAX)
        Fk = lambda kmin: 8*fsum(N[k]*f(k) for k in range(kmin, KMAX))
        Fc = lambda kmin, kmax=KMAX: 16*fsum(N[k]*(1 + mpf(k)/8)*f(k) for k in range(kmin, kmax))
        best = None
        for (bb, ab), tc in tcache.items():
                xFb, xCb, Yb, Ab = tc
                def fb(k, ab=ab):
                    if k == 4: return mpf(1)
                    if k == 6 and sharp6: return 36*F6*exp(6*xF)*ab**2
                    return F*k**2*exp(xF*k)*ab**(k-4)
                mkb = marked(xFb, Yb, Ab, fb, Pv, N, KMAX)
                if mkb is None: continue
                gFb, gCb, g2Fb, g2Cb, kFb, kCb = mkb
                S3 = 4*Pv*exp(4*xFb)*(kCb*Yb[2] + 4*(gCb*g2Fb + g2Cb*gFb)*Yb[1] + 16*kFb*Yb[0])
                val = (max if RMAX else (lambda *t: fsum(t)))((b/bb)**3*(a/ab)**2, (b/bb)**2*(a/ab)**3)*S3
                if best is None or val < best[0]: best = (val, (bb, ab))
                if LOW is not None:              # M_5 (Lemma mg:lem:jsums5): (n,E) = (3,2), (2,3), (2,7/2) exactly, the majorant from (4,2), (3,3), (2,4)
                    N3, N2, N2p = LOW[d]
                    val5 = exp(cT)*(2*N3*b**3*a**2 + N2*b**2*a**3 + N2p*b**2*a**(mpf(7)/2)) \
                           + (max if RMAX else (lambda *t: fsum(t)))((b/bb)**4*(a/ab)**2, (b/bb)**3*(a/ab)**3, (b/bb)**2*(a/ab)**4)*S3
                    if val5 < best[0]: best = (val5, (bb, ab, 'M5'))
                if EXACT is not None:            # the general corner majorant: exact (or upper-bound) counts of the pairs (J,p') in a down-closed set of (n,E),
                    # each term (count, n, E, |p'|, m) with weight count b^n a^E e^{m c_J} (f(|p'|) a^{4-|p'|} if |p'| > 4); the rest from the corners of the complement
                    vx = fsum(c*b**n*a**E*exp(m*cT)*(f(kp)/a**(kp-4) if kp > 4 else 1) for (c, n, E, kp, m) in EXACT[d]) \
                         + max((b/bb)**cn*(a/ab)**ce for (cn, ce) in CORNERS[d])*S3
                    if vx < best[0]: best = (vx, (bb, ab, 'corners'))
        for bb, tc in tcache11.items():
                xFb, xCb, Yb, Ab = tc
                mkb = marked(xFb, Yb, Ab, f, Pv, N, KMAX)
                if mkb is None: continue
                gFb, gCb, g2Fb, g2Cb, kFb, kCb = mkb
                S3 = 4*Pv*exp(4*xFb)*(kCb*Yb[2] + 4*(gCb*g2Fb + g2Cb*gFb)*Yb[1] + 16*kFb*Yb[0])
                val = (b/bb)**2*S3
                if best is None or val < best[0]: best = (val, (bb, 'a'))
        full = b + (4*(Pv-1)-KS)*b**2 + 48*binomial(18, 2)*(d-2)*(2*d-3)*b**3 + fsum(2*binomial(4*Pv+Tv*n, n)*(n+1)*(2*d-3)*b**n for n in range(4, NMAX)) \
               + fsum(2*binomial(4*Pv+Tv*n, n)*b**n*fsum(X(n, k)*f(k) for k in range(5, min(4*n+5, 199))) for n in range(1, NMAX)) \
               + 4*Pv*exp(4*xF)*(kC*Y[2] + 4*(gC*g2F + g2C*gF)*Y[1] + 16*kF*Y[0])
        # a downward exit (the U-turns: p' in the layer of p) costs eta = C a^4 (2 Ubar (1+delta') + Sigma_4)/(1-omega)  (step (c) of Thm 10.8)
        eta = C*a**4*(2*Ubar_hi*(1+dl) + full)/(1-omega) if uturn else mpf(1)
        parts = dict(two=(4*(Pv-1)-KS)*varpi*min(b**2, b2),
                     # chars (Theorem 14.2): the turns and jumps, widenings and narrowings by characters (Lemma 14.1), a in place of b
                     four=(10*(a**4 if chars else b**4) + 2*eta*b**4)*(d-3) + (fsum(ZSMALL[d][n]*b**n for n in ZSMALL[d])
                          + fsum(2*binomial(4*Pv+Tv*n, n)*(n+1)*(2*d-3)*b**n for n in range(max(ZSMALL[d]) + 1, NMAX)) if ZSMALL
                          else fsum(2*binomial(4*Pv+Tv*n, n)*(n+1)*(2*d-3)*b**n for n in range(5, NMAX))),
                     exits=KS*((a if chars else b)*f(6) + b*f(7)) + (4*(Pv-1)-2*KS)*b*f(8)
                           + (b**2*fsum(ZEXIT[d][k]*f(k) for k in ZEXIT[d]) if ZEXIT else fsum(2*binomial(4*Pv+Tv*2, 2)*b**2*X(2, k)*f(k) for k in range(5, 13)))
                           + fsum((b**n*fsum(c*f(k) for k, c in ZEXITN[d][n].items()) if (ZEXITN and d in ZEXITN and n in ZEXITN[d])
                                   else 2*binomial(4*Pv+Tv*n, n)*b**n*fsum(X(n, k)*f(k) for k in range(5, min(4*n+5, 199)))) for n in range(3, NMAX)),
                     one_plaq=4*Pv*b*(Fk(8)*Ak(4) + Fk(4)*Ak(8) + Fk(4)*em1x(Ak(4))),
                     multi=best[0])
        # the junction sum at a six-edge front-pair (two narrowings with one plaquette; crude counts for the rest)
        X6 = lambda n, k: 0 if k == 5 else min(binomial(4*n+6, k), (4*n+6)*N[k])
        sig6 = 2*(a if chars else b) + fsum(2*(2*d-3)*(mpf(2*n+3)/2)*binomial(6*Pv+Tv*n, n)*b**n for n in range(2, NMAX)) \
               + fsum(2*binomial(6*Pv+Tv*n, n)*b**n*fsum(X6(n, k)*f(k) for k in range(5, min(4*n+7, 199))) for n in range(1, NMAX)) \
               + 6*Pv*exp(6*xF)*(kC*Y[2] + 6*(gC*g2F + g2C*gF)*Y[1] + 36*kF*Y[0])
        kmin = 7 if (seven and sharp6) else 5                   # with f(6) separate and no five-edge fronts, the F-condition is needed for k >= 7 only
        return fsum(parts.values()), parts, Pv*(gC*Y[1]/kmin + gF*Y[0]), full, best[1], sig6
    F = mpf(1); F6 = mpf(1); dl = mpf(0)
    for _ in range(100 if fixed is None else 0):
        S, parts, need, full, bb, sig6 = solve(F, F6, dl); dln = S/Ubar_lo*(1+mpf(10)**-12)
        Fn = need/(Ubar_lo*(1+dln)-a)*(1+mpf(10)**-9)
        F6n = sig6/(36*exp(6*xF)*(Ubar_lo*(1+dln)-a))*(1+mpf(10)**-9) if sharp6 else F6
        done = abs(Fn-F) < mpf(10)**-12*F and abs(F6n-F6) < mpf(10)**-12*F6 and abs(dln-dl) < mpf(10)**-12*dln
        F, F6, dl = Fn, F6n, dln
        if done: break
    dl = dl*(1+mpf(10)**-9)                                      # final delta'; eta is evaluated at this delta'
    if fixed is not None: F, F6, dl = fixed['F'], fixed['F6'], fixed['dl']
    S, parts, need, full, bb, sig6 = solve(F, F6, dl)
    hyp = dict(sigma=S < Ubar_lo*dl, exits=need < F*(Ubar_lo*(1+dl)-a), omega=omega < 1, F6=F6 > 0,
               six=(sig6 < 36*F6*exp(6*xF)*(Ubar_lo*(1+dl)-a)) if sharp6 else True,
               old=omega + C*a**4*(2*Ubar_hi*(1+dl) + full)/(1-omega) <= 1, beyond=a**4*Ubar_hi*(1+dl) <= varpi, growth=a**4*Ubar_hi*(1+dl) <= 1,
               decorations=deco_ok)
    assert all(hyp.values()), hyp
    out_D = (8*(K+1)*th5 + exp(2*beta)*g['Delta'] + log1p(dl)) if g.get('bounds') else None   # bound mode: D <= 8(K+1)theta + e^{2beta} Delta + log(1+delta')
    return dict(delta=dl, parts=parts, b=b, a=a, F=F, F6=F6, sig6=sig6, Ubar=Ubar_hi, bars=bb, loss=8*(K+1)*th5 + log(Ubar_hi) + log1p(dl), D_bound=out_D,
                varpi=varpi, xF=xF, xC=xC, hJ=hJ, cT=cT)

GROUPS = (('SU(2)', su2), ('SU(3)', su3), ('U(1)', u1))
NIOTA = {'SU(2)': (2, 1), 'SU(3)': (3, 0), 'U(1)': (1, 0)}

def delta4_limit(d, N, iota, sharp6=True, varpi_pow=mpf(5)/2, uturn=True):
    """delta'_4 = lim delta'/beta^4 (Theorem 13.4): the moves of the tube at the sup norm.
    Shifts and jumps 8(4P_v-4-K_S)(2-iota)N^2u_1^2 if varpi = b^2 (O(beta^5) if varpi = b^3); turns and jumps through cubes 160(d-3)
    (and the U-turns 32(d-3), unless they are charged eta, uturn=True); bulges and hops 8K_S u_1^2 (two narrowings, sharp6=True),
    or (1872/5) K_S P_v^2 u_1^2 with the crude F of the first version."""
    Pv = 2*(d-1)*(d-2); KS = 4*(2*d-5); u1 = mpf(1+iota)/(2*N**2)
    shifts = 8*(4*Pv-4-KS)*(2-iota)*N**2*u1**2 if varpi_pow <= 2 else mpf(0)      # varpi = b^p with p > 2: the shifts and jumps are o(beta^4)
    cube = (160 if uturn else 192)*(d-3)
    exits = 8*KS*u1**2 if sharp6 else mpf(1872)/5*KS*Pv**2*u1**2
    parts = (shifts, mpf(cube), exits)
    return sum(parts), parts

def best_lower(beta, d, g):
    """the smallest admissible delta' of varpi = b^2, b^(5/2), b^3 (all are allowed by Theorem 13.4)"""
    out = None
    for vp in (2, mpf(5)/2, 3):
        try:
            lo = lower3e(beta, d, g, varpi_pow=vp); lo['varpi_pow'] = vp
        except AssertionError:
            continue
        if out is None or lo['delta'] < out['delta']: out = lo
    return out

def main():
    for gname, gf in GROUPS:
        for d in (3, 4):
            dl4, p4 = delta4_limit(d, *NIOTA[gname])
            dl42, _ = delta4_limit(d, *NIOTA[gname], varpi_pow=2)
            dl4c, _ = delta4_limit(d, *NIOTA[gname], sharp6=False, varpi_pow=2, uturn=False)
            print(f"{gname}, d={d}:  limit delta'_4 = {float(dl4):.6g} = {' + '.join(f'{float(x):.6g}' for x in p4)}  (varpi = b^(5/2) or b^3);"
                  f"  with varpi = b^2: {float(dl42):.6g};  first version (crude F, U-turns at b^4, varpi = b^2): {float(dl4c):.6g}")
            for beta in (mpf('1e-12'), mpf('1e-9'), mpf('1e-7'), mpf('1e-5'), mpf('1e-4')):
                g = gf(beta)
                lo = best_lower(beta, d, g)
                if lo is None:
                    print(f"  beta={float(beta):.0e}: hypotheses fail for varpi = b^2 and b^3"); continue
                lo_old = lower3(beta, d, g); up = upper2(beta, d, g)
                pc = ', '.join(f'{k}={float(v/beta**4):.3g}' for k, v in lo['parts'].items())
                print(f"  beta={float(beta):.0e}: delta'/beta^4 = {float(lo['delta']/beta**4):.4g} with varpi = b^{float(lo['varpi_pow']):g} (Section 11: {float(lo_old['delta']/beta**4):.3g});"
                      f"  window for m - 4log(1/u): [{float(-lo['loss']):+.7e}, {float(up['D2']):+.7e}], -log Ups = {float(-log(g['Ups'])):+.7e}")
                print(f"     parts/beta^4: {pc};  majorant at (bbar, abar) = {lo['bars'][0]}, {lo['bars'][1]};  Sigma^(6)/b = {float(lo['sig6']/lo['b']):.4g}")
                sys.stdout.flush()

if __name__ == "__main__":
    main()

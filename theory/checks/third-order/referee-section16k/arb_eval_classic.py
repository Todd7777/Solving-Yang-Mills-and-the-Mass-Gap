# Referee check (ELEVENTH pass on Section 16, Corollary mg:cor:classic), item 6: the lower bound of the corollary in BALL ARITHMETIC (Arb, python-flint).
# An adapted copy of the tenth referee's referee-section16j/arb_eval_cor.py (itself the ninth referee's referee-section16i/arb_eval.py). The method is theirs:
# the hypotheses of Theorem mg:thm:lower5 in bound mode (Lemma mg:lem:boundmode) at beta_* = beta_0(d)(1+1e-30), 400-bit balls, every series summed to KX
# fronts / NX plaquettes with its remainder bounded rigorously, the referee's own marked constants (six inequalities verified in Arb), and the authors'
# parameters (varpi, xi_F, xi_C, h_J, F, F_6, delta', bbar, abar) taken as exact binary numbers.
# Changes for the corollary (marked "16k" below):
#   - GROUP 2 = SU(2) with N = 2, iota = 1 (tau_0 = 1); GROUP 1 = U(1) with N = 1, iota = 0. The arb code before hard-coded iota = 0 and N = GROUP.
#   - u2bar = e^beta - 1 - beta (Lemma mg:lem:capelement(iv) / Corollary mg:cor:window), not Lemma mg:lem:u2sun.
#   - the authors' choice of varpi = b^p is reproduced (p = 2 or 5/2, the smaller delta' as in scaling_check.params; for U(1) the authors' search picks
#     p = 2), with their (bbar, abar) grid point, KMAX = 100 in d = 4 (third_order_effective.KMAX_DEFAULT; --kmax) and F_RAISE = 1e-7 for U(1) in d = 3 (--f-raise).
#   - the barred tree constants are taken from third_order_constants.tree_constants with the authors' KMAX of the case (100 in d = 4), not a fixed 60:
#     at abar = 2^-8 a truncation at 60 edges leaves them short by ~1e-13 relative, above their raise of 1e-15 (I check this below as well).
#   - KX = 300 fronts (not 200), as abar = 2^-8 makes the barred front series converge with ratio ~0.6 per edge; the remainder is bounded rigorously anyway.
#   - the comparisons are with the corollary's c_- and c'_- (1968, 570 for SU(2); 1.772e5, 8043 for U(1)).
#   - option --decision-d3 (a SUGGESTION, not the corollary's setting): in d = 3 the fronts beyond the exact counts (k > 12) are bounded by Lemma
#     mg:lem:frontsdecision for Z^2 (D = 2: N^v_k <= ((1+s)^4 - 1 - 4s)(y/s)^k, minimised over s on the 1/400 grid for k < 60, s = 0.325 beyond), as d = 4
#     already does for Z^3, instead of (d-1) binom((4d-6)k, k-1); the authors' pipeline is run with NV_DECISION = True and fronts to 100 edges, at the grid
#     point its full search picks (nv_decision_d3_probe.py): (bbar, abar) = (2^-7, 2^-6), varpi = b^(5/2), for both groups.
#   - addendum (the coordinator's follow-up): GROUP = 3, 4, 5 = SU(N) with iota = 0 and u2bar of Lemma mg:lem:u2sun (c_Lambda = 2 for N = 4), exactly as in
#     referee-section16j/arb_eval_cor.py, so that Corollary mg:cor:sectorscap in d = 3 can be certified with --decision-d3 too (and, without the option, the
#     tenth referee's certified values for d3_s9z.json are reproduced: a check of this extension).
# Usage: python3 arb_eval_classic.py GROUP d CMINUS [--decision-d3]   (GROUP = 1 (U(1)), 2 (SU(2)), 3, 4, 5 (SU(N)); d = 3 or 4; e.g. ... 1 3 177200)
import os, sys, math
here = os.path.dirname(os.path.abspath(__file__)); parent = os.path.join(here, '..')
sys.path.insert(0, parent)
from flint import arb, ctx
from mpmath import mp, mpf, nstr as mpnstr
import practical_beta as PB
import third_order_effective as T
import third_order_constants as TC

PREC = 400
KX, NX = 300, 110                                                    # 16k: fronts to 300 edges
ctx.prec = PREC
ONE = arb(1); E = ONE.exp()

B3 = '0.00040262143352281974401285396895067832618294529962141'   # beta_0(3)(1+1e-30), as in run_cap.sh
B4 = '0.00015376644804254827379172376907999414111262406766301'   # beta_0(4)(1+1e-30)
SETTINGS = {3: ('d3_s9z.json', {10: 10, 11: 11, 12: 24, 13: 0, 14: 70}, 15, B3),
            4: ('d4_s10x_92.json', {10: 60, 11: 66, 12: 320, 13: 0, 14: 1316, 15: 2520, 16: 12552}, 17, B4)}
# 16k: the authors' grid choices (orders/capz_su2_d3_scaling.out, cap_su2_d4_scaling.out, capz_u1_d3_scaling.out, cap_u1_d4_scaling.out)
GRID = {(2, 3): (2**-7, 2**-10), (2, 4): (2**-8, 2**-8), (1, 3): (2**-8, 2**-10), (1, 4): (2**-8, 2**-8),
        (3, 3): (2**-7, 2**-10), (4, 3): (2**-7, 2**-10), (5, 3): (2**-7, 2**-10),      # addendum: the tenth referee's grid points for d3_s9z.json (ABARS there),
        (3, 4): (2**-8, 2**-9), (4, 4): (2**-8, 2**-10), (5, 4): (2**-8, 2**-10)}       # and for d4_s10x_92.json
NIOTA = {2: (2, 1), 1: (1, 0), 3: (3, 0), 4: (4, 0), 5: (5, 0)}    # 16k: (N, iota)
DECISION_D3 = '--decision-d3' in sys.argv
if DECISION_D3:
    sys.argv.remove('--decision-d3')
    GRID[(2, 3)] = GRID[(1, 3)] = GRID[(3, 3)] = GRID[(4, 3)] = GRID[(5, 3)] = (2**-7, 2**-6)   # the point the full search picks (orders/capz_*_d3_scaling.out)

def A(x):
    """an exact arb from an mpmath number, an int or a decimal string"""
    if isinstance(x, int): return arb(x)
    if isinstance(x, str): return arb(arb(x).mid())
    return arb(arb(mpnstr(x, 70)).mid())

def comb(n, k): return arb(math.comb(n, k)) if 0 <= k <= n else arb(0)

def check(name, L, R, store):
    """L < R for the whole balls"""
    ok = bool((R - L) > 0)
    store[name] = (ok, (1 - L/R))
    return ok

def authors_params(G, d):
    """16k: the parameters at beta_* as the authors' scaling_check.params finds them (bound mode, certified counts), restricted to their grid point"""
    PB.GROUP = G; PB.BOUNDS = True; PB.CERTIFIED = True
    if d == 4: T.KMAX_DEFAULT = 100                                 # --kmax 100 (run_cap.sh classic)
    if G == 1 and d == 3: T.F_RAISE = mpf(10)**-7                 # --f-raise 1e-7
    terms, closed, cf, bstar = SETTINGS[d]
    PB.load_terms(d, os.path.join(parent, 'orders', terms), closed, cf)
    if DECISION_D3 and d == 3: T.NV_DECISION = True; T.KMAX_DEFAULT = 100       # load_terms sets NV_DECISION = False in d = 3
    beta = mpf(bstar); g = PB.group_data()(beta)
    bb, ab = GRID[(G, d)]
    best = None
    for vp in (2, mpf(5)/2):
        try: lo = T.lower3e(beta, d, g, varpi_pow=vp, chars=True, bbar_grid=[mpf(bb)], abar_grid=[mpf(ab)])
        except AssertionError: continue
        if best is None or lo['delta'] < best['delta']: best = lo; best['p'] = vp
    lo = best
    return dict(beta=beta, varpi=lo['varpi'], p=lo['p'], xF=lo['xF'], xC=lo['xC'], hJ=lo['hJ'], F=lo['F'], F6=lo['F6'], dl=lo['delta'], bb=lo['bars'][0],
                ab=lo['bars'][1], D_auth=lo['D_bound'], strong_auth=lo['strong'], g=g)

def run(G, d):
    mp.dps = 60
    P = authors_params(G, d)
    beta = A(P['beta'])
    K = 4*(2*d - 3)
    beta0 = (ONE + ONE/(E**2*K*(K + 2))).log()/2
    assert bool(beta - beta0 > 0), 'beta_* >= beta_0(d)'
    n_ = d - 1; Pv = 2*n_*(n_ - 1); KS = 4*(2*d - 5); Tv = 4*(Pv - 1)
    N, iota = NIOTA[G]                                             # 16k
    # --- group bounds (Lemmas mg:lem:char(iii), mg:lem:capelement(iv), mg:lem:dressed, mg:lem:upsbounds(i)) ---
    ubar = ((1 + iota)*beta + beta**2*beta.exp())/(2*N**2)
    ulow = (1 + iota)*beta*(-beta).exp()/(2*N**2)
    if N >= 3:                                                     # addendum: Lemma mg:lem:u2sun, as in referee-section16j/arb_eval_cor.py
        R = beta.expm1() - beta - beta**2/2
        dS, dA, dL = arb(N*(N + 1))/2, arb(N*N - 1), arb(N*(N - 1))/2
        u2 = ((beta**2/(8*N**2) + R/dS)/dS).max((beta**2/(4*N**2) + R/dA)/dA).max(R)
        if N >= 4: u2 = u2.max(((2 if N == 4 else 1)*beta**2/(8*N**2) + R/dL)/dL)
    else:
        u2 = beta.expm1() - beta                                   # 16k: Lemma mg:lem:capelement(iv)
    eps = (2*beta).expm1(); gam = eps + beta**2*beta.exp()
    s6 = (2 - iota)*N**2*ubar**6 + u2**5*beta.expm1()
    nB = (u2/ulow)**2*eps + arb(2 - iota).sqrt()*(N/ulow)*(u2**5*beta.expm1()).sqrt()
    c2 = (u2/ulow)**4*beta.exp()
    Ulo = (-2*beta).exp(); Uhi0 = beta.exp()
    assert bool(Ulo - c2 > 0) and bool(Ulo - ubar**4*(1 + gam) > 0), 'denominators of Lemma mg:lem:dressed'
    Delta = nB**2/(Ulo - c2) + 2*(2*beta).exp()*s6/(Ulo - ubar**4*(1 + gam))
    Uhi = Uhi0 + Delta
    # --- theta (Lemma mg:lem:theta7 with the closed polymers counted exactly), a, b ---
    closed, CF = SETTINGS[d][1], SETTINGS[d][2]
    x = E**2*K*eps
    th = 2*(d - 2)*E**6*s6 + sum((arb(c)*(E*eps)**n for n, c in closed.items()), arb(0)) + x**CF/(1 - x)
    lam = (2*(K + 1)*th).exp()
    a = lam*ubar; b = lam*gam
    varpi, xF, xC, hJ, F, F6, dl = (A(P[k]) for k in ('varpi', 'xF', 'xC', 'hJ', 'F', 'F6', 'dl'))
    bb, ab = A(P['bb']), A(P['ab'])
    # 16k: varpi must be b(beta_*)^p with the authors' p (Lemma mg:lem:scaling assumes varpi_* = bhat(beta_*)^p, 2 <= p <= 5/2); recheck it in Arb
    p = P['p']; pa = arb(2) if p == 2 else arb(5)/2
    varpi_check = bool(abs(varpi/b**pa - 1) < arb('1e-50'))
    C = 1 + (1/varpi).max((1 + gam)/Ulo)
    om = C*a**4*(1 + gam)
    assert bool(1 - om > 0)
    cJ = hJ.exp()/(1 - om)
    b2 = lam**2*(2*s6 + s6**2)/ulow**4
    # --- N^v_k: exact where counted; else (d-1) binom((4d-6)k, k-1), in d=4 also the decision bound (Lemma mg:lem:frontsdecision) on the 1/400 grid ---
    D_ = d - 1
    exact = dict(T.NOVR[d])
    def dec(s, k):
        s = arb(s)/400; y = ((1 + s)**(2*D_ - 1) - 1).max((1 + s)**(D_ - 1))
        return ((1 + s)**(2*D_) - 1 - 2*D_*s)*(y/s)**k
    s_fix = {2: 130, 3: 78}[D_]            # s = 0.325 (Z^2), 0.195 (Z^3): the bound used for the remainder beyond KX
    Nv = {4: arb(Pv), 5: arb(0)}
    KK = max(4*NX + 8, KX)
    for k in range(6, KK + 1):
        crude = n_*comb((4*d - 6)*k, k - 1)
        if k in exact: Nv[k] = arb(exact[k])
        elif (d == 4 or DECISION_D3) and k < 60: Nv[k] = crude.min(min((dec(s, k) for s in range(1, 400)), key=lambda v: v.upper()))
        elif d == 4 or DECISION_D3 or k >= 60: Nv[k] = crude.min(dec(s_fix, k))
        else: Nv[k] = crude
    s0 = arb(s_fix)/400; y0 = ((1 + s0)**(2*D_ - 1) - 1).max((1 + s0)**(D_ - 1)); C0 = (1 + s0)**(2*D_) - 1 - 2*D_*s0; Y0r = y0/s0
    def front_tail(w, j, Kc=KX):
        """sum_{k >= Kc} k^j N^v_k w^k <= C0 Kc^j r^Kc / (1 - (1+1/Kc)^j r), r = (y/s) w (Lemma mg:lem:frontsdecision)"""
        r = Y0r*w; q = (1 + ONE/Kc)**j*r
        assert bool(1 - q > 0)
        return C0*arb(Kc)**j*r**Kc/(1 - q)
    def cluster_tail(y, p_, M=NX):
        """sum_{m >= M} m^p binom(K_S m, m-1) y^m <= e y M^p (e K_S y)^(M-1) / (1 - (1+1/M)^p e K_S y)"""
        r = E*KS*y; q = (1 + ONE/M)**p_*r
        assert bool(1 - q > 0)
        return E*y*arb(M)**p_*r**(M - 1)/(1 - q)
    def poly_tail(c, bq, qM, qratio, M=NX):
        """sum_{n >= M} binom(c + Tv n, n) b^n q_n, with q_n <= qM qratio^(n-M): binom(c+Tn, n) <= (e(T + c/n))^n"""
        r = E*(Tv + arb(c)/M)*bq; q = r*qratio
        assert bool(1 - q > 0)
        return r**M*qM/(1 - q)
    def Ys(bc, z, j): return sum((arb(m)**j*comb(KS*m, m - 1)*(bc*z.exp())**m for m in range(1, NX)), arb(0)) + cluster_tail(bc*z.exp(), j)
    def As(aa, xf, j): return cJ.exp()*(sum((arb(k)**j*Nv[k]*aa**(arb(k)/2)*(xf*k).exp() for k in range(4, KX)), arb(0)) + front_tail(aa.sqrt()*xf.exp(), j))
    Fm = F.max(F6)
    def fval(k, aa):
        if k == 4: return arb(1)
        if k == 5: return arb(0)
        if k == 6: return 36*F6*(6*xF).exp()*aa**2
        return F*k**2*(xF*k).exp()*aa**(k - 4)
    def Bs(aa, xfB, j):       # sum_k k^j N_k f(k) e^{xfB k}; f(k) <= Fm k^2 e^{xF k} aa^(k-4) for k >= 7
        s = sum((arb(k)**j*Nv[k]*fval(k, aa)*(xfB*k).exp() for k in range(4, KX)), arb(0))
        return s + Fm*aa**-4*front_tail(aa*(xF + xfB).exp(), j + 2)
    res = {}
    res['varpi = bhat(beta_*)^p, p = ' + ('2' if p == 2 else '5/2')] = (varpi_check, arb(0))
    # --- decoration conditions (the authors' xi_F, xi_C, h_J) ---
    Y0 = Ys(b, xC, 0); A0 = As(a, xF, 0)
    check('deco xi_F', Pv*Y0 + A0, xF, res); check('deco xi_C', 4*Pv*Y0 + 8*A0, xC, res)
    w4 = a.sqrt()*xF.exp()
    assert bool(1 - 2*w4 > 0)       # (k+1) w^(k+1) <= k w^k for k >= 4 when w <= 1/2: the maximum over k is at k = 4
    check('deco h_J (k=4, the maximum)', C*4*Pv*Y0*a**2*(4*xF).exp(), hJ, res)
    # --- barred: the authors' xi-bar (tree_constants at their c_J, with their KMAX of the case), checked here ---
    mp.dps = 60
    kmax_auth = T.KMAX_DEFAULT                                     # 16k (100 in d = 4)
    Nv_up = {k: mpf(mpnstr(Nv[k].upper(), 60)) for k in range(4, max(kmax_auth, 200))}
    tc = TC.tree_constants(P['bb'], P['ab'], mpf(mpnstr(cJ.mid(), 60)), Pv, KS, Nv_up, kmax_auth, 40)
    xFb, xCb = A(tc[0]), A(tc[1])
    Yb0 = Ys(bb, xCb, 0); Ab0 = As(ab, xFb, 0)
    check('barred xi_F', Pv*Yb0 + Ab0, xFb, res); check('barred xi_C', 4*Pv*Yb0 + 8*Ab0, xCb, res)
    if d == 4:                                                    # 16k: what a truncation at 60 edges would have given (the authors' reason for --kmax 100)
        tc60 = TC.tree_constants(P['bb'], P['ab'], mpf(mpnstr(cJ.mid(), 60)), Pv, KS, Nv_up, 60, 40)
        xFb60 = A(tc60[0])
        res['_(info) barred xi_F with fronts truncated at 60 edges: room'] = (True, 1 - (Pv*Ys(bb, A(tc60[1]), 0) + As(ab, xFb60, 0))/xFb60)
    # --- marked constants (the ninth referee's construction) ---
    IOT = arb('1e-30')
    def up(x_): return arb(x_.upper())
    def marked(Y, Av, B):
        Y0_, Y1, Y2 = Y; A0_, A1, A2 = Av; B0, B1 = B
        det = (1 - A1)*(1 - 4*Pv*Y1) - Pv*Y1*8*A1
        assert bool(det > 0) and bool(1 - A1 > 0) and bool(1 - 4*Pv*Y1 > 0)
        g_ = [(B0*(1 - 4*Pv*Y1) + Pv*Y1*8*B0)/det, ((1 - A1)*8*B0 + 8*A1*B0)/det,
              (A0_*(1 - 4*Pv*Y1) + Pv*Y1*8*A0_)/det, ((1 - A1)*8*A0_ + 8*A1*A0_)/det]
        gF, gC, g2F, g2C = [up(v*(1 + IOT)) for v in g_]
        dt = (1 - 4*Pv*Y2)*(1 - A2) - 8*A2*Pv*Y2
        assert bool(dt > 0) and bool(1 - A2 > 0) and bool(1 - 4*Pv*Y2 > 0)
        r1 = gC*g2C + 8*(g2F*B1 + gF*A1); r2 = gF*g2F + g2F*B1 + gF*A1
        kC = up((r1*(1 - A2) + 8*A2*r2)/dt*(1 + IOT)); kF = up(((1 - 4*Pv*Y2)*r2 + Pv*Y2*r1)/dt*(1 + IOT))
        ok = [bool(gF - (Pv*gC*Y1 + B0 + gF*A1) > 0), bool(gC - (4*Pv*gC*Y1 + 8*B0 + 8*gF*A1) > 0),
              bool(g2F - (Pv*g2C*Y1 + A0_ + g2F*A1) > 0), bool(g2C - (4*Pv*g2C*Y1 + 8*A0_ + 8*g2F*A1) > 0),
              bool(kF - (gF*g2F + Pv*kC*Y2 + g2F*B1 + gF*A1 + kF*A2) > 0), bool(kC - (gC*g2C + 4*Pv*kC*Y2 + 8*(g2F*B1 + gF*A1 + kF*A2)) > 0)]
        return (gF, gC, g2F, g2C, kF, kC), ok
    Ya = [Ys(b, xC, j) for j in range(3)]; Aa = [As(a, xF, j) for j in range(3)]; Ba = [Bs(a, xF, j) for j in range(2)]
    Yb = [Ys(bb, xCb, j) for j in range(3)]; Abr = [As(ab, xFb, j) for j in range(3)]; Bbr = [Bs(ab, xFb, j) for j in range(2)]
    mk_a, ok_a = marked(Ya, Aa, Ba); mk_b, ok_b = marked(Yb, Abr, Bbr)
    res['six inequalities (actual)'] = (all(ok_a), arb(0)); res['six inequalities (barred)'] = (all(ok_b), arb(0))
    gF, gC, g2F, g2C, kF, kC = mk_a
    # --- the sums of the hypotheses ---
    f = lambda k: fval(k, a)
    fk = {k: f(k) for k in range(4, 4*NX + 8)}
    def X(n, k, c=4): return arb(0) if k == 5 else comb(4*n + c, k).min((4*n + c)*Nv[k])
    z = xF.exp()*a
    def ex_q(c):
        L = 4*NX + c
        q0 = Fm*a**-4*arb(L)**2*comb(L, 5)*z**5*(1 + z)**(L - 5)
        qr = (arb(L + 4)/L)**2*(arb(L + 4)/(L - 1))**5*(1 + z)**4
        return q0, qr
    def crude_ex(n, c=4): return 2*comb((4 if c == 4 else 6)*Pv + Tv*n, n)*b**n*sum((X(n, k, c)*fk[k] for k in range(6, 4*n + c + 1)), arb(0))
    def crude_4(n): return 2*comb(4*Pv + Tv*n, n)*(n + 1)*(2*d - 3)*b**n
    q0, qr = ex_q(4); tail_ex = 2*poly_tail(4*Pv, b, q0, qr)
    tail_4 = 2*(2*d - 3)*poly_tail(4*Pv, b, arb(NX + 1), arb(NX + 2)/(NX + 1))
    mark4 = 4*Pv*(4*xF).exp()*(kC*Ya[2] + 4*(gC*g2F + g2C*gF)*Ya[1] + 16*kF*Ya[0])
    full = (b + (4*(Pv - 1) - KS)*b**2 + 48*comb(18, 2)*(d - 2)*(2*d - 3)*b**3 + sum((crude_4(n) for n in range(4, NX)), arb(0)) + tail_4
            + sum((crude_ex(n) for n in range(1, NX)), arb(0)) + tail_ex + mark4)
    eta = C*a**4*(2*Uhi*(1 + dl) + full)/(1 - om)
    Z4 = dict(T.ZSMALL[d])
    if d == 3:
        four = arb(0)                         # Z_(n,4) = 0 for all n >= 3 in d = 3, infinite volume (Lemma mg:lem:planar(i)); (d-3) = 0
    else:
        four = (10*a**4 + 2*eta*b**4)*(d - 3) + sum((arb(c)*b**n for n, c in Z4.items()), arb(0)) \
               + sum((crude_4(n) for n in range(max(Z4) + 1, NX)), arb(0)) + tail_4
    ex = KS*(a*fk[6] + b*fk[7]) + (4*(Pv - 1) - 2*KS)*b*fk[8] + b**2*sum((arb(c)*fk[k] for k, c in T.ZEXIT[d].items()), arb(0))
    for n in range(3, NX):
        if n in T.ZEXITN.get(d, {}): ex += b**n*sum((arb(int(c))*fk[k] if isinstance(c, int) else A(c)*fk[k] for k, c in T.ZEXITN[d][n].items()), arb(0))
        else: ex += crude_ex(n)
    ex += tail_ex
    w_a = a.sqrt()
    Fk = lambda kmin: 8*(sum((Nv[k]*fk[k] for k in range(kmin, KX)), arb(0)) + Fm*a**-4*front_tail(a*xF.exp(), 2))
    Ak = lambda kmin: 8*cJ.exp()*(sum((Nv[k]*a**(arb(k)/2) for k in range(kmin, KX)), arb(0)) + front_tail(w_a, 0))
    AA = Ak(4)
    onep = 4*Pv*b*(Fk(8)*AA + Fk(4)*Ak(8) + Fk(4)*(AA.expm1() - AA))
    two = (4*(Pv - 1) - KS)*varpi*(b**2).min(b2)
    gFb, gCb, g2Fb, g2Cb, kFb, kCb = mk_b
    S3 = 4*Pv*(4*xFb).exp()*(kCb*Yb[2] + 4*(gCb*g2Fb + g2Cb*gFb)*Yb[1] + 16*kFb*Yb[0])
    exact_t = sum((arb(int(c))*b**n*a**A(str(E_))*(m*cJ).exp()*(fk[kp]/a**(kp - 4) if kp > 4 else 1) for (c, n, E_, kp, m) in T.EXACT[d]), arb(0))
    cf = None; cmax = None
    for (cn, ce) in T.CORNERS[d]:
        v = (b/bb)**cn*(a/ab)**A(str(ce))
        if cf is None or bool(v > cf): cmax = (cn, ce)
        cf = v if cf is None else cf.max(v)
    multi = exact_t + cf*S3
    S4 = two + four + ex + onep + multi
    q06, qr6 = ex_q(6)
    sig6 = (2*a + sum(((2*d - 3)*(2*n + 3)*comb(6*Pv + Tv*n, n)*b**n for n in range(2, NX)), arb(0))
            + (2*d - 3)*poly_tail(6*Pv, b, arb(2*NX + 3), arb(2*NX + 5)/(2*NX + 3))
            + sum((crude_ex(n, 6) for n in range(1, NX)), arb(0)) + 2*poly_tail(6*Pv, b, q06, qr6)
            + 6*Pv*(6*xF).exp()*(kC*Ya[2] + 6*(gC*g2F + g2C*gF)*Ya[1] + 36*kF*Ya[0]))
    need = Pv*(gC*Ya[1]/7 + gF*Ya[0])
    check('Sigma^(4)_4 < Ubar delta\'', S4, Ulo*dl, res)
    check('F-condition', need, F*(Ulo*(1 + dl) - a), res)
    check('F_6-condition', sig6, 36*F6*(6*xF).exp()*(Ulo*(1 + dl) - a), res)
    check('third hypothesis (omega)', om + C*a**4*(2*Uhi*(1 + dl) + full)/(1 - om), arb(1), res)
    check('a^4 Ubar (1+delta\') <= varpi', a**4*Uhi*(1 + dl), varpi, res)
    check('a^4 Ubar (1+delta\') <= 1', a**4*Uhi*(1 + dl), arb(1), res)
    kap = beta                                                    # Lemma mg:lem:scaling, kappa = beta_*
    check('stronger: F-condition, e^{-2b*} - a', need, F*(Ulo - a), res)
    check('stronger: F_6-condition, e^{-2b*} - a', sig6, 36*F6*(6*xF).exp()*(Ulo - a), res)
    check('stronger: third hypothesis <= 1/(1+kappa)', om + C*a**4*(2*Uhi*(1 + dl) + full)/(1 - om), 1/(1 + kap), res)
    check('stronger: a^4 Ubar (1+delta\') <= varpi/(1+kappa)', a**4*Uhi*(1 + dl), varpi/(1 + kap), res)
    check('delta\'_* <= kappa', dl, kap, res)
    res['_delta\'_*'] = (True, dl)
    res['_majorant/beta^4 (corner ' + str(cmax) + ')'] = (True, cf*S3/beta**4)
    res['_exact pairs/beta^4'] = (True, exact_t/beta**4)
    res['_Sigma4^(4)/beta^4'] = (True, S4/beta**4)
    res['_parts/beta^4: two, four, exits, one-plaquette, multi'] = (True, [v/beta**4 for v in (two, four, ex, onep, multi)])
    res['_D parts/beta^4: 8(K+1)theta, e^{2b}Delta, log(1+delta\')'] = (True, [v/beta**4 for v in (8*(K + 1)*th, (2*beta).exp()*Delta, dl.log1p())])
    Dd = 8*(K + 1)*th + (2*beta).exp()*Delta + dl.log1p()
    cbound = (1 + kap)**3*Dd/beta**4
    return P, res, Dd/beta**4, cbound

if __name__ == "__main__":
    G, d, cminus = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    P, res, Dn, cb = run(G, d)
    name = {2: 'SU(2) (N = 2, iota = 1)', 1: 'U(1) (N = 1, iota = 0)', 3: 'SU(3) (N = 3, iota = 0)', 4: 'SU(4) (N = 4, iota = 0)', 5: 'SU(5) (N = 5, iota = 0)'}[G]
    if DECISION_D3 and d == 3: name += ' [Lemma mg:lem:frontsdecision for Z^2 beyond 12 edges (--nv-decision-d3)]'
    print(f"{name}, d={d}, set {SETTINGS[d][0]}, varpi = b^{float(P['p']):g}, (bbar, abar) = ({P['bb']}, {P['ab']}), KMAX = {T.KMAX_DEFAULT}, F_RAISE = {mpnstr(T.F_RAISE, 3)};"
          f" beta_* = {mpnstr(P['beta'], 15)} (checked >= beta_0(d) in Arb); {PREC}-bit balls, sums to {KX} edges / {NX} plaquettes plus rigorous remainders:")
    allok = True
    for k, (ok, room) in res.items():
        if k.startswith('_'):
            v = room
            txt = ', '.join(x_.mid().str(6, radius=False) for x_ in v) if isinstance(v, list) else v.mid().str(8, radius=False)
            print(f"    {k[1:]:60s} {txt}"); continue
        allok = allok and ok
        print(f"    {k:60s} {'holds' if ok else 'NOT CERTIFIED'}" + (f"   (room {room.mid().str(4, radius=False)})" if room != 0 else ''))
    ok_c = bool(arb(cminus) - cb > 0)
    print(f"  authors' stronger-form rooms (mpmath): " + ', '.join(f"{k} {mpnstr(v, 4)}" for k, v in P['strong_auth'].items()))
    print(f"  D(beta_*)/beta_*^4 = {Dn.str(15, radius=False)} (authors: {mpnstr(P['D_auth']/P['beta']**4, 12)});  (1+kappa)^3 D/beta_*^4 <= {arb(cb.upper()).str(10, radius=False)} <= {cminus}: {ok_c}")
    print(f"  all certified: {allok and ok_c}")

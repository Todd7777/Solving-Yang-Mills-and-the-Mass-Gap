# Referee check (TENTH pass on Section 16, Corollary mg:cor:sectorscap), item 7: an adapted copy of the ninth referee's referee-section16i/arb_eval.py.
# The evaluation at beta_* of the hypotheses of Theorem mg:thm:lower5 in bound mode (Lemma mg:lem:boundmode), in BALL ARITHMETIC (Arb, python-flint, 400 bits),
# every series summed to KX/NX with its remainder bounded rigorously, the referee's own marked constants (six inequalities verified in Arb), and the
# authors' parameters (varpi, xi_F, xi_C, h_J, F, F_6, delta', bbar, abar) as exact binary numbers: see the ninth referee's header for the method.
# Changes for the corollary (marked "COR" below):
#   - the set S and its corners are read from a given JSON file: orders/d3_s9z.json (d=3, the set S''_9 u Z of the corollary), orders/d3_s9.json (d=3, the
#     set of Theorem mg:thm:thirdcap, for comparison) or orders/d4_s10x_92.json (d=4); the (bbar, abar) grid point is the one the authors' search picks
#     (orders/capz_*_scaling.out, cap_su4_d4_scaling.out, cap_*_scaling.out);
#   - GROUP = 3, 4 or 5 (for N = 4, Lemma mg:lem:u2sun with c_Lambda = 2; iota = 0 for all three);
#   - the STRONGER FORM of Lemma mg:lem:scaling is checked explicitly in Arb: the conditions on F and F_6 with e^{-2beta_*} - ahat in place of
#     e^{-2beta_*}(1+delta'_*) - ahat; the third hypothesis with 1/(1+kappa) in place of 1; the fourth with varpi_*/(1+kappa) in place of varpi_*; and delta'_* <= kappa;
#   - (1+kappa)^3 D(beta_*)/beta_*^4 is compared with the corollary's c'_- (or Theorem mg:thm:thirdcap's c_-), and the corner majorant
#     max_C (b/bbar)^n (a/abar)^E x (barred tree sum) of Lemma mg:lem:corners is reported (the corollary's "majorant").
# Usage: python3 arb_eval_cor.py GROUP d SET CMINUS      (SET = s9z, s9 (d=3) or x92 (d=4); e.g. python3 arb_eval_cor.py 4 3 s9z 0.134)
import os, sys, math
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..'))
from flint import arb, ctx
from mpmath import mp, nstr as mpnstr
import my_eval as ME
import third_order_effective as T

PREC = 400
KX, NX = 200, 110
ctx.prec = PREC
ONE = arb(1); E = ONE.exp()

def A(x):
    """an exact arb from an mpmath number, an int, a Fraction-like or a decimal string"""
    if isinstance(x, int): return arb(x)
    if isinstance(x, str): return arb(arb(x).mid())
    return arb(arb(mpnstr(x, 70)).mid())

def comb(n, k): return arb(math.comb(n, k)) if 0 <= k <= n else arb(0)

def check(name, L, R, store):
    """L < R for the whole balls"""
    ok = bool((R - L) > 0)
    store[name] = (ok, (1 - L/R))
    return ok

SETS = {'s9z': 'd3_s9z.json', 's9': 'd3_s9.json', 'x92': 'd4_s10x_92.json'}                                   # COR
ABARS = {(3, 3, 's9z'): 2**-10, (4, 3, 's9z'): 2**-10, (5, 3, 's9z'): 2**-10, (4, 4, 'x92'): 2**-10,          # COR: the authors' grid choices
         (3, 3, 's9'): 2**-11, (5, 3, 's9'): 2**-12, (3, 4, 'x92'): 2**-9, (5, 4, 'x92'): 2**-10}

def authors_params(G, d, setname):
    """COR: as my_eval.authors_params, with the terms file and abar of the case"""
    from mpmath import mpf
    PB = ME.PB
    PB.GROUP = G; PB.BOUNDS = True; PB.CERTIFIED = True
    _, closed, cf, bbs, bstar = ME.SETTINGS[d]
    terms = SETS[setname]
    ME.SETTINGS[d] = (terms, closed, cf, bbs, bstar)
    PB.load_terms(d, os.path.join(here, '..', 'orders', terms), closed, cf)
    beta = mpf(bstar); g = PB.group_data()(beta)
    lo = T.lower3e(beta, d, g, varpi_pow=mpf(5)/2, chars=True, bbar_grid=list(bbs), abar_grid=[mpf(ABARS[(G, d, setname)])])
    return dict(beta=beta, varpi=lo['varpi'], xF=lo['xF'], xC=lo['xC'], hJ=lo['hJ'], F=lo['F'], F6=lo['F6'], dl=lo['delta'], bb=lo['bars'][0],
                ab=lo['bars'][1], D_auth=lo['D_bound'], parts_auth=lo['parts'], sig6_auth=lo['sig6'], g=g)

def run(G, d, setname):
    mp.dps = 60
    P = authors_params(G, d, setname)           # parameters at beta_* (bound mode, certified counts), and the count tables in T
    beta = A(P['beta'])
    K = 4*(2*d - 3)
    beta0 = (ONE + ONE/(E**2*K*(K + 2))).log()/2
    assert bool(beta - beta0 > 0), 'beta_* >= beta_0(d)'
    n_ = d - 1; Pv = 2*n_*(n_ - 1); KS = 4*(2*d - 5); Tv = 4*(Pv - 1); N = G; iota = 0
    # --- group bounds (Lemmas mg:lem:char(iii), mg:lem:u2sun, mg:lem:capelement(iv), mg:lem:dressed, mg:lem:upsbounds(i)) ---
    ubar = ((1 + iota)*beta + beta**2*beta.exp())/(2*N**2)
    ulow = (1 + iota)*beta*(-beta).exp()/(2*N**2)
    R = beta.expm1() - beta - beta**2/2
    dS, dA, dL = arb(N*(N + 1))/2, arb(N*N - 1), arb(N*(N - 1))/2
    u2 = ((beta**2/(8*N**2) + R/dS)/dS).max((beta**2/(4*N**2) + R/dA)/dA).max(R)
    if N >= 4: u2 = u2.max(((2 if N == 4 else 1)*beta**2/(8*N**2) + R/dL)/dL)
    eps = (2*beta).expm1(); gam = eps + beta**2*beta.exp()
    s6 = (2 - iota)*N**2*ubar**6 + u2**5*beta.expm1()
    nB = (u2/ulow)**2*eps + arb(2 - iota).sqrt()*(N/ulow)*(u2**5*beta.expm1()).sqrt()
    c2 = (u2/ulow)**4*beta.exp()
    Ulo = (-2*beta).exp(); Uhi0 = beta.exp()
    assert bool(Ulo - c2 > 0) and bool(Ulo - ubar**4*(1 + gam) > 0), 'denominators of Lemma mg:lem:dressed'
    Delta = nB**2/(Ulo - c2) + 2*(2*beta).exp()*s6/(Ulo - ubar**4*(1 + gam))
    Uhi = Uhi0 + Delta
    # --- theta (Lemma mg:lem:theta7 with the closed polymers counted exactly), a, b ---
    closed, CF = ME.SETTINGS[d][1], ME.SETTINGS[d][2]
    x = E**2*K*eps
    th = 2*(d - 2)*E**6*s6 + sum((arb(c)*(E*eps)**n for n, c in closed.items()), arb(0)) + x**CF/(1 - x)
    lam = (2*(K + 1)*th).exp()
    a = lam*ubar; b = lam*gam
    varpi, xF, xC, hJ, F, F6, dl = (A(P[k]) for k in ('varpi', 'xF', 'xC', 'hJ', 'F', 'F6', 'dl'))
    bb, ab = A(P['bb']), A(P['ab'])
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
    KK = 4*NX + 8
    for k in range(6, KK + 1):
        crude = n_*comb((4*d - 6)*k, k - 1)
        if k in exact: Nv[k] = arb(exact[k])
        elif d == 4 and k < 60: Nv[k] = crude.min(min((dec(s, k) for s in range(1, 400)), key=lambda v: v.upper()))
        elif d == 4 or k >= 60: Nv[k] = crude.min(dec(s_fix, k))
        else: Nv[k] = crude
    s0 = arb(s_fix)/400; y0 = ((1 + s0)**(2*D_ - 1) - 1).max((1 + s0)**(D_ - 1)); C0 = (1 + s0)**(2*D_) - 1 - 2*D_*s0; Y0r = y0/s0
    def front_tail(w, j, Kc=KX):
        """sum_{k >= Kc} k^j N^v_k w^k <= C0 Kc^j r^Kc / (1 - (1+1/Kc)^j r), r = (y/s) w (Lemma mg:lem:frontsdecision)"""
        r = Y0r*w; q = (1 + ONE/Kc)**j*r
        assert bool(1 - q > 0)
        return C0*arb(Kc)**j*r**Kc/(1 - q)
    def cluster_tail(y, p, M=NX):
        """sum_{m >= M} m^p binom(K_S m, m-1) y^m <= e y M^p (e K_S y)^(M-1) / (1 - (1+1/M)^p e K_S y)"""
        r = E*KS*y; q = (1 + ONE/M)**p*r
        assert bool(1 - q > 0)
        return E*y*arb(M)**p*r**(M - 1)/(1 - q)
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
    # --- decoration conditions (the authors' xi_F, xi_C, h_J) ---
    Y0 = Ys(b, xC, 0); A0 = As(a, xF, 0)
    check('deco xi_F', Pv*Y0 + A0, xF, res); check('deco xi_C', 4*Pv*Y0 + 8*A0, xC, res)
    w4 = a.sqrt()*xF.exp()
    assert bool(1 - 2*w4 > 0)       # k w^k decreasing for k >= 4 when w <= 1/2... (k+1)w^(k+1)/(k w^k) = (1+1/k) w <= 5/4 w < 1
    check('deco h_J (k=4, the maximum)', C*4*Pv*Y0*a**2*(4*xF).exp(), hJ, res)
    # --- barred: the authors' xi-bar (tree_constants at their c_J), checked here ---
    import third_order_constants as TC
    from mpmath import mpf
    mp.dps = 60
    tc = TC.tree_constants(P['bb'], P['ab'], mpf(mpnstr(cJ.mid(), 60)), Pv, KS, {k: mpf(mpnstr(Nv[k].upper(), 60)) for k in range(4, 200)}, 60, 40)
    xFb, xCb = A(tc[0]), A(tc[1])
    Yb0 = Ys(bb, xCb, 0); Ab0 = As(ab, xFb, 0)
    check('barred xi_F', Pv*Yb0 + Ab0, xFb, res); check('barred xi_C', 4*Pv*Yb0 + 8*Ab0, xCb, res)
    # --- marked constants (referee's construction) ---
    IOT = arb('1e-30')
    def up(x): return arb(x.upper())
    def marked(Y, Av, B):
        Y0_, Y1, Y2 = Y; A0_, A1, A2 = Av; B0, B1 = B
        det = (1 - A1)*(1 - 4*Pv*Y1) - Pv*Y1*8*A1
        assert bool(det > 0) and bool(1 - A1 > 0) and bool(1 - 4*Pv*Y1 > 0)
        g = [(B0*(1 - 4*Pv*Y1) + Pv*Y1*8*B0)/det, ((1 - A1)*8*B0 + 8*A1*B0)/det,
             (A0_*(1 - 4*Pv*Y1) + Pv*Y1*8*A0_)/det, ((1 - A1)*8*A0_ + 8*A1*A0_)/det]
        gF, gC, g2F, g2C = [up(v*(1 + IOT)) for v in g]
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
    # exits majorant beyond NX: q_n <= Fm a^-4 L^2 binom(L,5) z^5 (1+z)^(L-5), ratio <= ((L+4)/L)^2 ((L+4)/(L-1))^5 (1+z)^4 (L = 4n+c >= 16)
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
        four = arb(0)                         # Z_(n,4) = 0 for all n >= 3 in d = 3, infinite volume (Remark mg:rem:nextorder(iv)); (d-3) = 0
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
    cf = None
    for (cn, ce) in T.CORNERS[d]:
        v = (b/bb)**cn*(a/ab)**A(str(ce))
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
    kap = beta                                                                                                    # COR: Lemma mg:lem:scaling
    check('stronger: F-condition, e^{-2b*} - a', need, F*(Ulo - a), res)
    check('stronger: F_6-condition, e^{-2b*} - a', sig6, 36*F6*(6*xF).exp()*(Ulo - a), res)
    check('stronger: third hypothesis <= 1/(1+kappa)', om + C*a**4*(2*Uhi*(1 + dl) + full)/(1 - om), 1/(1 + kap), res)
    check('stronger: a^4 Ubar (1+delta\') <= varpi/(1+kappa)', a**4*Uhi*(1 + dl), varpi/(1 + kap), res)
    check('delta\'_* <= kappa', dl, kap, res)
    res['_majorant/beta^4'] = (True, cf*S3/beta**4)
    res['_exact pairs/beta^4'] = (True, exact_t/beta**4)
    res['_Sigma4^(4)/beta^4'] = (True, S4/beta**4)
    res['_parts/beta^4: two, four, exits, one-plaquette, multi'] = (True, [v/beta**4 for v in (two, four, ex, onep, multi)])
    Dd = 8*(K + 1)*th + (2*beta).exp()*Delta + dl.log1p()
    cbound = (1 + kap)**3*Dd/beta**4
    return P, res, Dd/beta**4, cbound

if __name__ == "__main__":
    G, d, setname, cminus = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
    P, res, Dn, cb = run(G, d, setname)
    print(f"SU({G}), d={d}, set {SETS[setname]}, (bbar, abar) = ({P['bb']}, {P['ab']}); beta_* = {mpnstr(P['beta'], 15)} (checked >= beta_0(d) in Arb); {PREC}-bit balls,"
          f" sums to {KX}/{NX} plus rigorous remainders:")
    allok = True
    for k, (ok, room) in res.items():
        if k.startswith('_'):
            v = room
            txt = ', '.join(x.mid().str(6, radius=False) for x in v) if isinstance(v, list) else v.mid().str(8, radius=False)
            print(f"    {k[1:]:52s} {txt}"); continue
        allok = allok and ok
        print(f"    {k:52s} {'holds' if ok else 'NOT CERTIFIED'}" + (f"   (room {room.mid().str(4, radius=False)})" if room != 0 else ''))
    ok_c = bool(arb(cminus) - cb > 0)
    print(f"  D(beta_*)/beta_*^4 = {Dn.str(15, radius=False)} (authors: {mpnstr(P['D_auth']/P['beta']**4, 12)});  (1+kappa)^3 D/beta_*^4 <= {arb(cb.upper()).str(10, radius=False)} <= {cminus}: {ok_c}")
    print(f"  all certified: {allok and ok_c}")

# (Tenth referee: a verbatim copy of referee-section16i/my_eval.py, the ninth referee's evaluator; used here only for SETTINGS, PB and the group bounds.)
# Referee check (ninth pass on Section 16) of Theorem mg:thm:thirdcap(a)-(c): an evaluator of the hypotheses of Theorem mg:thm:lower5 in bound mode,
# written from the formulas of the paper (Definition mg:def:jconst, Lemmas mg:lem:jsums, mg:lem:jsums4, mg:lem:jsums6, mg:lem:corners, Theorems
# mg:thm:lower2, mg:thm:lower4, mg:thm:lower5, Lemma mg:lem:dressed, Lemma mg:lem:u2sun, Lemma mg:lem:capelement(iv), Lemma mg:lem:char(iii),
# Lemma mg:lem:upsbounds(i)), with every infinite series summed far beyond the authors' truncation (fronts to KX edges, clusters and polymers to NX
# plaquettes; the rest is below 1e-100 relative), and with its own solutions of the six inequalities of Definition mg:def:jconst.
#
# It uses the authors' code only (i) to obtain the parameters at beta_* (third_order_effective.lower3e, as scaling_check.params does), and (ii) to read
# the count tables (practical_beta.load_terms). It then:
#   1. checks the authors' decoration constants, barred tree constants and marked constants (third_order_constants.marked) against their defining
#      inequalities, with the extended sums;
#   2. computes its own valid marked constants (minimal solutions with the extended sums, raised correctly: the kappa's are solved with the raised g's);
#   3. evaluates every hypothesis and the deficit D with its own constants, and reports the room of each.
# Usage: python3 my_eval.py GROUP d          (GROUP = 3 or 5, d = 3 or 4; about 1-2 min)
import os, sys, json
here = os.path.dirname(os.path.abspath(__file__)); parent = os.path.join(here, '..')
sys.path.insert(0, parent)
import practical_beta as PB
import third_order_effective as T
import third_order_constants as TC
from mpmath import mp, mpf, e, exp, log, expm1, binomial, fsum, sqrt, nstr

KX, NX = 260, 130          # extended truncations (fronts, clusters/polymers)
IOTA = mpf(10)**-30        # raise of my own constants (covers the remainders beyond KX, NX, which are < 1e-100)

B3 = '0.00040262143352281974401285396895067832618294529962141'   # beta_0(3)(1+1e-30)
B4 = '0.00015376644804254827379172376907999414111262406766301'   # beta_0(4)(1+1e-30)
SETTINGS = {3: ('d3_s9.json', {10: 10, 11: 11, 12: 24, 13: 0, 14: 70}, 15, (mpf(2)**-7,), B3),
            4: ('d4_s10x_92.json', {10: 60, 11: 66, 12: 320, 13: 0, 14: 1316, 15: 2520, 16: 12552}, 17, (mpf(2)**-8,), B4)}
ABAR = {(3, 3): mpf(2)**-11, (5, 3): mpf(2)**-12, (3, 4): mpf(2)**-9, (5, 4): mpf(2)**-10}   # the (bbar, abar) the authors' search picks (cap_*.out)

# ---------------- the group data in bound mode, from the lemmas (my implementation) ----------------
def group_bounds(N, beta, iota=0):
    beta = mpf(beta)
    ubar = ((1 + iota)*beta + beta**2*exp(beta))/(2*N**2)                       # Lemma mg:lem:char(iii)
    ulow = (1 + iota)*beta*exp(-beta)/(2*N**2)
    R = expm1(beta) - beta - beta**2/2                                          # Lemma mg:lem:u2sun
    dS, dA, dL = mpf(N*(N + 1))/2, mpf(N*N - 1), mpf(N*(N - 1))/2
    cands = [(beta**2/(8*N**2) + R/dS)/dS, (beta**2/(4*N**2) + R/dA)/dA, R]
    if N >= 4: cands.append(((2 if N == 4 else 1)*beta**2/(8*N**2) + R/dL)/dL)
    u2 = max(cands)
    eps = expm1(2*beta); gam = expm1(2*beta) + beta**2*exp(beta)                  # gamma_g <= gamma_2
    s6 = (2 - iota)*N**2*ubar**6 + u2**5*expm1(beta)                             # Lemma mg:lem:capelement(iv)
    nB = (u2/ulow)**2*eps + sqrt(2 - iota)*(N/ulow)*sqrt(u2**5*expm1(beta))      # Lemma mg:lem:dressed, n_B with u2bar, ulow
    c2 = (u2/ulow)**4*exp(beta)
    lo, hi = exp(-2*beta), exp(beta)                                             # Lemma mg:lem:upsbounds(i)
    assert lo > c2 and lo > ubar**4*(1 + gam)
    Delta = nB**2/(lo - c2) + 2*exp(2*beta)*s6/(lo - ubar**4*(1 + gam))
    return dict(ubar=ubar, ulow=ulow, u2=u2, s6=s6, gam=gam, eps=eps, Delta=Delta, Ulo=lo, Uhi=hi + Delta)

# ---------------- counts ----------------
def front_counts(d, KMAXX):
    """N^v_k as the authors use them for k < 60 (exact where counted; d=3: (d-1) binom((4d-6)k,k-1); d=4: min with the decision bound over the
    1/400 grid), and for k >= 60 the bound of Lemma mg:lem:frontsdecision with s fixed (0.194 in Z^3, 0.325 in Z^2), valid for every k."""
    n = d - 1; Pv = 2*n*(n - 1)
    Nv = {4: mpf(Pv), 5: mpf(0)}
    exact = dict(T.NOVR[d])
    D = d - 1; s0 = {2: mpf('0.325'), 3: mpf('0.194')}[D]
    y0 = max((1 + s0)**(2*D - 1) - 1, (1 + s0)**(D - 1)); C0 = (1 + s0)**(2*D) - 1 - 2*D*s0
    for k in range(6, KMAXX + 1):
        crude = n*binomial((4*d - 6)*k, k - 1)
        if k in exact: Nv[k] = mpf(exact[k])
        elif k < 60:
            Nv[k] = crude if d == 3 else min(crude, T.decision_bound(D, k))
        else:
            Nv[k] = min(crude, C0*(y0/s0)**k)
    return Nv

def evaluate(G, d, P, verbose=True):
    """every hypothesis of Theorem mg:thm:lower5 at beta_* with the parameters P (from the authors' lower3e), my sums and my marked constants"""
    beta = P['beta']; gb = group_bounds(G, beta)
    K = 4*(2*d - 3); eps = gb['eps']; x = e**2*K*eps
    n_ = d - 1; Pv = 2*n_*(n_ - 1); KS = 4*(2*d - 5); Tv = 4*(Pv - 1)
    closed, CF = SETTINGS[d][1], SETTINGS[d][2]
    th = 2*(d - 2)*e**6*gb['s6'] + fsum(c*(e*eps)**n for n, c in closed.items()) + x**CF/(1 - x)   # Lemma mg:lem:theta7 with exact closed polymers
    lam = exp(2*(K + 1)*th)
    a = lam*gb['ubar']; b = lam*gb['gam']
    varpi, xF, xC, hJ, F, F6, dl = P['varpi'], P['xF'], P['xC'], P['hJ'], P['F'], P['F6'], P['dl']
    bb, ab = P['bb'], P['ab']
    Ulo, Uhi = gb['Ulo'], gb['Uhi']
    C = 1 + max(1/varpi, (1 + gb['gam'])/Ulo)
    om = C*a**4*(1 + gb['gam'])
    cJ = exp(hJ)/(1 - om)
    b2 = lam**2*(2*gb['s6'] + gb['s6']**2)/gb['ulow']**4
    Nv = front_counts(d, 4*NX + 8)
    out = dict(a=a, b=b, theta=th, cJ=cJ, omega=om, C=C)
    # --- series (extended) ---
    def Ys(bc, z, j, M=NX): return fsum(m**j*binomial(KS*m, m - 1)*(bc*exp(z))**m for m in range(1, M))
    def As(aa, xf, j, cj, M=KX): return exp(cj)*fsum(k**j*Nv[k]*aa**(mpf(k)/2)*exp(xf*k) for k in range(4, M))
    def fval(k, aa, F_, F6_, xf):     # f(k) with front weight aa (f-bar: aa = abar, the actual xi_F inside f)
        if k == 4: return mpf(1)
        if k == 5: return mpf(0)
        if k == 6: return 36*F6_*exp(6*xf)*aa**2
        return F_*k**2*exp(xf*k)*aa**(k - 4)
    def Bs(aa, xfB, j, M=KX): return fsum(k**j*Nv[k]*fval(k, aa, F, F6, xF)*exp(xfB*k) for k in range(4, M))
    # --- 1. decoration conditions (Lemma mg:lem:jsums with C_varpi; Theorem mg:thm:lower2), the authors' xi_F, xi_C, h_J ---
    Y0 = Ys(b, xC, 0); A0 = As(a, xF, 0, cJ)
    hreq = max(C*k*Pv*Y0*a**(mpf(k)/2)*exp(xF*k) for k in range(4, 200))
    rooms = dict(deco_xF=1 - (Pv*Y0 + A0)/xF, deco_xC=1 - (4*Pv*Y0 + 8*A0)/xC, deco_hJ=1 - hreq/hJ)
    # --- barred tree constants (authors') ---
    tc = TC.tree_constants(bb, ab, cJ, Pv, KS, {k: Nv[k] for k in range(4, 200)}, 60, 40)
    xFb, xCb = tc[0], tc[1]
    Yb0 = Ys(bb, xCb, 0); Ab0 = As(ab, xFb, 0, cJ)
    rooms.update(bar_xF=1 - (Pv*Yb0 + Ab0)/xFb, bar_xC=1 - (4*Pv*Yb0 + 8*Ab0)/xCb)
    # --- the six inequalities ---
    def six_lhs(g, Y, A, B):
        gF, gC, g2F, g2C, kF, kC = g
        return (Pv*gC*Y[1] + B[0] + gF*A[1], 4*Pv*gC*Y[1] + 8*B[0] + 8*gF*A[1],
                Pv*g2C*Y[1] + A[0] + g2F*A[1], 4*Pv*g2C*Y[1] + 8*A[0] + 8*g2F*A[1],
                gF*g2F + Pv*kC*Y[2] + g2F*B[1] + gF*A[1] + kF*A[2],
                gC*g2C + 4*Pv*kC*Y[2] + 8*(g2F*B[1] + gF*A[1] + kF*A[2]))
    def six_rooms(g, Y, A, B):
        L = six_lhs(g, Y, A, B)
        return [1 - l/r for l, r in zip(L, (g[0], g[1], g[2], g[3], g[4], g[5]))]
    def six_rooms_G(g, Y, A, B):
        """(5), (6) in the form Lemma mg:lem:marked needs: G_1 G_2 (the left sides of (1),(3), resp. (2),(4)) in place of g_F g'_F, g_C g'_C"""
        L = six_lhs(g, Y, A, B)
        L5 = L[4] - g[0]*g[2] + L[0]*L[2]; L6 = L[5] - g[1]*g[3] + L[1]*L[3]
        return [1 - L5/g[4], 1 - L6/g[5]]
    def my_marked(Y, A, B, iota=IOTA):
        """minimal solution of the six inequalities, made strict correctly: g's raised by (1+iota), the kappa's solved with the raised g's, then raised"""
        Y0_, Y1, Y2 = Y; A0_, A1, A2 = A; B0, B1 = B
        det = (1 - A1)*(1 - 4*Pv*Y1) - Pv*Y1*8*A1
        gF = (B0*(1 - 4*Pv*Y1) + Pv*Y1*8*B0)/det; gC = ((1 - A1)*8*B0 + 8*A1*B0)/det
        g2F = (A0_*(1 - 4*Pv*Y1) + Pv*Y1*8*A0_)/det; g2C = ((1 - A1)*8*A0_ + 8*A1*A0_)/det
        gF, gC, g2F, g2C = [(1 + iota)*v for v in (gF, gC, g2F, g2C)]
        dt = (1 - 4*Pv*Y2)*(1 - A2) - 8*A2*Pv*Y2
        r1 = gC*g2C + 8*(g2F*B1 + gF*A1); r2 = gF*g2F + g2F*B1 + gF*A1
        kC = (r1*(1 - A2) + 8*A2*r2)/dt; kF = ((1 - 4*Pv*Y2)*r2 + Pv*Y2*r1)/dt
        return (gF, gC, g2F, g2C, (1 + iota)*kF, (1 + iota)*kC)
    # actual (b, a, f), and barred (bbar, abar, fbar)
    Ya = [Ys(b, xC, j) for j in range(3)]; Aa = [As(a, xF, j, cJ) for j in range(3)]; Ba = [Bs(a, xF, j) for j in range(2)]
    Yb = [Ys(bb, xCb, j) for j in range(3)]; Abr = [As(ab, xFb, j, cJ) for j in range(3)]; Bbr = [Bs(ab, xFb, j) for j in range(2)]
    # the authors' marked constants (their truncated sums, INFLATE), checked against the extended sums
    def f_auth(k):
        if k == 4: return mpf(1)
        if k == 6: return 36*F6*exp(6*xF)*a**2
        return F*k**2*exp(xF*k)*a**(k - 4)
    def fb_auth(k):
        if k == 4: return mpf(1)
        if k == 6: return 36*F6*exp(6*xF)*ab**2
        return F*k**2*exp(xF*k)*ab**(k - 4)
    N60 = {k: Nv[k] for k in range(4, 200)}
    Ya60 = [TC.fsum(m**j*binomial(KS*m, m - 1)*(b*exp(xC))**m for m in range(1, 40)) for j in range(3)]
    Aa60 = [exp(cJ)*fsum(k**j*N60[k]*a**(mpf(k)/2)*exp(xF*k) for k in range(4, 60)) for j in range(3)]
    auth_a = TC.marked(xF, Ya60, Aa60, f_auth, Pv, N60, 60)
    auth_b = TC.marked(xFb, tc[2], tc[3], fb_auth, Pv, N60, 60)
    rooms_auth_a = six_rooms(auth_a, Ya, Aa, Ba); rooms_auth_b = six_rooms(auth_b, Yb, Abr, Bbr)
    rooms_auth_G = six_rooms_G(auth_a, Ya, Aa, Ba) + six_rooms_G(auth_b, Yb, Abr, Bbr)
    # the same at the authors' own truncated sums (shows the effect of the raise alone, without the tails)
    Ba60 = [fsum(k**j*N60[k]*f_auth(k)*exp(xF*k) for k in range(4, 60)) for j in range(2)]
    Bb60 = [fsum(k**j*N60[k]*fb_auth(k)*exp(xFb*k) for k in range(4, 60)) for j in range(2)]
    rooms_auth_a_trunc = six_rooms(auth_a, Ya60, Aa60, Ba60); rooms_auth_b_trunc = six_rooms(auth_b, tc[2], tc[3], Bb60)
    # mine
    mk_a = my_marked(Ya, Aa, Ba); mk_b = my_marked(Yb, Abr, Bbr)
    rooms_my_a = six_rooms(mk_a, Ya, Aa, Ba); rooms_my_b = six_rooms(mk_b, Yb, Abr, Bbr)
    out.update(auth_a=auth_a, auth_b=auth_b, mk_a=mk_a, mk_b=mk_b)
    # --- Sigma_4 (the full sum), eta, Sigma^(4)_4, Sigma^(6), need ---
    gF, gC, g2F, g2C, kF, kC = mk_a
    def X(n, k): return mpf(0) if k == 5 else min(binomial(4*n + 4, k), (4*n + 4)*Nv[k])
    def X6(n, k): return mpf(0) if k == 5 else min(binomial(4*n + 6, k), (4*n + 6)*Nv[k])
    f = lambda k: fval(k, a, F, F6, xF)
    crude_ex = lambda n: 2*binomial(4*Pv + Tv*n, n)*b**n*fsum(X(n, k)*f(k) for k in range(6, 4*n + 5))
    crude_4 = lambda n: 2*binomial(4*Pv + Tv*n, n)*(n + 1)*(2*d - 3)*b**n
    mark4 = 4*Pv*exp(4*xF)*(kC*Ya[2] + 4*(gC*g2F + g2C*gF)*Ya[1] + 16*kF*Ya[0])
    full = (b + (4*(Pv - 1) - KS)*b**2 + 48*binomial(18, 2)*(d - 2)*(2*d - 3)*b**3 + fsum(crude_4(n) for n in range(4, NX))
            + fsum(crude_ex(n) for n in range(1, NX)) + mark4)
    eta = C*a**4*(2*Uhi*(1 + dl) + full)/(1 - om)
    Z4 = dict(T.ZSMALL[d])
    Z4max = max(Z4) if Z4 else 4
    four = (10*a**4 + 2*eta*b**4)*(d - 3) + fsum(c*b**n for n, c in Z4.items()) + fsum(crude_4(n) for n in range(Z4max + 1, NX))
    ex = KS*(a*f(6) + b*f(7)) + (4*(Pv - 1) - 2*KS)*b*f(8) + b**2*fsum(c*f(k) for k, c in T.ZEXIT[d].items())
    for n in range(3, NX):
        if n in T.ZEXITN.get(d, {}): ex += b**n*fsum(c*f(k) for k, c in T.ZEXITN[d][n].items())
        else: ex += crude_ex(n)
    Fk = lambda kmin: 8*fsum(Nv[k]*f(k) for k in range(kmin, KX))
    Ak = lambda kmin: 8*exp(cJ)*fsum(Nv[k]*a**(mpf(k)/2) for k in range(kmin, KX))
    AA = Ak(4)
    onep = 4*Pv*b*(Fk(8)*AA + Fk(4)*Ak(8) + Fk(4)*(exp(AA) - 1 - AA))
    two = (4*(Pv - 1) - KS)*varpi*min(b**2, b2)
    # M_S (Lemma mg:lem:corners) with the barred constants
    gFb, gCb, g2Fb, g2Cb, kFb, kCb = mk_b
    S3 = 4*Pv*exp(4*xFb)*(kCb*Yb[2] + 4*(gCb*g2Fb + g2Cb*gFb)*Yb[1] + 16*kFb*Yb[0])
    exact = fsum(c*b**n*a**mpf(E)*exp(m*cJ)*(f(kp)/a**(kp - 4) if kp > 4 else 1) for (c, n, E, kp, m) in T.EXACT[d])
    corner = max((b/bb)**cn*(a/ab)**mpf(ce) for (cn, ce) in T.CORNERS[d])*S3
    multi = exact + corner
    S4 = two + four + ex + onep + multi
    sig6 = (2*a + fsum((2*d - 3)*(2*n + 3)*binomial(6*Pv + Tv*n, n)*b**n for n in range(2, NX))
            + fsum(2*binomial(6*Pv + Tv*n, n)*b**n*fsum(X6(n, k)*f(k) for k in range(6, 4*n + 7)) for n in range(1, NX))
            + 6*Pv*exp(6*xF)*(kC*Ya[2] + 6*(gC*g2F + g2C*gF)*Ya[1] + 36*kF*Ya[0]))
    need = Pv*(gC*Ya[1]/7 + gF*Ya[0])
    hyp = dict(sigma=(S4, Ulo*dl), exits=(need, F*(Ulo*(1 + dl) - a)), six=(sig6, 36*F6*exp(6*xF)*(Ulo*(1 + dl) - a)),
               old=(om + C*a**4*(2*Uhi*(1 + dl) + full)/(1 - om), mpf(1)), beyond=(a**4*Uhi*(1 + dl), varpi), growth=(a**4*Uhi*(1 + dl), mpf(1)))
    for k_, (l, r) in hyp.items(): rooms['H_' + k_] = 1 - l/r
    D = 8*(K + 1)*th + exp(2*beta)*gb['Delta'] + log(1 + dl)
    out.update(rooms_auth_G=rooms_auth_G)
    out.update(rooms=rooms, rooms_auth_a=rooms_auth_a, rooms_auth_b=rooms_auth_b, rooms_auth_a_trunc=rooms_auth_a_trunc,
               rooms_auth_b_trunc=rooms_auth_b_trunc, rooms_my_a=rooms_my_a, rooms_my_b=rooms_my_b, D=D, S4=S4,
               parts=dict(two=two, four=four, exits=ex, one_plaq=onep, multi=multi), sig6=sig6, need=need, full=full, eta=eta, gb=gb)
    return out

def authors_params(G, d):
    """the parameters at beta_* as scaling_check.params finds them (bound mode, certified counts), with the authors' (bbar, abar) grid point"""
    PB.GROUP = G; PB.BOUNDS = True; PB.CERTIFIED = True
    terms, closed, cf, bbs, bstar = SETTINGS[d]
    PB.load_terms(d, os.path.join(parent, 'orders', terms), closed, cf)
    beta = mpf(bstar); g = PB.group_data()(beta)
    lo = T.lower3e(beta, d, g, varpi_pow=mpf(5)/2, chars=True, bbar_grid=list(bbs), abar_grid=[ABAR[(G, d)]])
    return dict(beta=beta, varpi=lo['varpi'], xF=lo['xF'], xC=lo['xC'], hJ=lo['hJ'], F=lo['F'], F6=lo['F6'], dl=lo['delta'], bb=lo['bars'][0],
                ab=lo['bars'][1], D_auth=lo['D_bound'], parts_auth=lo['parts'], sig6_auth=lo['sig6'], g=g)

if __name__ == "__main__":
    G, d = int(sys.argv[1]), int(sys.argv[2])
    mp.dps = 60
    P = authors_params(G, d)
    mp.dps = 60
    r = evaluate(G, d, P)
    b4 = P['beta']**4
    print(f"SU({G}), d={d}, beta_* = {nstr(P['beta'], 12)}; parameters of the authors: varpi = b^(5/2), (bbar, abar) = ({P['bb']}, {P['ab']})")
    print(f"  group bounds: ubar/beta = {nstr(r['gb']['ubar']/P['beta'], 10)}, ulow/beta = {nstr(r['gb']['ulow']/P['beta'], 10)}, u2bar/beta^2 = "
          f"{nstr(r['gb']['u2']/P['beta']**2, 8)}, Delta/beta^4 = {nstr(r['gb']['Delta']/b4, 6)}, Ubar in [{nstr(r['gb']['Ulo'], 10)}, {nstr(r['gb']['Uhi'], 10)}]")
    print(f"  D/beta^4: mine {nstr(r['D']/b4, 10)}, authors {nstr(P['D_auth']/b4, 10)}")
    print("  parts of Sigma^(4)_4 / beta^4, mine vs authors: " + ', '.join(f"{k} {nstr(v/b4, 8)} / {nstr(P['parts_auth'][k]/b4, 8)}" for k, v in r['parts'].items()))
    print(f"  Sigma^(6)/b: mine {nstr(r['sig6']/r['b'], 12)}, authors {nstr(P['sig6_auth']/r['b'], 12)}")
    print("  rooms (1 - lhs/rhs) of the hypotheses with my constants and extended sums: " + ', '.join(f"{k} {nstr(v, 4)}" for k, v in r['rooms'].items()))
    lab = ('(1) gF', '(2) gC', "(3) g'F", "(4) g'C", '(5) kF', '(6) kC')
    for name, key in (("authors' marked constants, actual, extended sums", 'rooms_auth_a'), ("authors' marked constants, actual, their own truncated sums", 'rooms_auth_a_trunc'),
                      ("authors' marked constants, barred, extended sums", 'rooms_auth_b'), ("authors' marked constants, barred, their own truncated sums", 'rooms_auth_b_trunc'),
                      ("my marked constants, actual", 'rooms_my_a'), ("my marked constants, barred", 'rooms_my_b')):
        print(f"  six inequalities, {name}: " + ', '.join(f"{l} {nstr(v, 4)}" for l, v in zip(lab, r[key])))
    print("  (5), (6) in the form of Lemma mg:lem:marked (G_1 G_2 in place of g g'), the authors' marked constants, extended sums: actual "
          + ', '.join(nstr(v, 4) for v in r['rooms_auth_G'][:2]) + "; barred " + ', '.join(nstr(v, 4) for v in r['rooms_auth_G'][2:]))

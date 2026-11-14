# Referee check (ELEVENTH pass on Section 16, Corollary mg:cor:classic), items 1, 5, 6: the upper bounds of the corollary in BALL ARITHMETIC (Arb, python-flint).
# An adapted copy of the tenth referee's referee-section16j/upper_arb_sectors.py (same tails T'_{S_Q}, T'_{S_A} of Lemma mg:lem:uppertails with sigma5bar, the
# same counts as Theorem mg:thm:thirdcap(d): Z^3 to 12 plaquettes, Z^4 to 10 and the binomial bound at 11). Changes ("16k"):
#   - general (N, iota): SU(2) with N = 2, iota = 1, tau_0 = 1 (the eigenvector chi_rho, a_1 = |<f chi_rho>| = N u = 2u, C = + only), and U(1) with N = 1,
#     iota = 0, tau_0 = 1/2 (sqrt2 Re chi, a_1 = sqrt2 u; sqrt2 Im chi, a_1 = 0). The power-series check before hard-coded iota = 0 ((2-iota) = 2, 1/tau_0^2 = 4).
#   - my own derivation for SU(2) from Theorem mg:thm:upper2 (README section C): with psi = chi_rho, sqrt(2-iota) = 1, so
#       V_2 <= (1+T) Ups + T (1 + 2 a_1 + 3 a_1^2),   log(V_2/Ups) <= T (1 + e^{2b}(1 + 2a + 3a^2))      (Ups >= e^{-2b}, a_1 <= a := N ubar),
#       L_2/Ups^2 >= (1 - T_Q - eps^6)(1 - 2 e^{2b} T (1 + a) N ubar) - (1 + a'^2) e^{20b} N^4 T_Q        (tau_0 = 1),
#     i.e. B_mine = T_Q + eps^6 + 2 e^{2b} N ubar T (1+a) + (1+a'^2) e^{20b} N^4 T_Q, which is sector_upper.py's B with sqrt(2-iota) = 1; Corollary
#     mg:cor:sectorscap's displayed B_C has sqrt2 N ubar (sqrt2 > sqrt(2-iota) for iota = 1: weaker, also valid). All three variants are evaluated.
#   - the evaluation point is an exact binary number above beta_0(d), and the comparison is with the corollary's c_+ (SU(2): 0.0026, 0.0028; U(1): 0.0017,
#     0.0018).
# Usage: python3 upper_arb_classic.py
import math, os, sys
from flint import arb, ctx
from fractions import Fraction as Fr
ctx.prec = 300
ONE = arb(1); E = ONE.exp()
COUNTS = {3: {'A': (2, {6: 8, 9: 10, 10: 47, 11: 60, 12: 96}), 'Q': (6, {6: 12, 9: 30, 10: 93, 11: 144, 12: 162})},
          4: {'A': (4, {6: 32, 9: 60, 10: 506}), 'Q': (18, {6: 60, 8: 24, 9: 262, 10: 1143})}}
CPLUS = {('SU(2)', 3): '0.0026', ('SU(2)', 4): '0.0028', ('U(1)', 3): '0.0017', ('U(1)', 4): '0.0018'}
CASES = (('SU(2)', 2, 1, '+'), ('U(1)', 1, 0, '+-'))

def tails(beta, d, N, iota):
    """T'_{S_Q}, T'_{S_A} of Lemma mg:lem:uppertails with sigma5bar = (2-iota) N^2 ubar^5 + (e^b - 1 - b)^4 (e^b - 1), and eps, ubar"""
    K = 4*(2*d - 3); eps = (2*beta).expm1(); x = E**2*K*eps
    th = x**5/(1 - x)
    w = eps*((K + 1)*th).exp(); q = E*K*w
    ubar = ((1 + iota)*beta + beta**2*beta.exp())/(2*N**2)
    s5 = (2 - iota)*N**2*ubar**5 + (beta.expm1() - beta)**4*beta.expm1()
    T = {}
    for S, nC, PS, Jp in (('Q', 6, 24*d - 42, 24*(d - 2)), ('A', 1, 8*d - 11, K)):
        ncap, Nn = COUNTS[d][S]; Nn = dict(Nn)
        if d == 4: Nn[11] = math.comb(Jp + K*11, 11)
        n0 = max(Nn)
        assert bool(1 - q > 0)
        body = ncap*s5*(5*(K + 1)*th).exp() + sum((arb(c)*w**n for n, c in Nn.items()), arb(0)) + (1 + ONE/K)**Jp*q**(n0 + 1)/(1 - q)
        T[S] = (1 + eps)**nC*(th*PS).exp()*body
    return T['Q'], T['A'], eps, ubar

def E_sector(beta, d, N, iota, C, coef_a=4, sqrt_in_B=None):
    """the corollary's E_C (coef_a = 4); coef_a = 2 is my sharper case analysis; sqrt_in_B = sqrt2 reproduces Corollary mg:cor:sectorscap's displayed B_C"""
    TQ, TA, eps, ubar = tails(beta, d, N, iota)
    tau0 = arb(1 + iota)/2
    s = arb(2 - iota).sqrt() if sqrt_in_B is None else sqrt_in_B
    a = arb(2 - iota).sqrt()*N*ubar if C == '+' else arb(0)
    a2 = (1 + TA)*a + TA
    B = TQ + eps**6 + 2*(2*beta).exp()*s*N*ubar*TA*(1 + a) + (1 + a2**2)*(20*beta).exp()*(arb(N)**4/tau0**2)*TQ
    assert bool(1 - B > 0), 'B_C < 1'
    return TA*(1 + (2*beta).exp()*(1 + coef_a*a + 3*a**2)) + B/(1 - B), B, TQ, TA

# ---------- the power-series structure of E_C/beta^5, exactly (truncated power series in beta with Fraction coefficients) ----------
def ps_exp(c, K): return [Fr(c)**k/math.factorial(k) for k in range(K + 1)]
def ps_mul(A, B, K): return [sum(A[i]*B[k - i] for i in range(k + 1)) for k in range(K + 1)]
def ps_add(*As): return [sum(x) for x in zip(*As)]
def ps_scale(A, c): return [c*x for x in A]
def ps_pow(A, n, K):
    out = [Fr(1)] + [Fr(0)]*K
    for _ in range(n): out = ps_mul(out, A, K)
    return out
def ps_exp_of(A, K):
    out = [Fr(1)] + [Fr(0)]*K; term = [Fr(1)] + [Fr(0)]*K
    for m in range(1, K + 1):
        term = ps_scale(ps_mul(term, A, K), Fr(1, m)); out = ps_add(out, term)
    return out
def ps_inv1m(A, K):
    out = [Fr(1)] + [Fr(0)]*K; p = [Fr(1)] + [Fr(0)]*K
    for _ in range(K):
        p = ps_mul(p, A, K); out = ps_add(out, p)
    return out

def power_series_check(d, N, iota, C, K=14):
    """coefficients of E_C to beta^K, with e -> 2718/1000 and sqrt(2-iota) -> a rational stand-in (they enter only as positive multiplicative constants, so the
    signs of the coefficients do not depend on their values); 16k: general iota ((1+iota) in ubar, (2-iota) in sigma5bar and a, N^4/tau_0^2 = 4N^4/(1+iota)^2)"""
    Kc = 4*(2*d - 3); ee = Fr(2718, 1000); s2 = Fr(1414, 1000) if iota == 0 else Fr(1)
    one = [Fr(1)] + [Fr(0)]*K
    eps = ps_add(ps_exp(2, K), ps_scale(one, -1))
    x = ps_scale(eps, ee**2*Kc)
    th = ps_mul(ps_pow(x, 5, K), ps_inv1m(x, K), K)
    w = ps_mul(eps, ps_exp_of(ps_scale(th, Kc + 1), K), K); q = ps_scale(w, ee*Kc)
    eb = ps_exp(1, K)
    ubar = ps_scale(ps_add([Fr(0), Fr(1 + iota)] + [Fr(0)]*(K - 1), ps_mul([Fr(0), Fr(0), Fr(1)] + [Fr(0)]*(K - 2), eb, K)), Fr(1, 2*N*N))
    em1 = ps_add(eb, ps_scale(one, -1)); em1b = ps_add(em1, [Fr(0), Fr(-1)] + [Fr(0)]*(K - 1))
    s5 = ps_add(ps_scale(ps_pow(ubar, 5, K), (2 - iota)*N*N), ps_mul(ps_pow(em1b, 4, K), em1, K))
    T = {}
    for S, nC, PS, Jp in (('Q', 6, 24*d - 42, 24*(d - 2)), ('A', 1, 8*d - 11, Kc)):
        ncap, Nn = COUNTS[d][S]; Nn = dict(Nn)
        if d == 4: Nn[11] = math.comb(Jp + Kc*11, 11)
        n0 = max(Nn)
        body = ps_add(ps_scale(ps_mul(s5, ps_exp_of(ps_scale(th, 5*(Kc + 1)), K), K), ncap),
                      *[ps_scale(ps_pow(w, n, K), c) for n, c in Nn.items()],
                      ps_scale(ps_mul(ps_pow(q, n0 + 1, K), ps_inv1m(q, K), K), (1 + Fr(1, Kc))**Jp))
        T[S] = ps_mul(ps_mul(ps_pow(ps_add(one, eps), nC, K), ps_exp_of(ps_scale(th, PS), K), K), body, K)
    TQ, TA = T['Q'], T['A']
    a = ps_scale(ubar, s2*N) if C == '+' else [Fr(0)]*(K + 1)
    a2 = ps_add(ps_mul(ps_add(one, TA), a, K), TA)
    e2b, e20b = ps_exp(2, K), ps_exp(20, K)
    B = ps_add(TQ, ps_pow(eps, 6, K), ps_scale(ps_mul(ps_mul(e2b, ubar, K), ps_mul(TA, ps_add(one, a), K), K), 2*s2*N),
               ps_scale(ps_mul(ps_mul(ps_add(one, ps_pow(a2, 2, K)), e20b, K), TQ, K), Fr(4*N**4, (1 + iota)**2)))
    Ec = ps_add(ps_mul(TA, ps_add(one, ps_mul(e2b, ps_add(one, ps_scale(a, 4), ps_scale(ps_pow(a, 2, K), 3)), K)), K), ps_mul(B, ps_inv1m(B, K), K))
    return all(c == 0 for c in Ec[:5]), all(c >= 0 for c in Ec), all(c >= 0 for c in TQ) and all(c >= 0 for c in TA) and all(c >= 0 for c in B), Ec[5]

if __name__ == "__main__":
    ok = True
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    for d in (3, 4):
        K = 4*(2*d - 3)
        b0 = (ONE + ONE/(E**2*K*(K + 2))).log()/2                     # beta_0(d) as a ball
        bup = arb(arb(b0.upper())*(1 + arb('1e-30'))); bup = arb(bup.upper())   # an exact number above beta_0(d)
        assert bool(bup - b0 > 0)
        for name, N, iota, sectors in CASES:
            for C in sectors:
                Ec, B, TQ, TA = E_sector(bup, d, N, iota, C); Em, _, _, _ = E_sector(bup, d, N, iota, C, coef_a=2)
                Es2, _, _, _ = E_sector(bup, d, N, iota, C, sqrt_in_B=arb(2).sqrt())
                v = Ec/bup**4
                good = bool(arb(CPLUS[(name, d)]) - v > 0)
                ok = ok and good
                print(f"{name}, d={d}, C={C}: at beta = beta_0(d)(1+1e-30) (exact, above beta_0): E_C/beta^4 = {v.str(10)} (mine, 2a: {(Em/bup**4).str(10)};"
                      f" with sqrt2 in B as displayed in Cor. 16.20: {(Es2/bup**4).str(10)}), B_C = {B.str(4)}, T'_Q/beta^4 = {(TQ/bup**4).str(4)},"
                      f" T'_A/beta^4 = {(TA/bup**4).str(4)};  < c_+ = {CPLUS[(name, d)]}: {good}")
    print("all certified:", ok)
    # against the authors' sector_upper.E_sector (mpmath, 60 digits) at beta_0(d)
    from mpmath import mp, mpf, binomial
    import sector_upper as SUP, upper_exact as UE
    from constants_upper import beta0
    mp.dps = 60
    print("\nagainst sector_upper.E_sector at beta_0(d) (relative difference; the authors' printed values are E_C(beta_0)/beta_0^4):")
    for d in (3, 4):
        b0m = beta0(d); Kd = 4*(2*d - 3)
        if d == 3: cn = {3: UE.counts_to(3, 12)}
        else:
            cn = {4: {'Q': ({**{n: c for n, c in UE.COUNTS[4]['Q'][0].items() if n <= 10}, 11: binomial(48 + 11*Kd, 11)}, 18),
                      'A': ({**{n: c for n, c in UE.COUNTS[4]['A'][0].items() if n <= 10}, 11: binomial(20 + 11*Kd, 11)}, 4)}}
        for name, N, iota, sectors in CASES:
            for C in sectors:
                theirs = SUP.E_sector(b0m, d, N, iota, C, cn)
                mine, _, _, _ = E_sector(arb(mp.nstr(b0m, 60)), d, N, iota, C)
                rel = abs(mine.mid()/arb(mp.nstr(theirs, 60)) - 1)
                print(f"  {name}, d={d}, C={C}: authors {mp.nstr(theirs/b0m**4, 12)} beta^4, mine {(mine/arb(mp.nstr(b0m, 60))**4).str(12, radius=False)} beta^4, relative difference {rel.str(3, radius=False)}")
    print("\nthe power-series structure of E_C (coefficients of beta^0..beta^14 of the exact expansion; E_C = O(beta^5)):")
    for d in (3, 4):
        for name, N, iota, sectors in CASES:
            for C in sectors:
                low, nonneg, parts, c5 = power_series_check(d, N, iota, C)
                print(f"  {name}, d={d}, C={C}: coefficients of beta^0..4 vanish: {low};  all coefficients to beta^14 >= 0: {nonneg} (T'_Q, T'_A, B_C: {parts});"
                      f"  coefficient of beta^5 = {float(c5):.6g}")

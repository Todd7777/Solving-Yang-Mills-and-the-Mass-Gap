# Referee check (TENTH pass on Section 16, Corollary mg:cor:sectorscap), items 3 and 7: the upper bound for each sector of charge conjugation, in BALL
# ARITHMETIC (Arb, python-flint). An adapted copy of the ninth referee's referee-section16i/upper_arb.py (same tails T'_{S_Q}, T'_{S_A} of Lemma
# mg:lem:uppertails with sigma5bar, the same counts as Theorem mg:thm:thirdcap(d): Z^3 to 12 plaquettes, Z^4 to 10 and the binomial bound at 11), with the
# referee's own derivation of the sector bound from Theorem mg:thm:upper2 (README, section C):
#   psi = sqrt2 Re chi (C=+) or sqrt2 Im chi (C=-), eigenvalue U_C of V; a_1 = |<f psi>| = sqrt2 N u (C=+), 0 (C=-); only U_C >= e^{-2 beta}.
#   V_2 <= (1+T) U_C + T (1 + 2 a_1 + 3 a_1^2)                  (the case analysis of Lemma mg:lem:upperE),
#   so log(V_2/U_C) <= T (1 + e^{2b} (1 + 2a + 3a^2))         (log(1+x) <= x, U_C >= e^{-2b}, a_1 <= a);
#   L_2/U_C^2 >= (1 - T_Q - eps^6)(1 - 2 e^{2b} x) - (1 + a'^2) e^{16b} e^{4b} N^4/tau_0^2 T_Q >= 1 - B_mine,   x = T (1 + a_1) sqrt2 N u <= T (1 + a) sqrt2 N ubar,
#   B_mine := T_Q + eps^6 + 2 e^{2b} sqrt2 N ubar T (1 + a) + (1 + a'^2) e^{20b} N^4/tau_0^2 T_Q   (= the corollary's B_C),
#   E_mine := T (1 + e^{2b}(1 + 2a + 3a^2)) + B/(1 - B);  the corollary's E_C has 4a in place of 2a (weaker, also valid).
# Both E_C (the corollary's) and E_mine are evaluated at an exact binary number above beta_0(d), and compared with c'_+ of the table.
# The structure "E_C/beta^5 is a power series with non-negative coefficients" is checked separately (power_series_check below): every factor is expanded
# as a truncated power series with exact rational coefficients (sympy-free, by Fractions) and all coefficients of T'_S/beta^5, B_C/beta^5, E_C/beta^5 up to
# order 12 are checked to be >= 0.
# Usage: python3 upper_arb_sectors.py
import math
from flint import arb, ctx
ctx.prec = 300
ONE = arb(1); E = ONE.exp()
COUNTS = {3: {'A': (2, {6: 8, 9: 10, 10: 47, 11: 60, 12: 96}), 'Q': (6, {6: 12, 9: 30, 10: 93, 11: 144, 12: 162})},
          4: {'A': (4, {6: 32, 9: 60, 10: 506}), 'Q': (18, {6: 60, 8: 24, 9: 262, 10: 1143})}}
CPLUS = {(3, 3): '0.042', (3, 4): '0.047', (4, 3): '0.131', (4, 4): '0.148', (5, 3): '0.319', (5, 4): '0.361'}

def tails(beta, d, N, iota=0):
    """T'_{S_Q}, T'_{S_A} of Lemma mg:lem:uppertails with sigma5bar (as in the ninth referee's upper_arb.py), and eps, ubar"""
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

def E_sector(beta, d, N, C, coef_a=4, iota=0):
    TQ, TA, eps, ubar = tails(beta, d, N, iota)
    tau0 = arb(1 + iota)/2
    a = arb(2 - iota).sqrt()*N*ubar if C == '+' else arb(0)
    a2 = (1 + TA)*a + TA
    B = TQ + eps**6 + 2*(2*beta).exp()*arb(2 - iota).sqrt()*N*ubar*TA*(1 + a) + (1 + a2**2)*(20*beta).exp()*(arb(N)**4/tau0**2)*TQ
    assert bool(1 - B > 0), 'B_C < 1'
    return TA*(1 + (2*beta).exp()*(1 + coef_a*a + 3*a**2)) + B/(1 - B), B

def E_theorem(beta, d, N, iota=0):
    """the E of Lemma mg:lem:upperE (with Upsilon >= 1), for comparison"""
    TQ, TA, eps, ubar = tails(beta, d, N, iota)
    tau0 = arb(1 + iota)/2
    a = arb(2 - iota).sqrt()*N*ubar; a2 = (1 + TA)*a + TA
    B = TQ + eps**6 + 2*a*(1 + a)*TA + (1 + a2**2)*(16*beta).exp()*(arb(N)**4/tau0**2)*TQ
    return TA*(2 + 3*a + 3*a**2) + B/(1 - B)

# ---------- the power-series structure, exactly (truncated power series in beta with Fraction coefficients) ----------
from fractions import Fraction as Fr
def ps_exp(c, K):                      # e^{c beta}
    return [Fr(c)**k/math.factorial(k) for k in range(K + 1)]
def ps_mul(A, B, K): return [sum(A[i]*B[k - i] for i in range(k + 1)) for k in range(K + 1)]
def ps_add(*As): return [sum(x) for x in zip(*As)]
def ps_scale(A, c): return [c*x for x in A]
def ps_pow(A, n, K):
    out = [Fr(1)] + [Fr(0)]*K
    for _ in range(n): out = ps_mul(out, A, K)
    return out
def ps_exp_of(A, K):                   # e^{A(beta)} for A(0) = 0, via the exponential series
    out = [Fr(1)] + [Fr(0)]*K; term = [Fr(1)] + [Fr(0)]*K
    for m in range(1, K + 1):
        term = ps_scale(ps_mul(term, A, K), Fr(1, m)); out = ps_add(out, term)
    return out
def ps_inv1m(A, K):                    # 1/(1 - A) for A(0) = 0
    out = [Fr(1)] + [Fr(0)]*K; p = [Fr(1)] + [Fr(0)]*K
    for _ in range(K):
        p = ps_mul(p, A, K); out = ps_add(out, p)
    return out

def power_series_check(d, N, C, K=14):
    """the coefficients of E_C, with e and sqrt2 kept as exact rational upper/lower stand-ins is not possible; instead every constant factor is positive
    and every beta-dependence is a series with non-negative coefficients. We check this by expanding in beta with the constants e, sqrt2 replaced by
    rational numbers (e -> 2718/1000, sqrt2 -> 1414/1000): the signs of the coefficients do not depend on the values of these positive constants, since
    they enter only as positive multiplicative constants (e^2 K in x, e K in q, sqrt2 N in a)."""
    Kc = 4*(2*d - 3); ee, s2 = Fr(2718, 1000), Fr(1414, 1000)
    one = [Fr(1)] + [Fr(0)]*K
    eps = ps_add(ps_exp(2, K), ps_scale(one, -1))                          # e^{2b} - 1
    x = ps_scale(eps, ee**2*Kc)
    th = ps_mul(ps_pow(x, 5, K), ps_inv1m(x, K), K)                        # x^5/(1-x)
    w = ps_mul(eps, ps_exp_of(ps_scale(th, Kc + 1), K), K); q = ps_scale(w, ee*Kc)
    eb = ps_exp(1, K)
    ubar = ps_scale(ps_add([Fr(0), Fr(1)] + [Fr(0)]*(K - 1), ps_mul([Fr(0), Fr(0), Fr(1)] + [Fr(0)]*(K - 2), eb, K)), Fr(1, 2*N*N))
    em1 = ps_add(eb, ps_scale(one, -1)); em1b = ps_add(em1, [Fr(0), Fr(-1)] + [Fr(0)]*(K - 1))
    s5 = ps_add(ps_scale(ps_pow(ubar, 5, K), 2*N*N), ps_mul(ps_pow(em1b, 4, K), em1, K))
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
               ps_scale(ps_mul(ps_mul(ps_add(one, ps_pow(a2, 2, K)), e20b, K), TQ, K), Fr(N**4*4)))
    Ec = ps_add(ps_mul(TA, ps_add(one, ps_mul(e2b, ps_add(one, ps_scale(a, 4), ps_scale(ps_pow(a, 2, K), 3)), K)), K), ps_mul(B, ps_inv1m(B, K), K))
    low = all(c == 0 for c in Ec[:5])
    return low, all(c >= 0 for c in Ec), all(c >= 0 for c in TQ) and all(c >= 0 for c in TA) and all(c >= 0 for c in B), Ec[5]

if __name__ == "__main__":
    ok = True
    for d in (3, 4):
        K = 4*(2*d - 3)
        b0 = (ONE + ONE/(E**2*K*(K + 2))).log()/2                     # beta_0(d) as a ball
        bup = arb(arb(b0.upper())*(1 + arb('1e-30'))); bup = arb(bup.upper())   # an exact number above beta_0(d)
        assert bool(bup - b0 > 0)
        for N in (3, 4, 5):
            Eth = E_theorem(bup, d, N)/bup**4
            for C in ('+', '-'):
                Ec, B = E_sector(bup, d, N, C); Em, _ = E_sector(bup, d, N, C, coef_a=2)
                v, vm = Ec/bup**4, Em/bup**4
                good = bool(arb(CPLUS[(N, d)]) - v > 0)
                ok = ok and good
                print(f"SU({N}), d={d}, C={C}: at beta = beta_0(d)(1+1e-30) (exact, above beta_0): E_C/beta^4 = {v.str(10)} (mine, 2a: {vm.str(8)};"
                      f" Lemma upperE with Upsilon >= 1: {Eth.str(8)}), B_C = {B.str(4)};  < c'_+ = {CPLUS[(N, d)]}: {good}")
    print("all certified:", ok)
    print("\nthe power-series structure of E_C (coefficients of beta^0..beta^14 of the exact expansion; E_C = O(beta^5)):")
    for d in (3, 4):
        for N in (3, 4, 5):
            for C in ('+', '-'):
                low, nonneg, parts, c5 = power_series_check(d, N, C)
                print(f"  SU({N}), d={d}, C={C}: coefficients of beta^0..4 vanish: {low};  all coefficients to beta^14 >= 0: {nonneg} (T'_Q, T'_A, B_C: {parts});"
                      f"  coefficient of beta^5 = {float(c5):.6g}")

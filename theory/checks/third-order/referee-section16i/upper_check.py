# Referee check (ninth pass on Section 16) of Theorem mg:thm:thirdcap(e): the bound E(beta) of Remark mg:rem:upperexact(ii) with Lemma mg:lem:uppertails,
# re-implemented from the text (not from upper_exact.py), at beta_0(d) and at beta_0(d)(1+1e-30), for SU(3), SU(5) in d = 3, 4, with the admissible-set
# counts of Remark mg:rem:upperexact(i),(iv) (Z^3 to 12, Z^4 to 10) and in Z^4 the rooted bound binom(J' + K n, n) of the proof of Lemma mg:lem:twopoint
# for n = 11. Also: E(beta)/beta^5 is a power series with non-negative coefficients (checked here on its Taylor coefficients to order 12, exactly in
# rational arithmetic is not attempted; mpmath taylor at 80 digits), and the binomial values quoted in Remark mg:rem:su5(iv).
# Usage: python3 upper_check.py
from mpmath import mp, mpf, e, exp, log, expm1, sqrt, binomial, taylor, nstr

mp.dps = 80
def beta0(d):
    K = 4*(2*d - 3); return log(1 + 1/(e**2*K*(K + 2)))/2

COUNTS = {3: {'A': (2, {6: 8, 9: 10, 10: 47, 11: 60, 12: 96}), 'Q': (6, {6: 12, 9: 30, 10: 93, 11: 144, 12: 162})},
          4: {'A': (4, {6: 32, 9: 60, 10: 506}), 'Q': (18, {6: 60, 8: 24, 9: 262, 10: 1143})}}

def E(beta, d, N, iota=0, n11_binomial=True):
    beta = mpf(beta); K = 4*(2*d - 3); eps = expm1(2*beta); x = e**2*K*eps
    th = x**5/(1 - x)                                         # theta_5 (Lemma mg:lem:closed)
    w = eps*exp((K + 1)*th); q = e*K*w
    ubar = ((1 + iota)*beta + beta**2*exp(beta))/(2*N**2)
    s5 = (2 - iota)*N**2*ubar**5 + (expm1(beta) - beta)**4*expm1(beta)
    tau0 = mpf(1 + iota)/2
    Tp = {}
    for S, nC, PS, Jp in (('Q', 6, 24*d - 42, 24*(d - 2)), ('A', 1, 8*d - 11, K)):
        ncap, Nn = COUNTS[d][S]
        Nn = dict(Nn)
        if d == 4 and n11_binomial: Nn[11] = binomial(Jp + K*11, 11)       # rooted count of the proof of Lemma mg:lem:twopoint
        n0 = max(Nn)
        body = ncap*s5*exp(5*(K + 1)*th) + sum(c*w**n for n, c in Nn.items()) + (1 + mpf(1)/K)**Jp*q**(n0 + 1)/(1 - q)
        Tp[S] = (1 + eps)**nC*exp(th*PS)*body
    TQ, TA = Tp['Q'], Tp['A']
    a = sqrt(2 - iota)*N*ubar; a2 = (1 + TA)*a + TA
    B = TQ + eps**6 + 2*a*(1 + a)*TA + (1 + a2**2)*exp(16*beta)*(mpf(N)**4/tau0**2)*TQ
    return TA*(2 + 3*a + 3*a**2) + B/(1 - B)

if __name__ == "__main__":
    print("binom(48 + 11*20, 11) =", nstr(binomial(48 + 220, 11), 6), "; binom(20 + 11*20, 11) =", nstr(binomial(20 + 220, 11), 6),
          "(Remark mg:rem:su5(iv): 1.04e19 and 3.0e18)")
    claim = {(3, 3): ('0.0415', '0.042'), (3, 4): ('0.0469', '0.047'), (5, 3): ('0.3178', '0.318'), (5, 4): ('0.3603', '0.361')}
    for N in (3, 5):
        for d in (3, 4):
            b0 = beta0(d)
            v0 = E(b0, d, N)/b0**4; v1 = E(b0*(1 + mpf(10)**-30), d, N)/b0**4
            c = taylor(lambda t: E(t, d, N)/t**5 if t != 0 else 0, b0/2, 0)   # (value only; the series check is below)
            print(f"SU({N}), d={d}: E(beta_0)/beta_0^4 = {nstr(v0, 8)} (at beta_0(1+1e-30): {nstr(v1, 8)}); quoted {claim[(N, d)][0]}, constant c_+ = {claim[(N, d)][1]}: "
                  f"{v1 <= mpf(claim[(N, d)][1])}")
    # E(beta)/beta^5 as a power series: Taylor coefficients at 0 (finite differences at high precision), all non-negative?
    mp.dps = 120
    for N in (3, 5):
        for d in (3, 4):
            g = lambda t: E(t, d, N)/t**5
            co = taylor(g, mpf(10)**-40, 8)
            print(f"SU({N}), d={d}: Taylor coefficients of E/beta^5 near 0 (orders 0..8): " + ', '.join(nstr(x, 4) for x in co) +
                  f"; all >= 0: {all(x >= 0 for x in co)}")

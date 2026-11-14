# Part IV, Section 16 (Corollary mg:cor:sectorscap): the upper bound for each sector of charge conjugation on the whole range. Usage: python3 sector_upper.py
"""For the sector C, Theorem mg:thm:upper2 with psi the unit eigenvector of V in the sector C (sqrt2 Re chi for +, sqrt2 Im chi for -), eigenvalue U_C, and
a_1 = |<f psi>| = sqrt(2-iota) N u for C = + and 0 for C = - (f is C-invariant), gives m_C <= 4 log(1/u) + log(V_2/L_2) (the proof of Theorem mg:thm:sectors(ii)).
With only U_C >= e^{-2 beta} (Lemma mg:lem:upsbounds(i)), a >= a_1, a' := (1+T_p) a + T_p >= a_2:
    log(V_2/U_C) <= T_p (1 + e^{2beta} (1 + 4a + 3a^2)),
    -log(L_2/U_C^2) <= B_C/(1 - B_C),  B_C := T_Q + eps^6 + 2 e^{2beta} sqrt(2-iota) N ubar T_p (1 + a) + (1 + a'^2) e^{20 beta} N^4/tau_0^2 T_Q,
so m_C - 4 log(1/u) + log U_C <= E_C := T_p (1 + e^{2beta}(1 + 4a + 3a^2)) + B_C/(1 - B_C), beta^5 times a power series with non-negative coefficients
(as in the proof of Lemma mg:lem:upperE); it is also checked to increase on 401 points. The counts are those of Theorem mg:thm:thirdcap(d) (certify_upper.py)."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
for sub in ('second-order', 'upper-bound'): sys.path.insert(0, os.path.join(here, '..', sub))
from mpmath import mp, mpf, exp, expm1, sqrt, binomial
import upper_exact as UE
from constants_upper import beta0

def E_sector(beta, d, N, iota, C, counts):
    tau0 = mpf(1 + iota)/2
    ubar = ((1 + iota)*beta + beta**2*exp(beta))/(2*N**2)
    TQ, Tp, eps = UE.tails(beta, d, {'N': N, 'iota': iota, 'u': ubar}, counts=counts, kind='theta5', ubar=True)
    a = sqrt(2 - iota)*N*ubar if C == '+' else mpf(0)
    a2 = (1 + Tp)*a + Tp
    B = TQ + eps**6 + 2*exp(2*beta)*sqrt(2 - iota)*N*ubar*Tp*(1 + a) + (1 + a2**2)*exp(20*beta)*(mpf(N)**4/tau0**2)*TQ
    return Tp*(1 + exp(2*beta)*(1 + 4*a + 3*a**2)) + B/(1 - B)

if __name__ == "__main__":
    mp.dps = 60
    # the third-order terms of the two glueballs (Corollary mg:cor:sectorscap(i)); none for SU(4), SU(2), U(1), whose -log U are even in beta
    TERMS = {(3, '+'): mpf(1)/432, (3, '-'): mpf(1)/1296, (5, '+'): mpf(1)/6000, (5, '-'): mpf(1)/6000}
    CASES = [('SU(3)', 3, 0, '+-'), ('SU(4)', 4, 0, '+-'), ('SU(5)', 5, 0, '+-'),     # Corollary mg:cor:sectorscap
             ('SU(2)', 2, 1, '+'), ('U(1)', 1, 0, '+-')]                            # Corollary mg:cor:classic (iota = 1: the one plaquette state, C = +)
    if len(sys.argv) > 1: CASES = [c for c in CASES if c[0] in sys.argv[1:]]
    for d in (3, 4):
        b0 = beta0(d); K = 4*(2*d - 3)
        if d == 3: cn = {3: UE.counts_to(3, 12)}
        else:
            cn = {4: {'Q': ({**{n: c for n, c in UE.COUNTS[4]['Q'][0].items() if n <= 10}, 11: binomial(48 + 11*K, 11)}, 18),
                      'A': ({**{n: c for n, c in UE.COUNTS[4]['A'][0].items() if n <= 10}, 11: binomial(20 + 11*K, 11)}, 4)}}
        for name, N, iota, sectors in CASES:
            for C in sectors:
                E = E_sector(b0, d, N, iota, C, cn)
                pts = sorted([mpf(10)**(-30 + j*(30 + mp.log10(b0))/400) for j in range(401)]); pts[-1] = b0
                vals = [E_sector(p, d, N, iota, C, cn)/p**5 for p in pts]
                mono = all(x <= y for x, y in zip(vals, vals[1:]))
                cmp = f"(term/beta^4 >= {float(TERMS[(N, C)]/b0):.4f})" if (N, C) in TERMS else "(no third-order term)"
                lab = f"{name} d={d} sector {C}" if sectors != '+' else f"{name} d={d}"
                print(f"{lab}: E_C(beta_0)/beta^4 = {float(E/b0**4):.5f}  {cmp};  E/beta^5 increasing: {mono}")

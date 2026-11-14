# Part IV, Section 16 (Theorem mg:thm:thirdcap): the upper halves, certified. Usage: python3 certify_upper.py
"""The bound E(beta) of Remark mg:rem:upperexact(ii) for SU(3) and SU(5) in d = 3, 4, with only counts made by two independent programs: the admissible sets in
Z^3 to 12 plaquettes and in Z^4 to 10 (rel_closed_count.py and the fifth referee's two methods), and in Z^4 the bound binom(J' + K n, n) of the proof of
Lemma mg:lem:twopoint for n = 11 (no enumeration). E uses ubar >= u and Upsilon >= 1: Remark mg:rem:upperexact(ii) for SU(3), Lemma mg:lem:upsbounds(ii) for SU(5)
(rigorous_bounds.prove_upsilon_ge_1). E is beta^5 times a function that visibly increases (Remark mg:rem:upperexact(ii)); it is also checked on 401 points."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from mpmath import mp, mpf, binomial
import upper_exact as UE
import third_order_effective as T
from sun_group import sun
from rigorous_bounds import prove_upsilon_ge_1
from constants_upper import beta0

if __name__ == "__main__":
    mp.dps = 60
    r = prove_upsilon_ge_1(5, mpf('4.03e-4'))
    print(f"SU(5): Upsilon >= U_+ >= 1 on (0, 4.03e-4], which contains (0, beta_0(3)]: {r['ok']} (leading coefficient {r['lead']} of t^{r['k0']}, room {float(r['margin']):.6g})")
    TERM = {3: mpf(1)/432, 5: mpf(1)/6000}
    for N, gfun in ((3, dict(T.GROUPS)['SU(3)']), (5, sun(5))):
        for d in (3, 4):
            b0 = beta0(d); K = 4*(2*d - 3)
            if d == 3:
                C = {3: UE.counts_to(3, 12)}
                how = "counts to 12"
            else:
                C = {4: {'Q': ({**{n: c for n, c in UE.COUNTS[4]['Q'][0].items() if n <= 10}, 11: binomial(48 + 11*K, 11)}, 18),
                         'A': ({**{n: c for n, c in UE.COUNTS[4]['A'][0].items() if n <= 10}, 11: binomial(20 + 11*K, 11)}, 4)}}
                how = "counts to 10, binomial at 11"
            E = UE.E2(b0, d, gfun(b0), counts=C)
            pts = sorted([mpf(10)**(-30 + j*(30 + mp.log10(b0))/400) for j in range(401)]); pts[-1] = b0
            vals = [UE.E2(p, d, gfun(p), counts=C)/p**5 for p in pts]
            mono = all(x <= y for x, y in zip(vals, vals[1:]))
            print(f"SU({N}), d={d} ({how}): E(beta_0) = {float(E/b0**4):.6g} beta^4 against {float(TERM[N]/b0):.5g} beta^4: {E/b0**4 < TERM[N]/b0};"
                  f" E/beta^5 increasing on 401 points: {mono}")

# Referee (Section 16h), after the authors' count of the admissible sets with 11 plaquettes in Z^4 (orders/admissible_d4_n11.out: 760 for S_A, 2616 for S_Q):
# the upper-bound error E(beta_0(4)) of Remark mg:rem:upperexact(ii) by the referee's evaluator (upper_su5.py, formulas of the paper), for SU(5) and SU(3),
# with these counts, and E/beta^5 on 401 points.  Usage: python3 upper_final.py
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from mpmath import mp, mpf
import upper_su5 as U
mp.dps = 50
d = 4; b0 = U.beta0(d)
cnt = U.with11(U.COUNTS[4], 2616, 760)
pts = sorted([mpf(10)**(-30 + j*(30 + mp.log10(b0))/400) for j in range(401)]); pts[-1] = b0
for N, term, name in ((5, mpf(1)/6000, 'SU(5)'), (3, mpf(1)/432, 'SU(3)')):
    U.N = N; U.iota = 0; U.tau0 = mpf(1)/2
    for label, c in (('counts to 10', U.COUNTS[4]), ('counts to 11 (760, 2616)', cnt)):
        Ev = U.E(b0, d, c)
        vals = [U.E(p, d, c)/p**5 for p in pts]
        print(f"{name}, d=4, {label}: E(beta_0(4)) = {mp.nstr(Ev/b0**4, 6)} beta^4 against beta^3*{mp.nstr(term, 4)} = {mp.nstr(term/b0, 5)} beta^4;"
              f" E/beta^5 increasing on 401 points: {all(x <= y for x, y in zip(vals, vals[1:]))}")

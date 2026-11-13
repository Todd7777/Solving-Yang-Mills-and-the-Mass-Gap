# Section 16 of Part IV (Remark mg:rem:nextorder (viii)): the upper bound and the third-order term.  Usage: python3 upper_third_order.py
"""The error of the upper bound of Theorem mg:thm:upper2 for SU(3), D_2 + log(Ups) (D_2 = log(V_2/L_2), second_order_constants.upper2), against the
third-order term beta^3/432: the largest beta (bisection on a log scale) below which the upper bound places m - 4log(1/u) + log(Ups) below beta^3/432,
in d = 3, 4, and the size of the tails T_Q, T_p of Lemma mg:lem:twopoint at beta_0(3).  The error is O(beta^5) with a constant of order 10^12 to 10^13,
almost all from T_Q (crude rooted counts (1 + 1/K)^J' (eK)^n of the sets of Lemma mg:lem:five, activity eps^n)."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
for sub in ('second-order', 'upper-bound'): sys.path.insert(0, os.path.join(here, '..', sub))
sys.path.insert(0, here)
from mpmath import mpf, log, sqrt
import third_order_effective as T
from second_order_constants import upper2
from constants_upper import consts, beta0

su3 = dict(T.GROUPS)['SU(3)']

def err(beta, d):
    g = su3(beta); up = upper2(beta, d, g)
    return None if up['D2'] is None else up['D2'] + log(g['Ups'])

if __name__ == "__main__":
    for d in (3, 4):
        lo, hi = mpf('1e-12'), mpf('1e-4')
        for _ in range(80):
            mid = sqrt(lo*hi); e = err(mid, d)
            if e is not None and e < mid**3/432: lo = mid
            else: hi = mid
        print(f"d={d}: the upper bound places m - 4log(1/u) + log Ups below beta^3/432 for beta <= {float(lo):.3g}; error there {float(err(lo, d)/lo**5):.3g} beta^5")
        for be in ('1e-6', '1e-5', '1e-4'):
            beta = mpf(be); e = err(beta, d)
            print(f"   beta = {be}: error/beta^4 = {float(e/beta**4):.3g} against beta^3/(432 beta^4) = {float(1/(432*beta)):.4g}")
    b0 = beta0(3); c = consts(b0, 3, 3, 0)
    print(f"beta_0(3) = {float(b0):.6g}: T_Q = {float(c['TQ']):.3g} = {float(c['TQ']/b0**4):.3g} beta^4,  T_p = {float(c['Tp']):.3g} = {float(c['Tp']/b0**4):.3g} beta^4;"
          f"  error/beta^4 = {float(err(b0, 3)/b0**4):.3g}")

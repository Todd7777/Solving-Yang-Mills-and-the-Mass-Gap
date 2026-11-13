# Referee (Section 16, fourth pass): how much D at beta_0(3) depends on where the closed polymers are counted crudely (d = 3, S''_9).
# Usage: python3 theta_sensitivity.py
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import interp_frozen as I
from mpmath import mp, log
mp.dps = 60
b0 = I.beta0(3)
for cf, ce in ((15, {10: 10, 11: 11, 12: 24, 13: 0, 14: 70}), (14, {10: 10, 11: 11, 12: 24, 13: 0}), (13, {10: 10, 11: 11, 12: 24})):
    I.T.CLOSED_FROM = cf; I.T.CLOSED_EXACT = {3: ce}
    D, best, g = I.PB.deficit(3, b0)
    th = best['loss'] - log(best['Ubar']) - log(1 + best['delta'])
    print(f"closed polymers exact below {cf}, crude from {cf}: D/beta^4 = {float(D/b0**4):.4f} (8(K+1)theta {float(th/b0**4):.4f}) against {float(1/(432*b0)):.4f}")

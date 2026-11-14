# Referee sanity check (eleventh pass, Corollary mg:cor:classic; Numerical, not part of a proof): the deficit D at beta_* with the EXACT group data of SU(2)
# and U(1) (second_order_constants.su2 / u1: u, Upsilon, Ubar from Bessel functions; D = 8(K+1)theta + log(Ubar/Ups) + log(1+delta')) against bound mode
# (the corollary's D_bound). Bound mode must give the larger D; the difference shows what the proven bounds cost. Same sets, counts, KMAX, F_RAISE as run_cap.sh.
# Usage: python3 exact_vs_bound.py GROUP d
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); parent = os.path.join(here, '..'); sys.path.insert(0, parent)
from mpmath import mpf, nstr
import practical_beta as PB, third_order_effective as T
G, d = int(sys.argv[1]), int(sys.argv[2])
B = {3: '0.00040262143352281974401285396895067832618294529962141', 4: '0.00015376644804254827379172376907999414111262406766301'}[d]
PB.GROUP = G; PB.CERTIFIED = True
if d == 4: T.KMAX_DEFAULT = 100
if G == 1 and d == 3: T.F_RAISE = mpf(10)**-7
if d == 3: PB.load_terms(3, os.path.join(parent, 'orders', 'd3_s9z.json'), {10: 10, 11: 11, 12: 24, 13: 0, 14: 70}, 15)
else: PB.load_terms(4, os.path.join(parent, 'orders', 'd4_s10x_92.json'), {10: 60, 11: 66, 12: 320, 13: 0, 14: 1316, 15: 2520, 16: 12552}, 17)
beta = mpf(B); out = {}
for mode in (False, True):
    PB.BOUNDS = mode
    D, best, g = PB.deficit(d, beta)
    out[mode] = D/beta**4
    extra = f", Ups = {nstr(g['Ups'], 12)}, Ubar = {nstr(best['Ubar'], 12)}" if not mode else f", Delta/beta^4 = {nstr(g['Delta']/beta**4, 4)}"
    print(f"{'U(1)' if G == 1 else 'SU(2)'}, d={d}, {'bound mode' if mode else 'exact group data'}: D(beta_*)/beta_*^4 = {nstr(D/beta**4, 10)}, delta'/beta^4 = {nstr(best['delta']/beta**4, 10)}{extra}")
print(f"  bound mode >= exact: {out[True] >= out[False]};  relative cost of the bounds: {nstr(out[True]/out[False] - 1, 4)}")

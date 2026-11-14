# Referee probe (eleventh pass, Corollary mg:cor:classic): how much of the large lower constants c_- in d = 3 comes from the crude front count
# N^v_k <= (d-1) binom((4d-6)k, k-1) (~14.9^k in Z^2), which the d = 3 setting uses beyond the exact counts (k > 12), although Lemma mg:lem:frontsdecision proves
# N^v_k <= 0.783 * 4.081^k in Z^2 (D = 2) as well. NOT part of any proof: the authors' own pipeline (third_order_effective.lower3e in bound mode, certified
# counts, the full (bbar, abar) grid, varpi = b^2 or b^(5/2)) is rerun with NV_DECISION = True in d = 3, and D(beta_*)/beta_*^4 and the parts of delta' are
# printed, against the setting of the corollary (NV_DECISION = False in d = 3).
# Usage: python3 nv_decision_d3_probe.py GROUP [on|off]     (GROUP = 2 for SU(2), 1 for U(1), 3 for SU(3) for comparison)
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); parent = os.path.join(here, '..')
sys.path.insert(0, parent)
from mpmath import mp, mpf, log
import practical_beta as PB
import third_order_effective as T

B3 = '0.00040262143352281974401285396895067832618294529962141'
G = int(sys.argv[1]); mode = sys.argv[2] if len(sys.argv) > 2 else 'on'
PB.GROUP = G; PB.BOUNDS = True; PB.CERTIFIED = True
if G == 1: T.F_RAISE = mpf(10)**-7
PB.load_terms(3, os.path.join(parent, 'orders', 'd3_s9z.json'), {10: 10, 11: 11, 12: 24, 13: 0, 14: 70}, 15)
T.NV_DECISION = (mode == 'on')            # load_terms sets it False in d = 3; Lemma mg:lem:frontsdecision holds for D = 2
beta = mpf(B3); g = PB.group_data()(beta)
best = None
for vp in (2, mpf(5)/2):
    try: lo = T.lower3e(beta, 3, g, varpi_pow=vp, chars=True)
    except AssertionError: continue
    if best is None or lo['delta'] < best['delta']: best = lo; best['p'] = vp
k = beta
Ds = best['D_bound']/beta**4
parts = ', '.join(f"{kk} {float(v/beta**4):.4g}" for kk, v in best['parts'].items())
name = {1: 'U(1)', 2: 'SU(2)', 3: 'SU(3)'}[G]
print(f"{name}, d=3, NV_DECISION = {T.NV_DECISION}: D(beta_*)/beta_*^4 = {float(Ds):.6g}, (1+k)^3 D/beta_*^4 = {float((1 + k)**3*Ds):.6g};  varpi = b^{float(best['p']):g},"
      f" (bbar, abar) = ({float(best['bars'][0]):.6g}, {float(best['bars'][1]):.6g});  parts of delta'/beta^4: {parts};  stronger-form rooms: "
      + ', '.join(f"{kk} {float(v):.3g}" for kk, v in best['strong'].items()))

# Referee check: sensitivity of the authors' evaluation (third_order_effective.lower3e with fixed parameters, bound mode, certified counts) at beta_*
# to the working precision: 50 vs 100 digits. The parameters are those found at beta_* (as in scaling_check.params), rounded to 50 digits for both runs.
# Usage: python3 precision_check.py GROUP d
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..'))
from mpmath import mp, mpf, nstr
import my_eval as ME
import practical_beta as PB, third_order_effective as T
G, d = int(sys.argv[1]), int(sys.argv[2])
mp.dps = 60
P = ME.authors_params(G, d)
out = {}
for dps in (50, 100):
    mp.dps = 50
    fixed = {k: +mpf(P[k2]) for k, k2 in (('varpi', 'varpi'), ('xF', 'xF'), ('xC', 'xC'), ('hJ', 'hJ'), ('F', 'F'), ('F6', 'F6'), ('dl', 'dl'))}
    beta = +mpf(P['beta']); bb, ab = +mpf(P['bb']), +mpf(P['ab'])
    mp.dps = dps
    T._DEC.clear()
    g = PB.group_data()(beta)
    lo = T.lower3e(beta, d, g, chars=True, fixed=fixed, bbar_grid=[bb], abar_grid=[ab])
    out[dps] = lo
mp.dps = 100
a, b = out[50], out[100]
rel = lambda x, y: abs(x/y - 1) if y else abs(x)
print(f"SU({G}), d={d}, authors' lower3e with fixed parameters: |D_50/D_100 - 1| = {nstr(rel(a['D_bound'], b['D_bound']), 3)};  parts: "
      + ', '.join(f"{k} {nstr(rel(a['parts'][k], b['parts'][k]), 3)}" for k in b['parts']) + f"; Sigma^(6) {nstr(rel(a['sig6'], b['sig6']), 3)}")

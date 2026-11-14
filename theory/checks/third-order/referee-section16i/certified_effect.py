# Referee check: the effect of --certified in d=3 (exits with n = 4, 5 and k >= 17 made crude) on D(beta_*) in bound mode, SU(3) and SU(5).
# Usage: python3 certified_effect.py
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..'))
from mpmath import mp, mpf, nstr
import my_eval as ME
import practical_beta as PB, third_order_effective as T
mp.dps = 60
for G in (3, 5):
    res = {}
    for cert in (False, True):
        PB.GROUP = G; PB.BOUNDS = True; PB.CERTIFIED = cert
        terms, closed, cf, bbs, bstar = ME.SETTINGS[3]
        PB.load_terms(3, os.path.join(here, '..', 'orders', terms), closed, cf)
        beta = mpf(bstar); g = PB.group_data()(beta)
        lo = T.lower3e(beta, 3, g, varpi_pow=mpf(5)/2, chars=True, bbar_grid=list(bbs), abar_grid=[ME.ABAR[(G, 3)]])
        res[cert] = (lo['D_bound']/beta**4, lo['parts']['exits']/beta**4)
    print(f"SU({G}), d=3, beta_*: D/beta^4 = {nstr(res[False][0], 12)} (all exact exits) vs {nstr(res[True][0], 12)} (certified); exits part {nstr(res[False][1], 10)} vs {nstr(res[True][1], 10)}")

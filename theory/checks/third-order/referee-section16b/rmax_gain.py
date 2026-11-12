#!/usr/bin/env python3
"""Effect of replacing the sum R = (b/bb)^4(a/ab)^2 + (b/bb)^3(a/ab)^3 + (b/bb)^2(a/ab)^4 of Lemma mg:lem:jsums5 (and the two-term sum of M_4) by the
maximum of its terms, which is also a valid factor (every remaining pair (n,E) dominates one of the corners (4,2), (3,3), (2,4) componentwise, and
b <= bb, a <= ab).  Settings as in ../practical_beta.py (SU(3), theta_7, exact counts, M_5, Theorem 14.2 by characters).
Usage: python3 rmax_gain.py d beta [maxonly]"""
import sys
import toe_rmax as T
from mpmath import mpf
T.THETA7 = True
T.NOVR = {4: {6: 132, 7: 108, 8: 2034, 9: 4344}, 3: {6: 12, 7: 12, 8: 70, 9: 160, 10: 634}}
T.ZSMALL = {4: {5: 26, 6: 302, 7: 288}, 3: {5: 0, 6: 0, 7: 0}}
T.LOW = {4: (876, 320, 312), 3: (48, 24, 24)}
d, be = int(sys.argv[1]), sys.argv[2]
beta = mpf(be); g = dict(T.GROUPS)['SU(3)'](beta)
for rmax in ((True,) if len(sys.argv) > 3 else (False, True)):
    T.RMAX = rmax
    best = None
    for vp in (2, mpf(5)/2):
        try: lo = T.lower3e(beta, d, g, varpi_pow=vp, chars=True)
        except AssertionError: continue
        if best is None or lo['delta'] < best['delta']: best = lo
    print(f"d={d} beta={be} {'max' if rmax else 'sum'}: delta'/beta^4 = {float(best['delta']/beta**4):.4g}  (beta^3/432 = {float(1/(432*beta)):.4g} beta^4)", flush=True)

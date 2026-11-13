# Referee (Section 16, third pass): re-evaluate the deficit D of Remark mg:rem:practical(iii) / mg:rem:nextorder(iii) with the
# authors' evaluator (../practical_beta.py -> ../third_order_effective.py), and report in addition
#   * which majorant was selected (M_4, M_5, the corner majorant M_S, or Section 11's M_3) and at which (bbar, abar);
#   * whether the barred series are inside their domain of convergence at that point (the evaluator truncates them at
#     NMAX = 40 cluster sizes and KMAX = 60 front sizes, so a divergent series would be silently truncated);
#   * the parts of delta' and the theta-loss.
# Usage: python3 eval_thresholds.py MODE d beta [beta ...]   MODE in {s5, next, next7}
#   s5    : default of practical_beta.py (Lemma mg:lem:jsums5, theta_7)
#   next  : --next (Lemma mg:lem:corners with S_6, theta_10)
#   next7 : --next-theta7 (S_6 with theta_7)
#   noM5  : theta_7 and the exact counts of Remark (ii), without M_5 (the 'alone' thresholds of Remark (iii))
#   noM5sum: as noM5, with the sum of the two corner factors in M_4 (Lemma mg:lem:jsums4 as first stated)
#   next2 : --next2 (d=3 only: S_7', theta with the closed sets of 10, 11, 12 plaquettes counted exactly)
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
import practical_beta as PB
import third_order_effective as T
from mpmath import mpf, log, exp, sqrt

MODE = sys.argv[1]; d = int(sys.argv[2]); betas = sys.argv[3:]
if MODE == 'next':
    PB.ORDER = 6; T.CLOSED_FROM = 10
elif MODE == 'next7':
    PB.ORDER = 6
elif MODE == 'next2':          # the order after S_6 in d=3 (S_7', closed sets of 10-12 plaquettes exact), added by the authors on 13 Nov
    PB.ORDER = 7
elif MODE == 'noM5sum':        # as noM5, and the two corner factors of M_4 summed (no part of Lemma mg:lem:jsums5 at all)
    PB.ORDER = 5; T.LOW = None; T.RMAX = False
elif MODE == 'noM5':           # theta_7 and the exact counts, without Lemma mg:lem:jsums5 (M_4 with the maximum, or M_3)
    PB.ORDER = 5; T.LOW = None
else:
    PB.ORDER = 5
if PB.ORDER == 6: PB.order6(d)
if PB.ORDER == 7: assert d == 3; PB.order7_d3()

REC = {}
_orig_tc = T.tree_constants
def tc_rec(bc, a, cT, Pv, KS, N, KMAX, NMAX):
    out = _orig_tc(bc, a, cT, Pv, KS, N, KMAX, NMAX)
    REC[(bc, a)] = (out, cT, KS)
    return out
T.tree_constants = tc_rec

def conv_report(bars, d):
    """ratios to the radii of convergence of the barred series (crude N_k ~ c^k for large k)"""
    bb, ab = bars[0], bars[1]
    if ab == 'a':            # Section 11's majorant: abar = actual a; look up any record with this bbar
        keys = [k for k in REC if k[0] == bb and REC[k][0] is not None]
        key = min(keys, key=lambda k: k[1])
    else:
        key = (bb, ab)
    (xFb, xCb, Yb, Ab), cT, KS = REC[key]
    J = 4*d - 6
    rY = mpf(KS-1)**(KS-1)/mpf(KS)**KS            # binom(K_S m, m-1) ~ rY^{-m}
    rA = mpf(J-1)**(J-1)/mpf(J)**J                # binom(J k, k-1) ~ rA^{-k}
    abar = key[1]
    return dict(bbar=float(bb), abar=float(abar), xiF_bar=float(xFb), xiC_bar=float(xCb),
                Y_ratio=float(bb*exp(xCb)/rY), A_ratio=float(sqrt(abar)*exp(xFb)/rA))

for be in betas:
    beta = mpf(be); g = dict(T.GROUPS)['SU(3)'](beta); best = None
    for vp in (2, mpf(5)/2):
        REC.clear()
        try:
            lo = T.lower3e(beta, d, g, varpi_pow=vp, chars=True)
        except AssertionError as ex:
            print(f"  varpi=b^{float(vp)}: hypotheses fail ({ex})"); continue
        lo['vp'] = vp; lo['conv'] = conv_report(lo['bars'], d)
        if best is None or lo['delta'] < best['delta']: best = lo
    D = best['loss'] - log(g['Ups'])
    lth = best['loss'] - log(best['Ubar']) - log(1 + best['delta'])
    target = beta**3/432
    print(f"MODE={MODE} d={d} beta={be}: varpi=b^{float(best['vp'])}; D/beta^4 = {float(D/beta**4):.5g}"
          f" [theta-loss {float(lth/beta**4):.4g}, log(Ubar/Ups) {float(log(best['Ubar']/g['Ups'])/beta**4):.3g},"
          f" log(1+delta') {float(log(1+best['delta'])/beta**4):.5g}];  beta^3/432 = {float(target/beta**4):.5g} beta^4"
          f" -> {'resolved' if D < target else 'NOT resolved'}")
    print("   parts/beta^4: " + ', '.join(f"{k} {float(v/beta**4):.4g}" for k, v in best['parts'].items()))
    print(f"   majorant used: {best['bars']}; convergence at the selected point: {best['conv']}")
    print(f"   a = {float(best['a']):.6g}, b = {float(best['b']):.6g}, F = {float(best['F']):.4g}, F6 = {float(best['F6']):.4g}")
    sys.stdout.flush()

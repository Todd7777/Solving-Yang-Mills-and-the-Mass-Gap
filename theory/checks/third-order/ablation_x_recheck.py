# Section 16 of Part IV (Remark mg:rem:nextorder (x)): a recheck of the seventh referee's ablation at beta_0(4), with variants written for it.
# Usage: python3 ablation_x_recheck.py [VARIANT,VARIANT,... [JOBS]]      (a few minutes per variant)
"""D/beta^4 of the lower bound at beta_0(4) for SU(3), d=4, in the setting of (x) (orders/d4_s10x.json, closed polymers of 10..16 plaquettes exact),
with one ingredient removed, made crude or inflated (the threshold is 1/(432 beta_0(4)) = 15.054):
  baseline              the setting of (x)                                   min_set               S(8) u {(6,3),(7,2),(8,2)} (referee-section16g/s8_min_set.json,
  binomial_beyond_13    (d-1) binom((4d-6)k, k-1) in place of Lemma 16.10                          rebuilt by orders/assemble_set.py --base 8 --extra 6:3,7:2,8:2)
  decision_from_12      N^v_12, N^v_13 not used (Lemma 16.10 from 12 edges)  Z4_crude_10_12        Z_(10..12,4) crude
  exits4_crude          the exits with four plaquettes crude                 theta_crude_from_N    theta with the closed polymers exact below N plaquettes only (N = 12, 13, 14)
  orders9_10_x1e9       every pair count with n+E > 17/2 times 10^9          nE_ge_11_2_x1e6       every pair count with n+E >= 11/2 times 10^6
  all_pairs_x1e3        every pair count times 10^3"""
import sys, os
from fractions import Fraction
from multiprocessing import Pool
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
B0 = '1.53766448e-4'
CL = {10: 60, 11: 66, 12: 320, 13: 0, 14: 1316, 15: 2520, 16: 12552}
ALL = ['baseline', 'min_set', 'binomial_beyond_13', 'decision_from_12', 'Z4_crude_10_12', 'exits4_crude', 'theta_crude_from_12', 'theta_crude_from_13',
       'theta_crude_from_14', 'orders9_10_x1e9', 'nE_ge_11_2_x1e6', 'all_pairs_x1e3']
def run(variant):
    os.chdir(here)
    import practical_beta as PB, third_order_effective as T
    from mpmath import mpf
    terms, closed, cf = 'orders/d4_s10x.json', dict(CL), 17
    if variant == 'min_set': terms = 'referee-section16g/s8_min_set.json'
    if variant.startswith('theta_crude_from_'):
        cf = int(variant.split('_')[-1]); closed = {n: c for n, c in CL.items() if n < cf}
    PB.load_terms(4, terms, closed, cf)
    if variant == 'binomial_beyond_13': T.NV_DECISION = False
    if variant == 'decision_from_12': T.NOVR[4] = {k: v for k, v in T.NOVR[4].items() if k <= 11}
    if variant == 'Z4_crude_10_12': T.ZSMALL[4] = {n: c for n, c in T.ZSMALL[4].items() if n <= 9}
    if variant == 'exits4_crude': T.ZEXITN[4] = {k: v for k, v in T.ZEXITN[4].items() if k != 4}
    def scale(cond, f):
        T.EXACT[4] = [((t[0]*f,) + tuple(t[1:])) if cond(t[1], Fraction(t[2]).limit_denominator(2)) else t for t in T.EXACT[4]]
    if variant == 'orders9_10_x1e9': scale(lambda n, E: n + E > Fraction(17, 2), 10**9)
    if variant == 'nE_ge_11_2_x1e6': scale(lambda n, E: n + E >= Fraction(11, 2), 10**6)
    if variant == 'all_pairs_x1e3': scale(lambda n, E: True, 10**3)
    b = mpf(B0); D, best, g = PB.deficit(4, b)
    parts = sorted(((float(v/b**4), k) for k, v in best['parts'].items()), reverse=True)
    return variant, float(D/b**4), parts
if __name__ == '__main__':
    vs = sys.argv[1].split(',') if len(sys.argv) > 1 else ALL
    jobs = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    with Pool(jobs) as p:
        for v, r, parts in p.imap(run, vs):
            print(f"{v}: D/beta^4 = {r:.5g} ({'below' if r < 15.054 else 'ABOVE'} 15.054);  parts: " + ', '.join(f'{k} {a:.3g}' for a, k in parts[:5])); sys.stdout.flush()

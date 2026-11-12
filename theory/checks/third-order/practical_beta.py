# Section 16 of Part IV (Remark mg:rem:practical (iii)): the third-order lower bound at practical beta.  Usage: python3 practical_beta.py [d [beta ...]]
"""The smallest admissible delta' of Theorem mg:thm:lower5 (by characters, varpi = b^2 or b^(5/2)) for SU(3), with
  - theta_7 in place of theta_5 (Lemma mg:lem:theta7);
  - the exact counts of Remark mg:rem:practical (ii): N^v_k (fronts_count.py; d=4: k <= 8), Z_{n,4} for n = 5, 6, 7 (junction_count.py);
  - the junctions with three front-pairs of lowest order counted exactly, M_5 of Lemma mg:lem:jsums5 (N_3, N_2, N'_2 of Lemma mg:lem:threelow,
    three_front_count.py).
It prints delta'/beta^4 against the third-order term beta^3/432 = beta^4/(432 beta) of SU(3), and the three largest parts of delta'/beta^4.
Each beta takes a few minutes (the tree constants on the grid of (bbar, abar))."""
import sys
import third_order_effective as T
from mpmath import mpf
T.THETA7 = True
T.NOVR = {4: {6: 132, 7: 108, 8: 2034}, 3: {6: 12, 7: 12, 8: 70, 9: 160}}
T.ZSMALL = {4: {5: 26, 6: 302, 7: 288}, 3: {5: 0, 6: 0, 7: 0}}
T.LOW = {4: (876, 320, 312), 3: (48, 24, 24)}
DEFAULT = {4: ('1e-5', '3e-6', '1.4e-6', '1e-6'), 3: ('1e-4', '3e-5', '1.7e-5', '1e-5')}
def run(d, betas):
    for be in betas:
        beta = mpf(be); g = dict(T.GROUPS)['SU(3)'](beta); best = None
        for vp in (2, mpf(5)/2):
            try: lo = T.lower3e(beta, d, g, varpi_pow=vp, chars=True)
            except AssertionError: continue
            if best is None or lo['delta'] < best['delta']: best = lo
        target = 1/(432*beta)
        parts = sorted(((float(v/beta**4), k) for k, v in best['parts'].items()), reverse=True)[:3]
        print(f"d={d} beta={be}: delta'/beta^4 = {float(best['delta']/beta**4):.3g}, beta^3/432 = {float(target):.3g} beta^4 "
              f"({'resolved' if best['delta']/beta**4 < target else 'not resolved'});  largest parts: " + ', '.join(f'{k} {a:.3g}' for a, k in parts))
        sys.stdout.flush()
if __name__ == "__main__":
    ds = [int(sys.argv[1])] if len(sys.argv) > 1 else [4, 3]
    for d in ds: run(d, sys.argv[2:] if len(sys.argv) > 2 else DEFAULT[d])

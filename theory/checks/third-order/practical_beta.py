# Section 16 of Part IV (Remark mg:rem:practical (iii)): the third-order lower bound at practical beta.  Usage: python3 practical_beta.py [d [beta ...]]
"""The smallest admissible delta' of Theorem mg:thm:lower5 (by characters, varpi = b^2 or b^(5/2)) for SU(3), with
  - theta_7 in place of theta_5 (Lemma mg:lem:theta7);
  - the exact counts of Remark mg:rem:practical (ii): N^v_k (fronts_count.py, d=3: k <= 9, d=4: k <= 8; the referee's ref_fronts.py: N^v_10 = 634
    in d=3, N^v_9 = 4344 in d=4), Z_{n,4} for n = 5, 6, 7 (junction_count.py);
  - the junctions with three front-pairs of lowest order counted exactly, M_5 of Lemma mg:lem:jsums5 (N_3, N_2, N'_2 of Lemma mg:lem:threelow,
    three_front_count.py), with the corner factors of M_4 and R as maxima (the second referee's improvement).
It prints delta'/beta^4 against the third-order term beta^3/432 = beta^4/(432 beta) of SU(3), and the three largest parts of delta'/beta^4.
Each beta takes a few minutes (the tree constants on the grid of (bbar, abar))."""
import sys
import third_order_effective as T
from mpmath import mpf
T.THETA7 = True
T.NOVR = {4: {6: 132, 7: 108, 8: 2034, 9: 4344}, 3: {6: 12, 7: 12, 8: 70, 9: 160, 10: 634}}   # k = 9 (d=4), 10 (d=3): referee-section16/ref_fronts.py
T.ZSMALL = {4: {5: 26, 6: 302, 7: 288}, 3: {5: 0, 6: 0, 7: 0}}
T.LOW = {4: (876, 320, 312), 3: (48, 24, 24)}
T.RMAX = True      # the corner factors as maxima (Lemma mg:lem:jsums5)
DEFAULT = {4: ('1e-5', '3e-6', '2.4e-6', '1e-6'), 3: ('1e-4', '3e-5', '2.6e-5', '1e-5')}
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

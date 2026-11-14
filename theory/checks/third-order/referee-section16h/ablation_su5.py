# Referee (Section 16h): ablations with the authors' evaluator (practical_beta.deficit, imported read-only, files unchanged):
#   (1) the value of the (9,2) count (x0, x1, x1e9, x1e20) in d=4, SU(5) and SU(3);
#   (2) SU(5) with u_2 <= e^beta - 1 - beta (Lemma mg:lem:capelement(iv)) instead of Lemma mg:lem:u2sun, d=3 and d=4;
#   (3) the parts of delta' at beta_0(d) for SU(5) (all of them, not only the three largest).
# Usage: python3 ablation_su5.py     (about 3 minutes)
import os, sys, json
here = os.path.dirname(os.path.abspath(__file__)); TO = os.path.join(here, '..')
sys.path.insert(0, TO); sys.path.insert(0, os.path.join(TO, '..', 'upper-bound'))
from mpmath import mp, mpf, log
import practical_beta as PB
import third_order_effective as T
from constants_upper import beta0
from sun_group import sun
C3 = {10: 10, 11: 11, 12: 24, 13: 0, 14: 70}; C4 = {10: 60, 11: 66, 12: 320, 13: 0, 14: 1316, 15: 2520, 16: 12552}

def run(d, path, group, beta, mutate=None, crude_u2=False):
    PB.GROUP = group
    PB.load_terms(d, path, C3 if d == 3 else C4, 15 if d == 3 else 17)
    if mutate: T.EXACT[d] = mutate(list(T.EXACT[d]))
    if crude_u2:
        old = PB.group_data; PB.group_data = lambda: sun(group, sharp_u2=False)
    try:
        D, best, g = PB.deficit(d, beta)
    finally:
        if crude_u2: PB.group_data = old
    parts = {k: float(v/beta**4) for k, v in best['parts'].items()}
    return float(D/beta**4), float(log(best['Ubar']/g['Ups'])/beta**4), parts, best['bars']

if __name__ == "__main__":
    B3, B4 = beta0(3), beta0(4)
    p92 = os.path.join(TO, 'orders', 'd4_s10x_92.json'); p3 = os.path.join(TO, 'orders', 'd3_s9.json')
    for group in (5, 3):
        for fac in (0, 1, 10**9, 10**20):
            mut = lambda ts, f=fac: [((c*f if (n, E) == (9, 2.0) else c), n, E, k, m) for (c, n, E, k, m) in ts]
            D, lu, parts, bars = run(4, p92, group, B4, mut)
            print(f"SU({group}), d=4, S^x_10 u {{(9,2)}} with the (9,2) count x {fac:g}: D/beta^4 = {D:.6g}  (majorant part {parts['multi']:.4g}; bars {bars})")
            sys.stdout.flush()
    for d, path, B in ((3, p3, B3), (4, p92, B4)):
        D, lu, parts, bars = run(d, path, 5, B)
        Dc, luc, partsc, _ = run(d, path, 5, B, crude_u2=True)
        print(f"SU(5), d={d}: with Lemma 16.11: D/beta^4 = {D:.6g}, log(Ubar/Ups)/beta^4 = {lu:.3g};  with u_2 <= e^b-1-b: D/beta^4 = {Dc:.6g},"
              f" log(Ubar/Ups)/beta^4 = {luc:.4g};  target beta^3/6000 = {float(1/(6000*B)):.4g} beta^4")
        print(f"   parts of delta'/beta^4 (Lemma 16.11): " + ', '.join(f'{k} {v:.4g}' for k, v in sorted(parts.items(), key=lambda t: -t[1])))
        sys.stdout.flush()

# Section 16 of Part IV (Remark mg:rem:nextorder (vii)): the SU(3) third-order term on the whole range in d=3.
# Usage: python3 whole_range_d3.py [JOBS]          (about half an hour with 3 jobs)
"""The deficit D of the lower bound (Remark mg:rem:practical (iii)) for SU(3) in d=3 with
  - the set S_9'' of Remark mg:rem:nextorder (vii): orders/d3_s9.json, made by orders/assemble_d3_s9.py from S(8) and the classes of order 9;
  - theta with the closed polymers of 10 to 14 plaquettes counted exactly (10, 11, 24, 0, 70 through a plaquette of Z^3) and crudely from 15 on;
  - Z_(n,4) = 0 for n >= 3 and the exact exits of practical_beta.load_terms;
at a grid of beta from 1e-12 up to beta_0(3).

Interpolation. Fix the parameters found at beta_{i+1} (varpi, bbar, abar, F, F_6 and the decoration constants). Every sum in the hypotheses
of Theorem mg:thm:lower5 and in delta' then has terms that are products of factors increasing with beta, one of them a, b or b_2 (of order beta^2),
so it drops at least by the factor beta/beta_{i+1} (a/beta, b/beta and b_2/beta^2 increase with beta). The right sides change only through Ubar delta'
and Ubar(1 + delta') - a, with Ubar = 1 + beta/6 + O(beta^2), which vary much less. So the parameters and the same delta' remain admissible at every beta <= beta_{i+1}. As 8(K+1)theta and log(Ubar/Ups) also
increase with beta (checked below on a fine grid), D(beta) <= D(beta_{i+1}) there. Hence D < beta^3/432 on [beta_i, beta_{i+1}] once
D(beta_{i+1}) < beta_i^3/432, which is checked for each pair of neighbours.
Printed per point, in units of beta^4: D and its parts 8(K+1)theta, log(Ubar/Ups), log(1+delta'), and the majorant (the part of delta' from the corners)."""
import os, sys
from multiprocessing import Pool
from mpmath import mpf, log, expm1, e
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import practical_beta as PB
import third_order_effective as T
from third_order_effective import dressed

d = 3
CLOSED = {10: 10, 11: 11, 12: 24, 13: 0, 14: 70}
PB.load_terms(d, os.path.join(here, 'orders', 'd3_s9.json'), CLOSED, 15)
sys.path.insert(0, os.path.join(here, '..', 'upper-bound'))
from constants_upper import beta0
B0 = beta0(d)
GRID = [mpf(s) for s in ('1e-12', '1e-11', '1e-10', '1e-9', '1e-8', '1e-7', '1e-6', '3e-6', '1e-5', '2e-5', '5e-5', '9e-5', '1.2e-4', '1.6e-4',
                         '2e-4', '2.3e-4', '2.6e-4', '2.9e-4', '3.2e-4', '3.5e-4', '3.8e-4')] + [B0]

def one(beta):
    D, best, g = PB.deficit(d, beta)
    th = best['loss'] - log(best['Ubar']) - log(1 + best['delta'])          # 8(K+1)theta
    parts = dict(theta=th, ups=log(best['Ubar']/g['Ups']), delta=log(1 + best['delta']), majorant=best['parts']['multi'])
    return beta, D, parts, best['bars']

def theta_parts(beta):
    """8(K+1)theta and log(Ubar/Ups), the two parts of D that do not depend on the parameters"""
    K = 4*(2*d - 3); eps = expm1(2*beta); x = e**2*K*eps
    g = dict(T.GROUPS)['SU(3)'](beta); Ub, s6 = dressed(beta, g)
    th = 2*(d - 2)*e**6*s6 + x**T.CLOSED_FROM/(1 - x) + sum(c*(e*eps)**n for n, c in T.CLOSED_EXACT[d].items())
    return 8*(K + 1)*th, log(Ub/g['Ups'])

if __name__ == "__main__":
    jobs = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    fine = [mpf(10)**(-12 + k*mpf(8.6)/400) for k in range(401)] + [B0]
    fine = sorted(b for b in fine if b <= B0)
    tp = [theta_parts(b) for b in fine]
    mono_theta = all(p[0] <= q[0] for p, q in zip(tp, tp[1:])); mono_ups = all(p[1] <= q[1] for p, q in zip(tp, tp[1:]))
    print(f"8(K+1)theta and log(Ubar/Ups) increasing on {len(fine)} points of [1e-12, beta_0(3)]: {mono_theta}, {mono_ups}"); sys.stdout.flush()
    with Pool(jobs) as pool:
        res = sorted(pool.map(one, sorted(GRID, reverse=True), chunksize=1))
    ok = mono_theta and mono_ups; worst = 0
    for beta, D, parts, bars in res:
        pp = ',  '.join(f'{k} {float(v/beta**4):.3g}' for k, v in parts.items())
        print(f"beta = {float(beta):.6g}: D/beta^4 = {float(D/beta**4):.4g} ({pp}),  beta^3/432 = {float(1/(432*beta)):.4g} beta^4,  (bbar, abar) = {bars}")
    for (lo, Dlo, _, _), (hi, Dhi, _, _) in zip(res, res[1:]):
        r = Dhi/(lo**3/432); worst = max(worst, r); ok = ok and r < 1
        print(f"  [{float(lo):.3g}, {float(hi):.6g}]: D(hi)/(lo^3/432) = {float(r):.3g}")
    print(f"beta_0(3) = {float(B0):.10g}; {len(res)} points; D(beta_(i+1)) < beta_i^3/432 for all neighbours: {ok} (largest ratio {float(worst):.3f})")
    print("D/beta^4 increasing on the grid:", all(p[1]/p[0]**4 <= q[1]/q[0]**4 for p, q in zip(res, res[1:])))

# Section 16 of Part IV (Remark mg:rem:practical (iii)): the third-order lower bound at practical beta.  Usage: python3 practical_beta.py [--group N] [d [beta ...]]
"""The smallest admissible delta' of Theorem mg:thm:lower5 (by characters, varpi = b^2 or b^(5/2)) for SU(3), with
  - theta_7 in place of theta_5 (Lemma mg:lem:theta7);
  - the exact counts of Remark mg:rem:practical (ii): N^v_k (fronts_count.py, d=3: k <= 9, d=4: k <= 8; the referee's ref_fronts.py: N^v_10 = 634
    in d=3, N^v_9 = 4344 in d=4; fronts_decision.py beyond), Z_{n,4} for n = 5, 6, 7 (junction_count.py); with --terms in d=4 (Remark
    mg:rem:nextorder (ix)) also Z_{8,4}, Z_{9,4}, N^v_11 and the exits with three plaquettes (Z4_D4, NV_D4_EXTRA, ZEXITN_D4 below);
  - the junctions with three front-pairs of lowest order counted exactly, M_5 of Lemma mg:lem:jsums5 (N_3, N_2, N'_2 of Lemma mg:lem:threelow,
    three_front_count.py), with the corner factors of M_4 and R as maxima (the second referee's improvement).
It prints the whole deficit D of the lower bound (8(K+1)theta + log(Ubar/Ups) + log(1+delta')) over beta^4 against the third-order term
beta^3/432 = beta^4/(432 beta) of SU(3), and the three largest parts of delta'/beta^4.
Each beta takes a few minutes (the tree constants on the grid of (bbar, abar))."""
import sys
import third_order_effective as T
from mpmath import mpf, log, log1p, exp
T.THETA7 = True
T.NOVR = {4: {6: 132, 7: 108, 8: 2034, 9: 4344, 10: 43650}, 3: {6: 12, 7: 12, 8: 70, 9: 160, 10: 634, 11: 1728, 12: 6369}}   # fronts_decision.py
# The next order (Lemma mg:lem:corners with S_6): counts of the pairs (J,p') by exhaustive search (next_order_count.py), as (count, n, E, |p'|, m):
# weight count b^n a^E e^{m c_J} (f(|p'|) a^{4-|p'|} if |p'| > 4).  NEXT[d] = (N_4, N_{3,6}, N_{3,7}, N_{2,8}, N_{2,444}, N_{2,9}).
NEXT = {4: (42, 680, 344, 992, 4, 256), 3: (0, 24, 24, 48, 0, 16)}
# exact Z_(2,k), k >= 5 (zexit_count.py), in place of the crude counts of two-plaquette exits
T.ZEXIT = {4: {6: 20, 7: 12, 8: 182, 9: 296, 10: 606, 11: 456, 12: 458}, 3: {6: 4, 7: 4, 8: 22, 9: 36, 10: 42, 11: 24, 12: 18}}
def order6(d):
    N3, N2, N2p = T.LOW[d]; N4, N36, N37, N28, N2444, N29 = NEXT[d]
    h = mpf(1)/2
    T.EXACT = {d: [(2*N3, 3, 2, 4, 1), (N2, 2, 3, 4, 1), (N2p, 2, 3 + h, 4, 1),
                   (2*N4, 4, 2, 4, 1), (N36, 3, 3, 4, 1), (N37, 3, 3 + h, 4, 1),
                   (N28, 2, 4, 4, 1), (3*N2444, 2, 4, 4, 2), (N2, 2, 4, 6, 1), (N29, 2, 4 + h, 4, 1)]}
    T.CORNERS = {d: [(5, 2), (4, 3), (3, 4), (2, 5)]}
T.ZSMALL = {4: {5: 26, 6: 302, 7: 288}, 3: {5: 0, 6: 0, 7: 0}}
T.LOW = {4: (876, 320, 312), 3: (48, 24, 24)}
T.RMAX = True      # the corner factors as maxima (Lemma mg:lem:jsums5)
DEFAULT = {4: ('1e-5', '3e-6', '2.6e-6', '1e-6'), 3: ('1e-4', '3e-5', '2.3e-5', '1e-5')}
# The order after in d=3 (Remark mg:rem:nextorder (v)): S_7' = S_6 u {(5,2),(4,3),(4,7/2),(3,4),(3,9/2),(2,5),(2,11/2)}, counts by next_order_count.py;
# theta with the closed sets of 10, 11, 12 plaquettes counted exactly (10, 11, 24 through a plaquette of Z^3: closed_sets_ext.py) and x^13/(1-x) beyond;
# Z_(n,4) = 0 for 3 <= n <= 13 in d=3 (junction_count.py; the planar argument of the remark).
def order7_d3():
    d = 3; h = mpf(1)/2
    order6(d)
    N3, N2, N2p = T.LOW[d]; N4, N36, N37, N28, N2444, N29 = NEXT[d]
    T.EXACT[d] = list(T.EXACT[d]) + [(0, 5, 2, 4, 1), (276, 4, 3, 4, 1), (276, 4, 3 + h, 4, 1), (252, 3, 4, 4, 1), (0, 3, 4, 4, 2), (N36, 3, 4, 6, 1),
                                     (424, 3, 4 + h, 4, 1), (100, 2, 5, 4, 1), (0, 2, 5, 4, 2), (0, 2, 5, 6, 1), (N2p, 2, 5, 7, 1),
                                     (212, 2, 5 + h, 4, 1), (0, 2, 5 + h, 4, 2), (0, 2, 5 + h, 6, 1)]
    T.CORNERS[d] = [(6, 2), (5, 3), (4, 4), (3, 5), (2, 6)]
    T.CLOSED_FROM = 13; T.CLOSED_EXACT = {3: {10: 10, 11: 11, 12: 24}}
    T.ZSMALL = {4: {5: 26, 6: 302, 7: 288}, 3: {n: 0 for n in range(5, 14)}}
# Any order from order_counts.py: --terms FILE.json (terms and corners of S(m)), with --closed n:c,... and --closed-from N for theta (closed
# polymers of n plaquettes counted exactly, crude from N on); in d=3, Z_(n,4) = 0 for all n >= 3 (Remark mg:rem:nextorder (iv)).
Z4_D4 = {5: 26, 6: 302, 7: 288, 8: 1988, 9: 406, 10: 3193, 11: 5911, 12: 26956}   # junction_count.py 3 9 (Z_8: 5 min, Z_9: 25 min); junctions.c 3 12 4 (13 s)
ZEXITN_D4 = {3: {6: 288, 7: 264, 8: 862, 9: 760, 10: 3490, 11: 6596, 12: 14838, 13: 18444, 14: 23672, 15: 15296, 16: 10868},   # zexit_count_n.py 3 3 (90 s)
             4: {6: 152, 7: 40, 8: 3460, 9: 6164, 10: 19664, 11: 28668, 12: 84115, 13: 159020, 14: 356520, 15: 552112, 16: 828023, 17: 862000,
                 18: 843880, 19: 479680, 20: 265315}}                # zexits.c 3 4 (3 min; both layers; as the sixth referee)
NV_D4_EXTRA = {11: 147888, 12: 1121083, 13: 4773186}   # N^v_k of Z^3: fronts_decision.py 3 12 (97 min), fronts_decision.c 3 13 (8 min)
IX_SETTING = False   # --ix: the d=4 setting of Remark mg:rem:nextorder (ix) (N^v_k exact to 11 and crude beyond, Z_(n,4) exact to 9), to reproduce its numbers
def load_terms(d, path, closed=None, closed_from=None):
    import json
    J = json.load(open(path))
    T.EXACT = {d: [tuple(t) for t in J['terms']]}
    T.CORNERS = {d: [tuple(c) for c in J['corners']]}
    if closed: T.CLOSED_EXACT = {d: closed}
    if closed_from: T.CLOSED_FROM = closed_from
    if d == 3: T.ZSMALL = {4: {5: 26, 6: 302, 7: 288}, 3: {n: 0 for n in range(5, 40)}}; T.NV_DECISION = False   # d=3 as refereed (Remark 16.7(vii))
    if d == 4:                     # exact Z_(n,4) in d=4 (junction_count.py, junctions.c), crude beyond; exact N^v_k beyond 10 where counted, and beyond that
        T.ZSMALL = {4: dict(Z4_D4), 3: {5: 0, 6: 0, 7: 0}}     # the decision bound of Lemma mg:lem:frontsdecision in place of (d-1) binom((4d-6)k, k-1)
        T.NOVR = {4: {**T.NOVR[4], **NV_D4_EXTRA}, 3: T.NOVR[3]}
        T.ZEXITN = {4: dict(ZEXITN_D4)}
        T.NV_DECISION = True
        if IX_SETTING:
            T.ZSMALL[4] = {n: c for n, c in Z4_D4.items() if n <= 9}
            T.NOVR = {4: {**{k: v for k, v in T.NOVR[4].items() if k <= 10}, 11: 147888}, 3: T.NOVR[3]}
            T.ZEXITN = {4: {3: dict(ZEXITN_D4[3])}}
            T.NV_DECISION = False
    if d == 3: T.ZEXITN = {3: {3: {6: 20, 7: 20, 8: 46, 9: 52, 10: 150, 11: 260, 12: 382, 13: 360, 14: 288, 15: 144, 16: 76},
                               4: {8: 76, 9: 148, 10: 328, 11: 484, 12: 1143, 13: 2020, 14: 3114, 15: 3672, 16: 3651, 17: 2760, 18: 1780, 19: 780, 20: 315},
                               5: {8: 12, 10: 376, 11: 880, 12: 2271, 13: 4096, 14: 8562, 15: 15532, 16: 25227, 17: 33440, 18: 38064, 19: 35740,
                                   20: 28695, 21: 18600, 22: 10212, 23: 3972, 24: 1296}}}     # zexit_count_n.py
    if d == 3 and CERTIFIED:       # the counts with n = 4, 5 and k >= 17 (one program) -> crude 2 binom(4P_v + T_v n, n) X(n,k) (Theorem mg:thm:thirdcap)
        from junction_constants import counts
        from mpmath import binomial
        Pv, KS, Tv, Nc = counts(3, 'infinite')
        for k in T.NOVR[3]: Nc[k] = mpf(T.NOVR[3][k])
        for n in (4, 5):
            for k in range(17, 4*n + 5):
                X = 0 if k == 5 else min(binomial(4*n + 4, k), (4*n + 4)*Nc[k])
                T.ZEXITN[3][n][k] = max(T.ZEXITN[3][n].get(k, 0), 2*binomial(4*Pv + Tv*n, n)*X)
ORDER = 5
# the group: SU(3) (default; second_order_constants.su3, as refereed) or SU(N) with --group N (sun_group.sun: exact Haar moments, and u_2 by sun_group.u2_bound
# in Lemma mg:lem:dressed); TERM[N] is the size of the third-order term of the gap, |beta^3 coefficient of -log Upsilon| (Corollary mg:cor:third)
GROUP = 3
TERM = {3: mpf(1)/432, 5: mpf(1)/6000}
BOUNDS = False      # --bounds: the group data replaced by proven bounds (rigorous_bounds.bound_group; Theorem mg:thm:thirdcap), D by its bound D_bound
CERTIFIED = False   # --certified: only counts made by two independent programs; in d=3 the exits with 4, 5 plaquettes and k >= 17 edges (one program) become crude
def group_data():
    if BOUNDS:
        from rigorous_bounds import bound_group
        return bound_group(GROUP, 0, sharp_u2=True)
    if GROUP == 3: return dict(T.GROUPS)['SU(3)']
    from sun_group import sun
    return sun(GROUP, sharp_u2=True)
def deficit(d, beta):
    """the whole deficit D of the lower bound at beta for the group (the smaller delta' of varpi = b^2, b^(5/2)), the evaluation, and the group data"""
    g = group_data()(beta); best = None
    for vp in (2, mpf(5)/2):
        try: lo = T.lower3e(beta, d, g, varpi_pow=vp, chars=True)
        except AssertionError: continue
        if best is None or lo['delta'] < best['delta']: best = lo
    if BOUNDS: return best['D_bound'], best, g
    return best['loss'] - log(g['Ups']), best, g
def run(d, betas):
    """The criterion is the whole deficit of the lower bound, D = 8(K+1)theta + log(Ubar/Ups) + log(1+delta'), against beta^3/432: then
    m - 4log(1/u) + log(Ups) >= -D > -beta^3/432."""
    if ORDER == 6: order6(d)
    if ORDER == 7: order7_d3()
    for be in betas:
        beta = mpf(be); D, best, g = deficit(d, beta)
        target = TERM[GROUP]*beta**3
        lth = best['loss'] - log(best['Ubar']) - log1p(best['delta'])
        if BOUNDS: g = dict(g, Ups=best['Ubar']/(1 + exp(2*beta)*g['Delta']/best['Ubar']))   # for the printout only: log(Ubar/Ups) -> the bound e^{2beta} Delta
        parts = sorted(((float(v/beta**4), k) for k, v in best['parts'].items()), reverse=True)[:3]
        print(f"d={d} beta={be}: D/beta^4 = {float(D/beta**4):.4g} = theta-loss {float(lth/beta**4):.3g} + log(Ubar/Ups) {float(log(best['Ubar']/g['Ups'])/beta**4):.2g}"
              f" + log(1+delta') {float(log1p(best['delta'])/beta**4):.4g};  beta^3/{int(round(float(1/TERM[GROUP])))} = {float(target/beta**4):.4g} beta^4 ({'resolved' if D < target else 'not resolved'});"
              f"  largest parts of delta': " + ', '.join(f'{k} {a:.3g}' for a, k in parts))
        sys.stdout.flush()
DEFAULT6 = {4: ('1.3e-5', '1.2e-5', '1e-5', '3e-6'), 3: ('1e-4', '9.5e-5', '9e-5', '3e-5')}
if __name__ == "__main__":
    # default: Remark mg:rem:practical(iii) (Lemma mg:lem:jsums5, theta_7);  --next: the next order S_6 of Lemma mg:lem:corners with theta_10;
    # --next-theta7: S_6 with theta_7
    args = sys.argv[1:]
    if '--group' in args: i = args.index('--group'); GROUP = int(args[i + 1]); del args[i:i + 2]
    if '--bounds' in args: BOUNDS = True; args.remove('--bounds')
    if '--certified' in args: CERTIFIED = True; args.remove('--certified')
    if '--ix' in args: IX_SETTING = True; args.remove('--ix')
    if '--next' in args: ORDER = 6; T.CLOSED_FROM = 10; args.remove('--next')
    elif '--next-theta7' in args: ORDER = 6; args.remove('--next-theta7')
    elif '--next2' in args: ORDER = 7; args.remove('--next2')     # d=3 only
    elif '--terms' in args:
        i = args.index('--terms'); path = args[i + 1]; del args[i:i + 2]
        closed = None; cf = None
        if '--closed' in args:
            i = args.index('--closed'); closed = {int(x.split(':')[0]): int(x.split(':')[1]) for x in args[i + 1].split(',')}; del args[i:i + 2]
        if '--closed-from' in args:
            i = args.index('--closed-from'); cf = int(args[i + 1]); del args[i:i + 2]
        ORDER = 0; load_terms(int(args[0]), path, closed, cf)
    else: ORDER = 5
    ds = [int(args[0])] if args else [4, 3]
    if ORDER == 7: ds = [3]
    for d in ds: run(d, args[1:] if len(args) > 1 else {5: DEFAULT, 6: DEFAULT6, 7: {3: ('2e-4', '1.8e-4', '1e-4')}, 0: {3: ('3e-4',), 4: ('1e-5',)}}[ORDER][d])

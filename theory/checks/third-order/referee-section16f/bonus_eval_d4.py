# Referee (Section 16, sixth pass), for the next step only (not a check of the paper): D/beta^4 in d=4 with the revised set (orders/d4_s8b.json) and, in
# addition, the referee's N^v_12 = 1121083, N^v_13 = 4773186 (fronts_enum.c) and Z_(10,4) = 3193, Z_(11,4) = 5911, Z_(12,4) = 26956 (junc_partition.c),
# and then also N^v_14 = 31228116, N^v_15 = 149251700 (fronts_enum.c 3 15).
# Usage: python3 bonus_eval_d4.py beta [beta ...]
import os, sys, copy
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
from mpmath import mpf
import practical_beta as PB
import third_order_effective as T
NOVR0 = copy.deepcopy(T.NOVR)
for be in sys.argv[1:]:
    beta = mpf(be)
    for name, extra in (("as in the paper", 0), ("+ N^v_12, N^v_13, Z_(10..12,4)", 1), ("+ also N^v_14, N^v_15", 2)):
        T.NOVR = copy.deepcopy(NOVR0); T.CLOSED_EXACT = None
        PB.load_terms(4, os.path.join(here, '..', 'orders', 'd4_s8b.json'), {10: 60, 11: 66}, 12)
        if extra:
            T.NOVR[4].update({12: 1121083, 13: 4773186})
            T.ZSMALL[4].update({10: 3193, 11: 5911, 12: 26956})
        if extra == 2:
            T.NOVR[4].update({14: 31228116, 15: 149251700})
        D, best, g = PB.deficit(4, beta)
        p = {k: float(v/beta**4) for k, v in best['parts'].items()}
        print(f"beta = {be}: {name:32s} D/beta^4 = {float(D/beta**4):.3f} against {float(1/(432*beta)):.2f};  multi {p['multi']:.3g}, one_plaq {p['one_plaq']:.3g}, "
              f"exits {p['exits']:.3g}, four {p['four']:.3g}"); sys.stdout.flush()

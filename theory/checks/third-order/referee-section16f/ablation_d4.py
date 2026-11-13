# Referee (Section 16, sixth pass): which ingredients of Remark mg:rem:nextorder(ix) carry the d=4 threshold, with the authors' evaluator.
# Usage: python3 ablation_d4.py [beta] [--terms FILE]        (default 4.3e-5 and orders/d4_s8.json; each variant about 20 s)
"""D/beta^4 for SU(3), d=4, at beta (default 4.3e-5), with the setting of practical_beta.py --terms orders/d4_s8.json --closed 10:60,11:66 --closed-from 12,
and with one ingredient removed (replaced by the crude bound used before) or inflated (to test robustness against undercounts).  The criterion is
D < beta^3/432."""
import os, sys, copy
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
from mpmath import mpf, log
import practical_beta as PB
import third_order_effective as T

TERMS = os.path.join(here, '..', 'orders', 'd4_s8.json'); M7 = os.path.join(here, '..', 'orders', 'd4_m7.json')
if '--terms' in sys.argv:
    i = sys.argv.index('--terms'); TERMS = sys.argv[i + 1]; del sys.argv[i:i + 2]
beta = mpf(sys.argv[1]) if len(sys.argv) > 1 else mpf('4.3e-5')
NOVR0 = copy.deepcopy(T.NOVR)

def setup(terms=TERMS, closed={10: 60, 11: 66}, cf=12):
    T.NOVR = copy.deepcopy(NOVR0); T.CLOSED_EXACT = None; T.CLOSED_FROM = 8
    PB.load_terms(4, terms, closed, cf)

def run(name, tweak=None, **kw):
    setup(**kw)
    if tweak: tweak()
    D, best, g = PB.deficit(4, beta)
    parts = {k: float(v/beta**4) for k, v in best['parts'].items()}
    lth = float((best['loss'] - log(best['Ubar']) - log(1 + best['delta']))/beta**4)
    tgt = float(1/(432*beta))
    print(f"{name:58s} D/beta^4 = {float(D/beta**4):9.3f}  (target {tgt:.2f}: {'resolved' if D/beta**4 < tgt else 'NOT resolved'});"
          f"  theta-loss {lth:.3g}, multi {parts['multi']:.3g}, one_plaq {parts['one_plaq']:.3g}, exits {parts['exits']:.3g}, four {parts['four']:.3g}, two {parts['two']:.2g}")
    sys.stdout.flush()

def drop_z(*ns):
    def f():
        for n in ns: T.ZSMALL[4].pop(n, None)
    return f
def no_nv11():
    T.NOVR[4].pop(11, None)
def no_zexit3():
    T.ZEXITN = None
def scale_exact(fac):
    def f():
        T.EXACT[4] = [(c*fac, n, E, k, m) for (c, n, E, k, m) in T.EXACT[4]]
    return f
def inflate_counts(fac):
    def f():
        for n in (8, 9): T.ZSMALL[4][n] *= fac
        T.NOVR[4][11] *= fac
        T.ZEXITN = {4: {3: {k: c*fac for k, c in T.ZEXITN[4][3].items()}}}
    return f

print(f"SU(3), d=4, beta = {float(beta):g}, terms {os.path.basename(TERMS)}")
run("baseline (as in the paper)")
run("theta crude from 10 plaquettes (no exact 10, 11)", closed=None, cf=10)
run("theta: 10 exact, crude from 11", closed={10: 60}, cf=11)
run("theta: 10, 11 counted 1000x, crude from 12", closed={10: 60000, 11: 66000}, cf=12)
run("Z_(9,4) crude", tweak=drop_z(9))
run("Z_(8,4), Z_(9,4) crude", tweak=drop_z(8, 9))
run("N^v_11 crude", tweak=no_nv11)
run("Z_(3,k) crude", tweak=no_zexit3)
run("S(7) instead (corners (2,6),(3,5),(4,4),(5,3),(6,2))", terms=M7)
run("all exact pair counts x 10", tweak=scale_exact(10))
run("all exact pair counts x 100", tweak=scale_exact(100))
run("all exact pair counts x 1e6", tweak=scale_exact(10**6))
run("Z_8, Z_9, N^v_11, Z_(3,k) all x 1000", tweak=inflate_counts(1000))

# Referee (Section 16, seventh pass): which ingredients of Remark mg:rem:nextorder(x) carry the bound at beta_0(4), with the authors' evaluator.
# Usage: python3 ablation_x.py [beta] [variant ...]       (default beta_0(4) = 1.53766448e-4 and all variants; about 20 s per variant)
"""D/beta^4 for SU(3), d = 4, in the setting of practical_beta.py --terms orders/d4_s10x.json --closed 10:60,...,16:12552 --closed-from 17, with one
ingredient removed (replaced by what was used before), inflated (robustness against undercounts) or changed.  The terms are rebuilt from the class
caches by the referee's own assembler (assemble_check_x.py), so that single classes or groups of classes can be scaled or left out.  The criterion
at beta_0(4) is D/beta^4 < 1/(432 beta) = 15.05.  Also printed: the chosen (bbar, abar), and which corner sets the maximum in M_S."""
import os, sys, json, copy
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..')); sys.path.insert(0, here)
from fractions import Fraction as Fr
from mpmath import mpf, log
import practical_beta as PB
import third_order_effective as T
import assemble_check_x as A

SCR = os.environ.get('REFSCR', '/tmp')
beta = mpf(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1][0].isdigit() else mpf('1.53766448e-4')
only = [a for a in sys.argv[1:] if not a[0].isdigit()]
CLOSED = {10: 60, 11: 66, 12: 320, 13: 0, 14: 1316, 15: 2520, 16: 12552}
NOVR0 = copy.deepcopy(T.NOVR); BEST, _ = A.load_caches(); DB0 = T.decision_bound

def write_terms(inS=A.in_S, scale=None, name='terms'):
    need = A.needed_classes(inS)
    T_ = A.terms_of(BEST, inS, need, scale)
    C = [(n, float(E)) for n, E in A.corners(inS)]
    path = os.path.join(SCR, f'abl_{name}.json'); json.dump({'terms': T_, 'corners': C}, open(path, 'w'))
    return path

BASE = write_terms(name='base')

def run(name, tweak=None, terms=BASE, closed=CLOSED, cf=17):
    if only and not any(o in name for o in only): return
    T.NOVR = copy.deepcopy(NOVR0); T.CLOSED_EXACT = None; T.CLOSED_FROM = 8; T._DEC.clear(); T.decision_bound = DB0
    PB.load_terms(4, terms, dict(closed) if closed else None, cf)
    if tweak: tweak()
    D, best, g = PB.deficit(4, beta)
    parts = {k: float(v/beta**4) for k, v in best['parts'].items()}
    lth = float((best['loss'] - log(best['Ubar']) - log(1 + best['delta']))/beta**4)
    tgt = float(1/(432*beta))
    bars = best['bars']; corner = ''
    if len(bars) == 3 and bars[2] == 'corners':
        b, a = best['b'], best['a']
        fac = sorted(((float((b/bars[0])**n*(a/bars[1])**E), (n, E)) for n, E in T.CORNERS[4]), reverse=True)
        corner = f"  max corner {fac[0][1]} ({fac[0][0]:.2g}; next {fac[1][1]} {fac[1][0]:.2g})"
    bstr = '(' + ', '.join(f'2^{float(log(x, 2)):.0f}' if not isinstance(x, str) else x for x in bars) + ')'
    print(f"{name:62s} D/b^4 = {float(D/beta**4):9.3f} ({'ok' if D/beta**4 < tgt else 'FAILS'} vs {tgt:.2f}); multi {parts['multi']:.3g}, exits {parts['exits']:.3g},"
          f" one_plaq {parts['one_plaq']:.3g}, four {parts['four']:.3g}, two {parts['two']:.2g}, theta {lth:.2g}; bars {bstr}{corner}")
    sys.stdout.flush()

def drop_novr(*ks):
    def f():
        for k in ks: T.NOVR[4].pop(k, None)
    return f
def no_decision(): T.NV_DECISION = False
def no_decision_no_novr():
    T.NV_DECISION = False
    for k in (12, 13): T.NOVR[4].pop(k, None)
def fixed_s(s=mpf('0.194')):
    def f():
        def db(D, k, grid=400):
            y = max((1 + s)**(2*D - 1) - 1, (1 + s)**(D - 1))
            return ((1 + s)**(2*D) - 1 - 2*D*s)*(y/s)**k
        T.decision_bound = db; T._DEC.clear()
    return f
def drop_z(*ns):
    def f():
        for n in ns: T.ZSMALL[4].pop(n, None)
    return f
def scale_z(fac, ns=(10, 11, 12)):
    def f():
        for n in ns: T.ZSMALL[4][n] *= fac
    return f
def drop_zexit(n):
    def f(): T.ZEXITN[4].pop(n, None)
    return f
def scale_zexit(n, fac):
    def f(): T.ZEXITN[4][n] = {k: c*fac for k, c in T.ZEXITN[4][n].items()}
    return f
def scale_novr(fac, ks=(12, 13)):
    def f():
        for k in ks: T.NOVR[4][k] *= fac
    return f
def scale_exact(fac):
    def f(): T.EXACT[4] = [(c*fac, n, E, k, m) for (c, n, E, k, m) in T.EXACT[4]]
    return f

def S_layers(layers, extra=((3, Fr(6)),), exclude=((5, Fr(11, 2)),)):
    return lambda n, E: A.in_S(n, E, layers=layers, extra=extra, exclude=exclude)
def order_of(K, n):   # the smallest order m with a pair of class K at n: n + E_min <= m + 1/2
    import math
    return math.ceil(n + min(A.E_of(K, k) for k in set(K)) - Fr(1, 2))

if __name__ == "__main__":
    print(f"SU(3), d=4, beta = {float(beta):.9g}; target 1/(432 beta) = {float(1/(432*beta)):.4f}")
    run("baseline (as in the paper)")
    run("N^v_k: decision bound off (binomial beyond 13)", no_decision)
    run("N^v_12, N^v_13 exact dropped (decision bound from 12)", drop_novr(12, 13))
    run("N^v_k: exact to 11, binomial beyond (setting of (ix))", no_decision_no_novr)
    run("N^v_k: decision bound with s = 0.194 fixed (as quoted)", fixed_s())
    run("N^v_12, N^v_13 x 100", scale_novr(100))
    run("Z_(n,4), n = 10..12, crude", drop_z(10, 11, 12))
    run("Z_(n,4), n = 10..12, x 100", scale_z(100))
    run("Z_(n,4), n = 8..12, x 100", scale_z(100, (8, 9, 10, 11, 12)))
    run("Z_(4,k) crude", drop_zexit(4))
    run("Z_(4,k) x 10", scale_zexit(4, 10))
    run("Z_(4,k) x 100", scale_zexit(4, 100))
    run("Z_(3,k) x 100", scale_zexit(3, 100))
    run("theta: closed crude from 12 (10, 11 exact)", closed={10: 60, 11: 66}, cf=12)
    run("theta: closed crude from 14 (10..13 exact)", closed={10: 60, 11: 66, 12: 320, 13: 0}, cf=14)
    run("theta: closed crude from 15", closed={10: 60, 11: 66, 12: 320, 13: 0, 14: 1316}, cf=15)
    run("theta: closed crude from 16", closed={10: 60, 11: 66, 12: 320, 13: 0, 14: 1316, 15: 2520}, cf=16)
    run("theta: closed counts 12..16 x 100", closed={10: 60, 11: 66, 12: 32000, 13: 0, 14: 131600, 15: 252000, 16: 1255200})
    run("theta: closed counts 10..16 x 10^4", closed={n: c*10**4 for n, c in CLOSED.items()})
    run("set: S(8) u L9 (no (3,6), no L10) = orders/d4_s9.json", terms=write_terms(S_layers(((8, 2), (9, 4)), extra=(), exclude=()), name='s9'))
    run("set: S(8) u L9 u (3,6) (no L10)", terms=write_terms(S_layers(((8, 2), (9, 4)), exclude=()), name='s9x'))
    run("set: S^x_10 without (3,6)", terms=write_terms(S_layers(((8, 2), (9, 4), (10, 5)), extra=()), name='nox36'))
    run("set: S(8) u L10 (no L9 beyond, i.e. L9 n>=5 only) u (3,6)", terms=write_terms(S_layers(((8, 2), (9, 5), (10, 5))), name='noL9n4'))
    run("set: S^x_10 without (8,2)", terms=write_terms(lambda n, E: A.in_S(n, E) and (n, E) != (8, 2), name='no82'))
    run("set: S^x_10 without (7,3),(8,2)", terms=write_terms(lambda n, E: A.in_S(n, E) and (n, E) not in ((8, 2), (7, 3)), name='no7382'))
    run("all exact pair counts x 10", scale_exact(10))
    run("all exact pair counts x 10^3", scale_exact(10**3))
    run("all exact pair counts x 10^6", scale_exact(10**6))
    run("order-10 classes (pairs at n+E = 10, 10.5) x 10^3", terms=write_terms(scale=lambda K, n: 1000 if order_of(K, n) >= 10 else 1, name='o10x'))
    run("order-9 classes x 10^3", terms=write_terms(scale=lambda K, n: 1000 if order_of(K, n) == 9 else 1, name='o9x'))
    run("order <= 8 classes x 10^3", terms=write_terms(scale=lambda K, n: 1000 if order_of(K, n) <= 8 else 1, name='o8x'))
    run("order-9 and order-10 classes x 10^9", terms=write_terms(scale=lambda K, n: 10**9 if order_of(K, n) >= 9 else 1, name='o910x'))
    run("set: S(8) u {(7,2), (8,2)} (L9, L10 reduced to the class 4,4 at n = 7, 8)", terms=write_terms(lambda n, E: A.in_S(n, E, layers=((8, 2),), extra=((7, Fr(2)), (8, Fr(2))), exclude=()), name='s8_7282'))
    run("set: S(8) u {(6,3), (7,2), (8,2)}", terms=write_terms(lambda n, E: A.in_S(n, E, layers=((8, 2),), extra=((6, Fr(3)), (7, Fr(2)), (8, Fr(2))), exclude=()), name='s8_637282'))
    def scale_high(fac):
        def f(): T.EXACT[4] = [(c*fac if n + E >= 5.5 else c, n, E, k, m) for (c, n, E, k, m) in T.EXACT[4]]
        return f
    run("exact pairs with n + E >= 11/2 x 10^6 (all but the lowest order, Lemma 16.3)", scale_high(10**6))

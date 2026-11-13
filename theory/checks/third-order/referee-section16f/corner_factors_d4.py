# Referee (Section 16, sixth pass): the corner factors of S''_8 in d=4 at the evaluated point, and what moving one corner would give.
# Usage: python3 corner_factors_d4.py [beta] [--scanA beta1,beta2,...]          (default 4.3e-5)
"""With the authors' evaluator (setting of Remark mg:rem:nextorder(ix)): b, a and the chosen (bbar, abar); the factor (b/bbar)^n (a/abar)^E of each corner of
S''_8, i.e. which corner sets the majorant.  Then two hypothetical sets (not claims of the paper; for the next step only):
  (A) S''_8 plus the pairs at (3,5), (3,11/2), with the referee's counts of the order-8 classes with three plaquettes (junc_partition_d4_order8_n3.txt);
  (B) S''_8 plus (2,6), (2,13/2) with no pairs at all (a lower bound for D after counting them, if their counts were zero)."""
import os, sys, copy, json
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
from mpmath import mpf
from fractions import Fraction as Fr
import practical_beta as PB
import third_order_effective as T
SCAN = []
if '--scanA' in sys.argv:
    i = sys.argv.index('--scanA'); SCAN = sys.argv[i + 1].split(','); del sys.argv[i:i + 2]
beta = mpf(sys.argv[1]) if len(sys.argv) > 1 else mpf('4.3e-5')
TERMS = os.path.join(here, '..', 'orders', 'd4_s8.json')
NOVR0 = copy.deepcopy(T.NOVR)

def corners_of(inS):
    EXP = [(n, Fr(e, 2)) for n in range(2, 14) for e in [4] + list(range(6, 30)) if not (n == 2 and e < 6)]
    out = [(n, E) for n, E in EXP if not inS(n, E)]
    return sorted((n, float(E)) for n, E in out if not any((n2, E2) != (n, E) and n2 <= n and E2 <= E for n2, E2 in out))

def run(name, extra_terms=(), corners=None):
    T.NOVR = copy.deepcopy(NOVR0); T.CLOSED_EXACT = None
    PB.load_terms(4, TERMS, {10: 60, 11: 66}, 12)
    T.EXACT[4] = list(T.EXACT[4]) + list(extra_terms)
    if corners is not None: T.CORNERS[4] = corners
    D, best, g = PB.deficit(4, beta)
    b, a, bars = best['b'], best['a'], best['bars']
    print(f"{name}: D/beta^4 = {float(D/beta**4):.3f} (target {float(1/(432*beta)):.2f}); multi {float(best['parts']['multi']/beta**4):.3g}; b = {float(b):.4g}, a = {float(a):.4g}, "
          f"(bbar, abar) = (2^{int(round(float(__import__('mpmath').log(bars[0], 2))))}, 2^{int(round(float(__import__('mpmath').log(bars[1], 2))))})" if len(bars) == 3 else f"{name}: {bars}")
    if len(bars) == 3:
        bb, ab = bars[0], bars[1]
        fac = sorted(((float((b/bb)**n*(a/ab)**E), (n, E)) for (n, E) in T.CORNERS[4]), reverse=True)
        print("   corner factors (b/bbar)^n (a/abar)^E, largest first: " + ", ".join(f"{c}: {v:.3g}" for v, c in fac))
    sys.stdout.flush()

S8 = lambda n, E: n + E <= Fr(15, 2) or (n + E <= Fr(17, 2) and n >= 4)
if not SCAN: run("S''_8 (as in the paper)")
# (A): the pairs at n = 3, E = 5, 11/2 (classes with further fronts K; p' of size k; e^{(|K|-1) c_J}; multiplicity = number of fronts of size k in K)
h = mpf(1)/2
A_terms = [(45700, 3, 5, 4, 1), (2*3028, 3, 5, 6, 1), (2*9768, 3, 5, 4, 2), (344, 3, 5, 7, 1),
           (51640, 3, 5 + h, 4, 1), (5712, 3, 5 + h, 6, 1), (2*9600, 3, 5 + h, 4, 2)]
SA = lambda n, E: S8(n, E) or (n == 3 and E <= Fr(11, 2))
if SCAN:
    for be in SCAN:
        beta = mpf(be); run(f"(A) at beta = {be}", A_terms, corners_of(SA))
    sys.exit(0)
run("(A) S''_8 + (3,5), (3,11/2) counted", A_terms, corners_of(SA))
SB = lambda n, E: S8(n, E) or (n == 2 and E <= Fr(13, 2))
run("(B) S''_8 + (2,6), (2,13/2) as if empty", (), corners_of(SB))
SAB = lambda n, E: SA(n, E) or SB(n, E)
run("(A)+(B)", A_terms, corners_of(SAB))

# Referee (Section 16, seventh pass): which exponents of S^x_10 are needed at beta_0(4), and how much each class may be undercounted.
# Usage: python3 essential_x.py [beta]          (default beta_0(4); about 20 s per variant)
"""(1) For every exponent (n, E) of S^x_10 with n + E >= 8 (and for the whole layers), D/beta^4 with that exponent removed from S (Lemma mg:lem:corners
allows any finite S; the corners are recomputed), against 1/(432 beta).  (2) For every class K counted at order 9 or 10 (n + E_min >= 9), the factor by
which its counts may be multiplied before D/beta^4 reaches 1/(432 beta) (bisection on a logarithmic scale, up to 10^9): a measure of how wrong a
count could be without changing the conclusion."""
import os, sys, copy
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..')); sys.path.insert(0, here)
from fractions import Fraction as Fr
from mpmath import mpf
import practical_beta as PB
import third_order_effective as T
import assemble_check_x as A
import ablation_x as AB      # reuses its setup (prints the baseline once)

beta = AB.beta; tgt = 1/(432*beta)
def D_of(terms):
    T.NOVR = copy.deepcopy(AB.NOVR0); T.CLOSED_EXACT = None; T.CLOSED_FROM = 8; T._DEC.clear(); T.decision_bound = AB.DB0
    PB.load_terms(4, terms, dict(AB.CLOSED), 17)
    D, best, g = PB.deficit(4, beta)
    return D/beta**4

S_all = [(n, E) for n, E in A.all_exponents() if A.in_S(n, E)]
print(f"S^x_10 has {len(S_all)} exponents; testing those with n + E >= 8 one at a time")
for n, E in sorted(S_all, key=lambda t: (t[0] + t[1], t[0])):
    if n + E < 8: continue
    path = AB.write_terms(lambda m, F, n=n, E=E: A.in_S(m, F) and (m, F) != (n, E), name=f'drop_{n}_{float(E)}')
    v = D_of(path)
    print(f"  without ({n}, {E}): D/beta^4 = {float(v):8.3f}  {'needed' if v >= tgt else 'not needed'}"); sys.stdout.flush()

def order_min(K, n): return n + min(A.E_of(K, k) for k in set(K))
need = A.needed_classes(A.in_S)
for K in sorted(need, key=lambda K: (len(K), K)):
    ns = [n for n in range(2, need[K] + 1) if order_min(K, n) >= 9]
    if not ns: continue
    lo, hi = 0.0, 9.0
    for _ in range(7):
        mid = (lo + hi)/2
        path = AB.write_terms(scale=lambda K2, n2, K=K, ns=ns, f=10**mid: f if (K2 == K and n2 in ns) else 1, name='scl')
        if D_of(path) < tgt: lo = mid
        else: hi = mid
    print(f"  class {list(K)} at n = {ns}: counts may be multiplied by up to ~10^{lo:.1f} before D reaches the term"); sys.stdout.flush()

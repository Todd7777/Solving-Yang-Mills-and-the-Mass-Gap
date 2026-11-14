# Referee (Section 16, seventh pass): the set S^x_10 of Remark mg:rem:nextorder(x), its corners and its terms, rebuilt from the definition of E
# without the authors' assemble_set.py / order_counts.classes.  Usage: python3 assemble_check_x.py [--dump OUT.json]
"""Exponents: n >= 2, E in {2, 3, 7/2, 4, 9/2, ...}, E >= 3 if n = 2 (Lemma mg:lem:corners).  A junction with further front-pairs of sizes K (a multiset,
|K| >= 2, sizes 4 or >= 6) and n plaquettes gives, for each further front-pair p' of size k, a pair (J, p') with E = (k - 4) + (sum K - k)/2.
S^x_10 = S(8) u {n+E <= 19/2, n >= 4} u {(3,6)} u {n+E <= 21/2, n >= 5} minus {(5,11/2)}.
Checks: (1) the corners are the minimal exponents outside S; (2) every class K with a pair in S at n is counted (in one of the caches) to at least n;
(3) the terms (count, n, E, |p'|, |K|-1) rebuilt here equal orders/d4_s10x.json; (4) the caches agree with each other where they overlap."""
import json, os, sys, itertools
from fractions import Fraction as Fr
from collections import Counter, defaultdict

here = os.path.dirname(os.path.abspath(__file__)); orders = os.path.join(here, '..', 'orders')
CACHES = [os.path.join(orders, f) for f in ('d4_m8_classes.json', 'd4_m9_classes.json', 'd4_m10_classes.json')]

def in_E(n, E):
    return n >= 2 and (E == 2 or E >= 3) and (E * 2).denominator == 1 and (n > 2 or E >= 3)

def in_S(n, E, layers=((8, 2), (9, 4), (10, 5)), extra=((3, Fr(6)),), exclude=((5, Fr(11, 2)),)):
    if not in_E(n, E): return False
    if (n, E) in exclude: return False
    if (n, E) in extra: return True
    return any(n + E <= Fr(2 * m + 1, 2) and n >= nmin for m, nmin in layers)

def E_of(K, k):
    return Fr(k - 4) + Fr(sum(K) - k, 2)

def all_exponents(nmax=14, emax=12):
    for n in range(2, nmax + 1):
        for twoE in range(4, 2 * emax + 1):
            E = Fr(twoE, 2)
            if in_E(n, E): yield (n, E)

def corners(inS):
    out = [(n, E) for n, E in all_exponents() if not inS(n, E)]
    return sorted((n, E) for n, E in out if not any((m, F) != (n, E) and m <= n and F <= E for m, F in out))

def needed_classes(inS, maxsize=24, maxr=7, nmax=12):
    """every multiset K (sizes 4, 6..maxsize; 2 <= |K| <= maxr) and the largest n with a pair (n, E(K,k)) in S"""
    sizes = [4] + list(range(6, maxsize + 1))
    need = {}
    for r in range(2, maxr + 1):
        for K in itertools.combinations_with_replacement(sizes, r):
            if min(E_of(K, k) for k in set(K)) > 9: continue
            ns = [n for n in range(1, nmax + 1) for k in set(K) if inS(n, E_of(K, k))]
            if ns: need[K] = max(ns)
    return need

def load_caches(paths=CACHES):
    best = {}; clash = []
    for p in paths:
        for key, cnt in json.load(open(p)).items():
            K = tuple(json.loads(key.split('|')[0])); nm = int(key.split('|')[1]); cnt = {int(n): v for n, v in cnt.items()}
            if K in best:
                nm0, c0 = best[K]
                for n in range(1, min(nm, nm0) + 1):
                    if c0.get(n, 0) != cnt.get(n, 0): clash.append((K, n, c0.get(n, 0), cnt.get(n, 0)))
                if nm > nm0: best[K] = (nm, cnt)
            else: best[K] = (nm, cnt)
    return best, clash

def terms_of(best, inS, need, scale=None):
    """the aggregated terms; scale(K, n) -> factor for the count of class K at n (ablation)"""
    agg = Counter()
    for K in need:
        nm, cnt = best[K]
        for n, N in cnt.items():
            if not N: continue
            f = scale(K, n) if scale else 1
            for k, mult in Counter(K).items():
                E = E_of(K, k)
                if inS(n, E): agg[(n, float(E), k, len(K) - 1)] += mult * N * f
    return [[c, n, E, k, r] for (n, E, k, r), c in sorted(agg.items())]

if __name__ == '__main__':
    C = corners(in_S)
    print('corners of S^x_10:', [(n, str(E)) for n, E in C])
    need = needed_classes(in_S)
    best, clash = load_caches()
    print(f'{len(need)} classes have a pair in S^x_10; cache clashes: {clash or "none"}')
    missing = [(K, n, best.get(K, (-1,))[0]) for K, n in sorted(need.items()) if best.get(K, (-1,))[0] < n]
    print('classes not counted far enough:', missing or 'none')
    for K, n in sorted(need.items(), key=lambda t: (len(t[0]), t[0])):
        pairs = sorted({(m, str(E_of(K, k))) for m in range(2, n + 1) for k in set(K) if in_S(m, E_of(K, k))})
        print(f'  K={list(K)}: needed to n={n}, counted to {best[K][0]}; pairs in S: {pairs}')
    T = terms_of(best, in_S, need)
    ref = json.load(open(os.path.join(orders, 'd4_s10x.json')))
    same_terms = sorted(map(tuple, T)) == sorted(tuple(t) for t in ref['terms'])
    same_corners = [(n, float(E)) for n, E in C] == [tuple(c) for c in ref['corners']]
    print(f'{len(T)} terms; equal to orders/d4_s10x.json: terms {same_terms}, corners {same_corners}')
    if '--dump' in sys.argv:
        json.dump({'terms': T, 'corners': [(n, float(E)) for n, E in C]}, open(sys.argv[sys.argv.index('--dump') + 1], 'w'))

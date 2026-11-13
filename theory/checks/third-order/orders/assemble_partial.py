# Section 16 of Part IV (Lemma mg:lem:corners): a set S(m) plus the pairs of the next layer with at least NMIN plaquettes, and its corners.
# Usage: python3 orders/assemble_partial.py base.json cache.json m NMIN out.json
"""base.json: the terms of S(m) = {n + E <= m + 1/2} (order_counts.py or order_counts_par.py D m).  cache.json: counts of classes of order m+1, keyed
"[k1, k2, ...]|nmax" -> {n: count} as in order_counts_par.py.  The new set is S = S(m) u {(n,E): m + 1/2 < n + E <= m + 3/2, n >= NMIN}.  The script
checks that the classes in the cache are all the classes of order m+1 (order_counts.classes) with such a pair, each counted far enough, adds their
pairs (count, n, E, |p'|, |K| - 1) as order_counts.py does, and computes the corners as the minimal exponents outside S (the exponent set of Lemma
mg:lem:corners: n >= 2, E in {2, 3, 7/2, 4, ...}, E >= 3 if n = 2).  Used for S''_9 in d=3 (as assemble_d3_s9.py) and S''_8 in d=4."""
import json, os, sys
from collections import Counter
from fractions import Fraction
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..'))
import order_counts as OC

def assemble(base, cache, m, nmin):
    lo, hi = Fraction(2*m + 1, 2), Fraction(2*m + 3, 2)
    counts, counted_to = {}, {}
    for key, cnt in cache.items():
        K = tuple(json.loads(key.split('|')[0])); counts[K] = {int(n): v for n, v in cnt.items()}; counted_to[K] = int(key.split('|')[1])
    need = {}
    for K, nmax in OC.classes(m + 1).items():
        ns = [n for n in range(nmin, nmax + 1) for k in set(K) if lo < n + OC.E_of(K, k) <= hi]
        if ns: need[K] = max(ns)
    missing = [K for K in need if K not in counts or counted_to[K] < need[K]]
    assert not missing, missing
    agg = Counter()
    for t in base['terms']: agg[(t[1], t[2], t[3], t[4])] += t[0]
    for K in need:
        for n, N in counts[K].items():
            if n < nmin or not N: continue
            for k, mult in Counter(K).items():
                Ek = OC.E_of(K, k)
                if lo < n + Ek <= hi: agg[(n, float(Ek), k, len(K) - 1)] += mult*N
    terms = [(c, n, E, k, r) for (n, E, k, r), c in sorted(agg.items())]
    inS = lambda n, E: n + E <= lo or (n + E <= hi and n >= nmin)
    EXP = [(n, Fraction(e, 2)) for n in range(2, 2*m + 6) for e in [4] + list(range(6, 4*m + 10)) if not (n == 2 and e < 6)]
    out = [(n, E) for n, E in EXP if not inS(n, E)]
    corners = sorted((n, E) for n, E in out if not any((n2, E2) != (n, E) and n2 <= n and E2 <= E for n2, E2 in out))
    return terms, [(n, float(E)) for n, E in corners], sorted(need.items())

if __name__ == "__main__":
    base = json.load(open(sys.argv[1])); cache = json.load(open(sys.argv[2])); m = int(sys.argv[3]); nmin = int(sys.argv[4])
    terms, corners, need = assemble(base, cache, m, nmin)
    json.dump({'terms': terms, 'corners': corners}, open(sys.argv[5], 'w'))
    print(f"classes of order {m+1} used: {need}"); print(len(terms), 'terms; corners', corners)

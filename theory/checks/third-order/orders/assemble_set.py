# Section 16 of Part IV (Lemma mg:lem:corners): the exact terms and the corners of a set S of exponents, from cached class counts.
# Usage: python3 orders/assemble_set.py --base M0 [--layer M:NMIN ...] [--extra N:E,...] [--exclude N:E,...] --cache C1.json [--cache C2.json ...] --out OUT.json
"""S = S(M0) u {(n,E): n + E <= M + 1/2, n >= NMIN for each layer M:NMIN} u {the extra exponents}, minus the excluded exponents (of the layers), inside the exponent set of Lemma mg:lem:corners
(n >= 2, E in {2, 3, 7/2, 4, ...}, E >= 3 if n = 2).  The classes (sizes K of the further fronts, |K| >= 2) are those of order_counts.classes(Mmax),
Mmax the largest order used; a class with a pair (n, E) in S must be counted to that n in one of the caches (format of order_counts_par.py:
"[k1, k2, ...]|nmax" -> {n: count}; for each class the entry with the largest nmax is used), which the script asserts.  Each pair (J, p') with
(n, E) in S gives a term (count, n, E, |p'|, |K| - 1), as order_counts.py; the corners are the minimal exponents outside S.
With --base 8 alone it rebuilds orders/d4_m8.json; with --base 8 --layer 9:4 it rebuilds orders/d4_s9.json."""
import json, os, sys
from collections import Counter
from fractions import Fraction
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..'))
import order_counts as OC

def parse(argv):
    a = dict(base=None, layers=[], extra=set(), exclude=set(), caches=[], out=None)
    i = 0
    while i < len(argv):
        k = argv[i]; v = argv[i + 1]; i += 2
        if k == '--base': a['base'] = int(v)
        elif k == '--layer': m, nmin = v.split(':'); a['layers'].append((int(m), int(nmin)))
        elif k == '--extra': a['extra'] |= {(int(x.split(':')[0]), Fraction(x.split(':')[1])) for x in v.split(',')}
        elif k == '--exclude': a['exclude'] |= {(int(x.split(':')[0]), Fraction(x.split(':')[1])) for x in v.split(',')}
        elif k == '--cache': a['caches'].append(v)
        elif k == '--out': a['out'] = v
        else: raise SystemExit(f'unknown option {k}')
    return a

def in_exponent_set(n, E): return n >= 2 and E >= 2 and (n > 2 or E >= 3)

def assemble(base, layers, extra, caches, exclude=frozenset()):
    assert all(n + E > Fraction(2 * base + 1, 2) for n, E in exclude), 'only exponents of the layers can be excluded'
    def inS(n, E):
        if not in_exponent_set(n, E): return False
        if n + E <= Fraction(2 * base + 1, 2): return True
        if (n, E) in exclude: return False
        if any(n + E <= Fraction(2 * m + 1, 2) and n >= nmin for m, nmin in layers): return True
        return (n, E) in extra
    mmax = max([base] + [m for m, _ in layers] + [int(n + E + Fraction(1, 2)) for n, E in extra])
    counts, counted_to = {}, {}
    for path in caches:
        for key, cnt in json.load(open(path)).items():
            K = tuple(json.loads(key.split('|')[0])); nm = int(key.split('|')[1])
            if nm >= counted_to.get(K, -1): counts[K] = {int(n): v for n, v in cnt.items()}; counted_to[K] = nm
    need = {}
    for K, nmax in OC.classes(mmax).items():
        ns = [n for n in range(2, nmax + 1) for k in set(K) if inS(n, OC.E_of(K, k))]
        if ns: need[K] = max(ns)
    missing = [(K, n) for K, n in need.items() if counted_to.get(K, -1) < n]
    assert not missing, f'classes not counted far enough: {missing}'
    agg = Counter()
    for K in need:
        for n, N in counts[K].items():
            if not N: continue
            for k, mult in Counter(K).items():
                Ek = OC.E_of(K, k)
                if inS(n, Ek): agg[(n, float(Ek), k, len(K) - 1)] += mult * N
    terms = [(c, n, E, k, r) for (n, E, k, r), c in sorted(agg.items())]
    EXP = [(n, Fraction(e, 2)) for n in range(2, 2 * mmax + 6) for e in [4] + list(range(6, 4 * mmax + 10)) if in_exponent_set(n, Fraction(e, 2))]
    out = [(n, E) for n, E in EXP if not inS(n, E)]
    corners = sorted((n, E) for n, E in out if not any((n2, E2) != (n, E) and n2 <= n and E2 <= E for n2, E2 in out))
    return terms, [(n, float(E)) for n, E in corners], sorted(need.items())

if __name__ == "__main__":
    a = parse(sys.argv[1:])
    terms, corners, need = assemble(a['base'], a['layers'], a['extra'], a['caches'], a['exclude'])
    json.dump({'terms': terms, 'corners': corners}, open(a['out'], 'w'))
    print(f"{len(need)} classes used; {len(terms)} terms; corners {corners}")

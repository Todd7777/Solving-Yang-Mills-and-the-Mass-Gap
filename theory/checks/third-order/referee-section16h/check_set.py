# Referee (Section 16h): rebuild S^x_10 u {(9,2)} (Remark mg:rem:su5(iii)) from its definition, without order_counts.py / assemble_set.py:
# the exponent set, the classes with a pair in it (from E = (k-4) + (sum K - k)/2), the terms from the cached class counts, and the corners;
# compare with orders/d4_s10x_92.json and orders/d4_s10x.json.  Usage: python3 check_set.py
import json, os, itertools
from fractions import Fraction as F
from collections import Counter
here = os.path.dirname(os.path.abspath(__file__)); OD = os.path.join(here, '..', 'orders')
H = F(1, 2)

def in_calE(n, E):        # the exponents of Lemma mg:lem:corners: n >= 2, E in {2, 3, 7/2, 4, ...}, E >= 3 if n = 2
    return n >= 2 and (E == 2 or (E >= 3 and (2*E).denominator == 1)) and (n > 2 or E >= 3)

def S_x10(n, E, with92):
    if not in_calE(n, E): return False
    if n + E <= 8 + H: return True                                    # S(8)
    if n + E <= 9 + H and n >= 4: return True                         # order 9 with n >= 4
    if (n, E) == (3, 6): return True
    if n + E <= 10 + H and n >= 5 and (n, E) != (5, F(11, 2)): return True
    return with92 and (n, E) == (9, 2)

def E_of(K, k): return F(k - 4) + F(sum(K) - k, 2)

def classes_needed(inS, nmax_scan=30):
    sizes = [4] + list(range(6, 25))
    need = {}
    for r in range(2, 7):                                   # |K| <= 6 suffices: E >= 2(|K|-1) and every exponent of S has E <= 9
        for K in itertools.combinations_with_replacement(sizes, r):
            if min(E_of(K, k) for k in K) > 9: continue
            ns = [n for n in range(2, nmax_scan) for k in set(K) if inS(n, E_of(K, k))]
            if ns: need[K] = max(ns)
    return need

def load_caches():
    counts, upto = {}, {}
    for f in ('d4_m8_classes.json', 'd4_m9_classes.json', 'd4_m10_classes.json', 'd4_m11_classes.json'):
        for key, cnt in json.load(open(os.path.join(OD, f))).items():
            K = tuple(json.loads(key.split('|')[0])); nm = int(key.split('|')[1])
            c = {int(n): v for n, v in cnt.items()}
            if K in counts:                                     # consistency of overlapping caches
                for n in set(c) & set(counts[K]): assert c[n] == counts[K][n], (K, n)
            if nm >= upto.get(K, -1): counts[K] = c; upto[K] = nm
    return counts, upto

def build(with92):
    inS = lambda n, E: S_x10(n, E, with92)
    need = classes_needed(inS)
    counts, upto = load_caches()
    missing = {K: n for K, n in need.items() if upto.get(K, -1) < n}
    agg = Counter()
    for K, nmax in need.items():
        for n in range(2, nmax + 1):
            N = counts.get(K, {}).get(n, 0)
            for k, mult in Counter(K).items():
                if inS(n, E_of(K, k)) and N: agg[(n, E_of(K, k), k, len(K) - 1)] += mult*N
    terms = sorted((c, n, float(E), k, r) for (n, E, k, r), c in agg.items())
    EXP = [(n, F(e, 2)) for n in range(2, 30) for e in range(4, 40) if in_calE(n, F(e, 2))]
    out = [x for x in EXP if not inS(*x)]
    corners = sorted((n, float(E)) for n, E in out if not any((m, D) != (n, E) and m <= n and D <= E for m, D in out))
    return need, missing, terms, corners

if __name__ == "__main__":
    for with92, fname in ((False, 'd4_s10x.json'), (True, 'd4_s10x_92.json')):
        need, missing, terms, corners = build(with92)
        J = json.load(open(os.path.join(OD, fname)))
        jt = sorted(tuple(t) for t in J['terms']); jc = sorted(tuple(c) for c in J['corners'])
        print(f"{'S^x_10 u {(9,2)}' if with92 else 'S^x_10'}: {len(need)} classes with a pair in S; not counted far enough: {missing if missing else 'none'}")
        print(f"   {len(terms)} terms (file {fname}: {len(jt)}), identical: {terms == [tuple(t) for t in jt]};  corners {corners} (file: {jc}), identical: {corners == jc}")
        e2 = sorted({K for K in need if any(E_of(K, k) == 2 for k in K)})
        print(f"   classes with a pair of exponent E = 2: {e2}; needed to n = {[need[K] for K in e2]}")
    # E = 2 forces K = (4,4): over all multisets with sizes in {4,6,7,...,39}, |K| <= 6
    sizes = [4] + list(range(6, 25))
    allK = [K for r in range(2, 7) for K in itertools.combinations_with_replacement(sizes, r) if any(E_of(K, k) == 2 for k in K)]
    print("multisets K (|K| <= 6, sizes <= 24) with some E(K,k) = 2:", allK)

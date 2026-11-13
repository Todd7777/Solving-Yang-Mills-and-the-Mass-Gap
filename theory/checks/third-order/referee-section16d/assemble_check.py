# Referee (Section 16, fourth pass): S(8) and S''_9 rebuilt from the referee's own counts and class enumeration, compared with orders/d3_m8.json
# and orders/d3_s9.json.  Usage: python3 assemble_check.py
"""Classes: a class is (n, K), K the multiset of the sizes of the further front-pairs (|K| >= 2: at least three front-pairs), sizes 4 or >= 6 (no
fronts of five edges). A pair (J, p') with |p'| = k in K has exponent E = (k - 4) + (sum K - k)/2 (Lemma mg:lem:jsums4: E collects the factors
a^{|psi|/2} of the further front-pairs other than p', and f(k) = O(a^{k-4})), weight count * b^n a^E e^{(|K|-1) c_J} (times f(k) a^{4-k} if k > 4),
and multiplicity = the number of front-pairs of size k among the further ones. (n, E) must lie in the exponent set of Lemma mg:lem:corners: n >= 2,
and E >= 3 if n = 2.  The classes are enumerated here directly (not with order_counts.classes); the counts are read from the referee's outputs
(outputs/junctions_*.txt, by junc_faces.c).  Printed: the classes and terms of S(8) and of S''_9, and whether they equal the authors' JSON files."""
import os, re, json, itertools
from collections import Counter, defaultdict
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__))

def my_counts():
    cnt = {}
    for fn in ('junctions_validation_order8.txt', 'junctions_order9.txt'):
        for line in open(os.path.join(here, 'outputs', fn)):
            m = re.match(r'n=(\d+) R=(\d+) sizes=([\d,]+): (\d+) junctions', line)
            if m:
                n, K, c = int(m.group(1)), tuple(int(x) for x in m.group(3).split(',')), int(m.group(4))
                assert cnt.get((n, K), c) == c
                cnt[(n, K)] = c
    return cnt

def E_of(K, k): return Fr(k - 4) + Fr(sum(K) - k, 2)
def in_exp(n, E): return n >= 2 and (n > 2 or E >= 3)

def classes(inS, nmax=12, kmax=24, rmax=6):
    """all (n, K) with some pair (n, E) in S (inS a predicate), with the k's giving such pairs"""
    sizes = [4] + list(range(6, kmax + 1))
    out = {}
    for r in range(2, rmax + 1):
        for K in itertools.combinations_with_replacement(sizes, r):
            if E_of(K, K[0]) > 12: continue
            for n in range(2, nmax + 1):
                ks = sorted(k for k in set(K) if in_exp(n, E_of(K, k)) and inS(n, E_of(K, k)))
                if ks: out[(n, K)] = ks
    return out

def terms(cls, cnt):
    agg = Counter(); missing = []
    for (n, K), ks in cls.items():
        if (n, K) not in cnt: missing.append((n, K)); continue
        for k in ks:
            if cnt[(n, K)]: agg[(n, float(E_of(K, k)), k, len(K) - 1)] += Counter(K)[k]*cnt[(n, K)]
    return sorted([c, n, E, k, m] for (n, E, k, m), c in agg.items()), missing

if __name__ == "__main__":
    cnt = my_counts()
    S8 = lambda n, E: n + E <= Fr(17, 2)
    S9pp = lambda n, E: n + E <= Fr(17, 2) or (n + E <= Fr(19, 2) and n >= 4)
    c8 = classes(S8)
    print(f"S(8): {len({K for (n, K) in c8})} classes (multisets K), {len(c8)} (n, K) pairs of class and size")
    t8, miss8 = terms(c8, cnt)
    print(f"  terms: {len(t8)} kinds of pairs; (n, K) not counted by the referee: {miss8}")
    J8 = json.load(open(os.path.join(here, '..', 'orders', 'd3_m8.json')))
    print("  equal to orders/d3_m8.json:", sorted(map(list, J8['terms'])) == t8)
    c9 = classes(S9pp)
    new = {nk: ks for nk, ks in c9.items() if nk not in c8 or any(not S8(nk[0], E_of(nk[1], k)) for k in ks)}
    print(f"S''_9: classes with a pair outside S(8): {sorted({K for (n, K) in new})} ({len({K for (n, K) in new})} multisets)")
    t9, miss9 = terms(c9, cnt)
    print(f"  terms: {len(t9)}; (n, K) not counted by the referee: {miss9}")
    J9 = json.load(open(os.path.join(here, '..', 'orders', 'd3_s9.json')))
    print("  equal to orders/d3_s9.json:", sorted(map(list, J9['terms'])) == t9)
    if sorted(map(list, J9['terms'])) != t9:
        A = {tuple(x[1:]): x[0] for x in J9['terms']}; B = {tuple(x[1:]): x[0] for x in t9}
        for k in sorted(set(A) | set(B)):
            if A.get(k) != B.get(k): print('   differ at', k, A.get(k), B.get(k))
    # all order-9 classes (m = 9, any n), for the count "36 classes"
    c9all = classes(lambda n, E: n + E <= Fr(19, 2))
    print(f"S(9): {len({K for (n, K) in c9all})} classes (multisets K)")

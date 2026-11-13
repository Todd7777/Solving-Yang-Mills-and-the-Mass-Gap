# Section 16 of Part IV (Remark mg:rem:nextorder (vii)): the set S_9'' in d=3 and its corners.  Usage: python3 orders/assemble_d3_s9.py
"""S_9'' = S(8) u {pairs of the layer 17/2 < n + E <= 19/2 with n >= 4}; the pairs of that layer with n = 2, 3 (classes with fronts of 12 to 15
edges, slow to count) are left out, which Lemma mg:lem:corners allows: the corners become (2,7), (3,6), (5,5), (6,4), (7,3), (8,2).
The counts of the classes used, by n (order_counts.count_junctions in d=3), are listed below: the first eight from order_counts_par.py, the last
five counted separately; S(8) is orders/d3_m8.json (order_counts.py 2 8).  Output: orders/d3_s9.json."""
import json, os, sys
from collections import Counter
from fractions import Fraction
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..'))
import order_counts as OC
CLASSES = {
 "[4, 4]|7": {
  "3": 48
 },
 "[4, 6]|6": {
  "2": 24,
  "3": 24,
  "4": 276
 },
 "[4, 7]|6": {
  "2": 24,
  "3": 24,
  "4": 276
 },
 "[4, 8]|5": {
  "1": 4,
  "2": 48,
  "3": 252,
  "4": 604,
  "5": 1152
 },
 "[4, 9]|5": {
  "1": 8,
  "2": 16,
  "3": 424,
  "4": 704,
  "5": 2288
 },
 "[4, 10]|4": {
  "1": 20,
  "2": 100,
  "3": 732,
  "4": 2432
 },
 "[4, 11]|4": {
  "1": 16,
  "2": 212,
  "3": 856,
  "4": 4556
 },
 "[4, 12]|3": {
  "1": 12,
  "2": 380,
  "3": 1872
 },
 "[4, 4, 7]|4": {
  "3": 168,
  "4": 168
 },
 "[4, 4, 6]|4": {
  "3": 168,
  "4": 168
 },
 "[6, 7]|4": {
  "3": 128,
  "4": 128
 },
 "[6, 6]|4": {
  "3": 64,
  "4": 64
 },
 "[4, 4, 4]|5": {
  "4": 292,
  "5": 4
 }
}
counts = {}; counted_to = {}
for key, cnt in CLASSES.items():
    K = tuple(json.loads(key.split('|')[0])); counts[K] = {int(n): v for n, v in cnt.items()}; counted_to[K] = int(key.split('|')[1])
lo, hi = Fraction(17, 2), Fraction(19, 2)
# the classes of order 9 (order_counts.classes) with a pair in the layer 17/2 < n + E <= 19/2 with n >= 4, and the largest such n: each must be counted that far
need = {}
for K, nmax in OC.classes(9).items():
    ns = [n for n in range(4, nmax + 1) for k in set(K) if lo < n + OC.E_of(K, k) <= hi]
    if ns: need[K] = max(ns)
missing = [K for K in need if K not in counts or counted_to[K] < need[K]]
assert not missing, missing
m8 = json.load(open(os.path.join(here, 'd3_m8.json')))
agg = Counter()
for t in m8['terms']: agg[(t[1], t[2], t[3], t[4])] += t[0]
for K in need:
    for n, N in counts[K].items():
        if n < 4 or not N: continue
        for k, mult in Counter(K).items():
            Ek = OC.E_of(K, k)
            if lo < n + Ek <= hi: agg[(n, float(Ek), k, len(K) - 1)] += mult * N
terms = [(c, n, E, k, r) for (n, E, k, r), c in sorted(agg.items())]
corners = [(2, 7), (3, 6), (5, 5), (6, 4), (7, 3), (8, 2)]
# check: these are the minimal elements of the complement of S_9'' in the exponent set (n >= 2; E in {2, 3, 7/2, 4, ...}, E >= 3 if n = 2)
EXP = [(n, Fraction(e, 2)) for n in range(2, 12) for e in [4] + list(range(6, 24)) if not (n == 2 and e < 6)]
inS = lambda n, E: n + E <= lo or (n + E <= hi and n >= 4)
out = [(n, E) for n, E in EXP if not inS(n, E)]
minimal = sorted((n, E) for n, E in out if not any((n2, E2) != (n, E) and n2 <= n and E2 <= E for n2, E2 in out))
assert minimal == [(n, Fraction(E)) for n, E in corners], minimal
json.dump({'terms': terms, 'corners': corners}, open(os.path.join(here, 'd3_s9.json'), 'w'))
print(len(terms), 'terms; corners', corners)

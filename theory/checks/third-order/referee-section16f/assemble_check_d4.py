# Referee (Section 16, sixth pass): S(7) and S''_8 in d=4 rebuilt from the referee's own counts (junc_partition.c), compared with
# ../orders/d4_m7.json, ../orders/d4_s8.json and ../orders/d4_s8b.json; the classes are enumerated here directly from the definition of E, and the corners of S''_8
# are computed as the minimal elements of the exponent set outside it.          Usage: python3 assemble_check_d4.py
"""For a junction with n plaquettes and further front-pairs of sizes K (|K| >= 2), each choice of p' in K gives a pair (J, p') with
E = (|p'| - 4) + (sum K - |p'|)/2 and weight count * b^n a^E e^{(|K|-1) c_J} (f(|p'|) a^{4-|p'|} if |p'| > 4) (proof of Lemma mg:lem:jsums4).
Exponents: n >= 2, E in {2, 3, 7/2, 4, ...}, E >= 3 if n = 2 (Lemma mg:lem:corners); no fronts with five edges."""
import json, os, itertools
from fractions import Fraction as Fr
from collections import Counter
here = os.path.dirname(os.path.abspath(__file__))

# the referee's counts (outputs/junc_partition_d4.txt, junc_partition_d4_order8_n3.txt): (sorted K) -> {n: junctions}
MY = {(4, 4): {3: 876, 4: 42, 5: 912, 6: 1866}, (4, 6): {2: 320, 3: 680, 4: 18276, 5: 11744}, (4, 7): {2: 312, 3: 344, 4: 17028, 5: 3440},
      (4, 8): {2: 992, 3: 9220, 4: 56960}, (4, 9): {2: 256, 3: 15576, 4: 50920}, (4, 10): {2: 2260, 3: 45700}, (4, 11): {2: 4636, 3: 51640},
      (6, 6): {2: 64, 3: 3028}, (6, 7): {2: 0, 3: 5712}, (4, 4, 4): {2: 4, 3: 16, 4: 24314}, (4, 4, 6): {2: 0, 3: 9768}, (4, 4, 7): {2: 0, 3: 9600}}
COUNTED_TO = {(4, 4): 6, (4, 6): 5, (4, 7): 5, (4, 8): 4, (4, 9): 4, (4, 10): 3, (4, 11): 3, (6, 6): 3, (6, 7): 3, (4, 4, 4): 4, (4, 4, 6): 3, (4, 4, 7): 3}

def Eof(K, k): return Fr(k - 4) + Fr(sum(K) - k, 2)
def in_exponent_set(n, E): return n >= 2 and E >= 2 and (n > 2 or E >= 3)

def pairs_needed(inS, maxsize=30, maxfronts=6, nmax=12):
    """all (K, n, k) with (n, E) in S, K a multiset of sizes in {4, 6, 7, ...}; S must be contained in {n + E <= 9}"""
    sizes = [4] + list(range(6, maxsize + 1)); out = []
    for r in range(2, maxfronts + 1):
        for K in itertools.combinations_with_replacement(sizes, r):
            if Eof(K, K[0]) > 9: continue
            for k in sorted(set(K)):
                E = Eof(K, k)
                for n in range(2, nmax + 1):
                    if in_exponent_set(n, E) and inS(n, E): out.append((K, n, k))
    return out

def terms_from(needed):
    agg = Counter(); missing = []
    for K, n, k in needed:
        if K not in MY or COUNTED_TO[K] < n: missing.append((K, n, k)); continue
        N = MY[K].get(n, 0)
        if N: agg[(n, float(Eof(K, k)), k, len(K) - 1)] += Counter(K)[k] * N
    return sorted([[c, n, E, k, r] for (n, E, k, r), c in agg.items()], key=lambda t: (t[1], t[2], t[3], t[4])), missing

def corners(inS):
    EXP = [(n, Fr(e, 2)) for n in range(2, 14) for e in [4] + list(range(6, 30)) if in_exponent_set(n, Fr(e, 2))]
    out = [(n, E) for n, E in EXP if not inS(n, E)]
    return sorted((n, E) for n, E in out if not any((n2, E2) != (n, E) and n2 <= n and E2 <= E for n2, E2 in out))

S7 = lambda n, E: n + E <= Fr(15, 2)
S8pp = lambda n, E: n + E <= Fr(15, 2) or (n + E <= Fr(17, 2) and n >= 4)
S8b = lambda n, E: S8pp(n, E) or (n, E) == (3, 5)        # the revised (ix): S''_8 with (3,5) added
for name, inS, path in (("S(7)", S7, "d4_m7.json"), ("S''_8", S8pp, "d4_s8.json"), ("S''_8 u {(3,5)}", S8b, "d4_s8b.json")):
    need = pairs_needed(inS)
    terms, missing = terms_from(need)
    J = json.load(open(os.path.join(here, '..', 'orders', path)))
    theirs = sorted([list(t) for t in J['terms']], key=lambda t: (t[1], t[2], t[3], t[4]))
    print(f"{name}: classes with a pair in the set: {sorted(set(K for K, n, k in need))}")
    print(f"  pairs not covered by the referee's counts: {missing}")
    print(f"  {len(terms)} terms from the referee's counts; equal to orders/{path}: {terms == theirs}")
    if terms != theirs:
        print("   ours  :", terms); print("   theirs:", theirs)
    cs = corners(inS)
    print(f"  corners (minimal exponents outside the set): {[(n, float(E)) for n, E in cs]};  file: {J['corners']}")

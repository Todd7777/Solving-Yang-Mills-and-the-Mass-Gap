# Referee check for Part IV, Section 16, Remark mg:rem:practical(ii): N^v_k, the number of fronts with k edges through a vertex of Z^n
# (a front: a non-empty connected set of edges in which every vertex has degree >= 2; Section 7 of Part IV). Written independently of
# fronts_count.py (which uses Redelmeier's method); here a decision search:
#   * first decide the 2n edges at the origin v: a front through v contains a subset of >= 2 of them (the others are excluded);
#   * while some vertex has degree 1, take the one with fewest undecided edges and branch over the non-empty subsets of its undecided
#     edges (the rest excluded);
#   * when every vertex has degree >= 2, branch on the smallest undecided edge sharing a vertex with the current set (include / exclude);
#     the set is recorded at the end of this chain.
# Every front through v is recorded exactly once (same argument as in ref_closed.py).
# Usage: python3 ref_fronts.py n KMAX
import sys, itertools
from collections import Counter
n = int(sys.argv[1]); KMAX = int(sys.argv[2])
U = [tuple(1 if k == i else 0 for k in range(n)) for i in range(n)]
def add(x, y): return tuple(a + b for a, b in zip(x, y))
def sub(x, y): return tuple(a - b for a, b in zip(x, y))
def ends(e): x, k = e; return (x, add(x, U[k]))
def at(v): return [(v, k) for k in range(n)] + [(sub(v, U[k]), k) for k in range(n)]
found = Counter(); nodes = 0
cycles = Counter()     # how many of the fronts are cycles (all degrees 2), by size

def rec(IN, OUT, deg):
    global nodes
    nodes += 1
    defi = [w for w, c in deg.items() if c == 1]
    if len(IN) + (len(defi) + 1) // 2 > KMAX: return
    if defi:
        best = None
        for w in defi:
            c = [e for e in at(w) if e not in IN and e not in OUT]
            if best is None or len(c) < len(best[1]):
                best = (w, c)
                if len(c) <= 1: break
        w, cands = best
        room = KMAX - len(IN)
        for r in range(1, min(len(cands), room) + 1):
            for T in itertools.combinations(cands, r):
                deg2 = Counter(deg)
                for e in T:
                    for y in ends(e): deg2[y] += 1
                rec(IN | frozenset(T), OUT | (frozenset(cands) - frozenset(T)), deg2)
        return
    if len(IN) < KMAX:
        und = sorted({e for w in deg for e in at(w) if e not in IN and e not in OUT})
        if und:
            c = und[0]; deg2 = Counter(deg)
            for y in ends(c): deg2[y] += 1
            rec(IN | {c}, OUT, deg2)
            rec(IN, OUT | {c}, deg)
            return
    found[len(IN)] += 1
    if all(c == 2 for c in deg.values()): cycles[len(IN)] += 1

v = tuple([0] * n); E0 = at(v)
for r in range(2, len(E0) + 1):
    for S in itertools.combinations(E0, r):
        deg = Counter()
        for e in S:
            for y in ends(e): deg[y] += 1
        rec(frozenset(S), frozenset(E0) - frozenset(S), deg)
print(f"Z^{n}: fronts through a vertex by number of edges (k <= {KMAX}): {dict(sorted(found.items()))}; of which cycles: {dict(sorted(cycles.items()))}; nodes {nodes}")

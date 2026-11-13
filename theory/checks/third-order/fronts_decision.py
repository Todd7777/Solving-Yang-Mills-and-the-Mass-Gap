# Section 16 of Part IV (Remark mg:rem:practical (ii)): exact N^v_k by a decision search.  Usage: python3 fronts_decision.py D K
"""Fronts through the origin of Z^D with at most K edges (connected edge sets, every vertex of degree >= 2), by a decision search: decide all edges at one vertex at a time (first the origin, with at
least two edges; then always an undecided vertex of the current set, preferring degree one, then fewest undecided edges, ties by the vertex). For a
given front the sequence of decisions is forced, so every front is produced exactly once. Links are (x, i): the edge from x to x + e_i.
Results (agree with fronts_count.py, a Redelmeier enumeration, for k <= 9 in Z^2 and k <= 8 in Z^3, and with the referee's ref_fronts.py for k <= 10, 9):
  D=2, K=12 (9 s):  {4: 4, 6: 12, 7: 12, 8: 70, 9: 160, 10: 634, 11: 1728, 12: 6369}
  D=3, K=10 (2.5 min): {4: 12, 6: 132, 7: 108, 8: 2034, 9: 4344, 10: 43650}"""
import sys, itertools
from collections import Counter
def make(D):
    E = [tuple(1 if k == i else 0 for k in range(D)) for i in range(D)]
    add = lambda x, y: tuple(a + b for a, b in zip(x, y))
    sub = lambda x, y: tuple(a - b for a, b in zip(x, y))
    inc = lambda v: [(v, i) for i in range(D)] + [(sub(v, E[i]), i) for i in range(D)]
    ends = lambda l: (l[0], add(l[0], E[l[1]]))
    return E, add, sub, inc, ends
def fronts_through_origin(D, K):
    E, add, sub, inc, ends = make(D)
    O = (0,) * D
    out = []
    def rec(S, X, deg, decided, v):
        U = [l for l in inc(v) if l not in S and l not in X]
        for r in range(0, len(U) + 1):
            if len(S) + r > K: break
            if deg[v] + r < 2: continue
            for T in itertools.combinations(U, r):
                S2 = S | frozenset(T); X2 = X | frozenset(l for l in U if l not in T)
                deg2 = Counter(deg)
                for l in T:
                    for w in ends(l): deg2[w] += 1
                dec2 = decided | {v}
                cand = [w for w, c in deg2.items() if c >= 1 and w not in dec2]
                if not cand:
                    out.append(S2); continue
                ones = sum(1 for w in cand if deg2[w] == 1)
                if len(S2) + (ones + 1) // 2 > K: continue
                def key(w):
                    return (0 if deg2[w] == 1 else 1, sum(1 for l in inc(w) if l not in S2 and l not in X2), w)
                rec(S2, X2, deg2, dec2, min(cand, key=key))
    rec(frozenset(), frozenset(), Counter(), frozenset(), O)
    return out
if __name__ == "__main__":
    D, K = int(sys.argv[1]), int(sys.argv[2])
    F = fronts_through_origin(D, K)
    assert len(F) == len(set(F))
    print(f"D={D}, K={K}: {dict(sorted(Counter(len(f) for f in F).items()))}")

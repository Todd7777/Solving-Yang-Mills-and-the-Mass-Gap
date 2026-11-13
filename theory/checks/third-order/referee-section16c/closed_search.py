# Referee (Section 16, third pass): closed polymers through a plaquette.  Usage: python3 closed_search.py n NMAX
"""All link-connected sets Y of plaquettes of Z^n containing p0, |Y| <= NMAX, in which every link of E(Y) lies in at least two plaquettes
of Y (the polymers with non-zero activity, proof of Lemma mg:lem:closed).
Search (independent of closed_sets.py, closed_sets_ext.py and ref_closed.py): one plaquette per step.  If some link of the current set S
lies in exactly one plaquette, Y contains a further plaquette through the smallest such link: branch over all of them.  If S is closed and
S != Y, Y (link-connected) contains a plaquette sharing a link with S: branch over all of them.  There are no exclusion sets, so the result
of the search from S depends on S only, and memoising visited sets is sound.  Prune: each link lying in one plaquette of S needs a new
plaquette through it, and a new plaquette has 4 links."""
import sys
from collections import Counter
n = int(sys.argv[1]); NMAX = int(sys.argv[2])
E = [tuple(1 if k == i else 0 for k in range(n)) for i in range(n)]
add = lambda x, y: tuple(a + b for a, b in zip(x, y))
sub = lambda x, y: tuple(a - b for a, b in zip(x, y))
def links(p):
    x, i, j = p
    return ((x, i), (add(x, E[i]), j), (add(x, E[j]), i), (x, j))
def plaqs_through(l):
    x, k = l; out = []
    for m in range(n):
        if m == k: continue
        i, j = min(k, m), max(k, m)
        out.append((x, i, j)); out.append((sub(x, E[m]), i, j))
    return out
p0 = ((0,) * n, 0, 1)
found = set(); visited = set()
sys.setrecursionlimit(100000)
def rec(S, cnt):
    if S in visited: return
    visited.add(S)
    defi = [l for l, c in cnt.items() if c == 1]
    if len(defi) > 4 * (NMAX - len(S)): return
    if not defi:
        found.add(S)
        if len(S) == NMAX: return
        cands = {p for q in S for l in links(q) for p in plaqs_through(l)} - S
    else:
        cands = [p for p in plaqs_through(min(defi)) if p not in S]
    for p in cands:
        cnt2 = Counter(cnt)
        for l in links(p): cnt2[l] += 1
        rec(S | {p}, cnt2)
rec(frozenset([p0]), Counter(links(p0)))
by = Counter(len(Y) for Y in found)
print(f"Z^{n}, NMAX={NMAX}: closed link-connected sets through p0 by size: {dict(sorted(by.items()))}  ({len(visited)} states)")
# the 2(n-2) unit cubes through p0 are the only closed sets of at most 9 plaquettes iff by[k] == 0 for 7 <= k <= 9 and by[6] == 2(n-2)

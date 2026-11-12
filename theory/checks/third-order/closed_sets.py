# Section 16 of Part IV (Lemma mg:lem:seven and Remark mg:rem:practical): a compass, not part of the proofs. Usage: python3 closed_sets.py D NMAX
"""Enumerate connected sets Y of plaquettes of Z^D containing a fixed plaquette p0, |Y| <= NMAX, in which every link of E(Y) lies in
at least two plaquettes of Y ("closed"). Deficiency-driven search: take the smallest link lying in exactly one plaquette of the current
set S, and branch over the non-empty sets of further plaquettes through it (the others through it are excluded for good). Every closed Y
containing p0 is produced exactly once, provided no closed proper subset S of Y can be completed (true for NMAX <= 9: a closed S has >= 6
plaquettes, and <= 3 further plaquettes cannot close their own new links)."""
import sys, itertools
from collections import Counter
D = int(sys.argv[1]); NMAX = int(sys.argv[2])
E = [tuple(1 if k == i else 0 for k in range(D)) for i in range(D)]
def add(x, y): return tuple(a + b for a, b in zip(x, y))
def links(p):
    x, i, j = p
    return ((x, i), (add(x, E[i]), j), (add(x, E[j]), i), (x, j))
def plaqs_through(l):
    x, k = l; out = []
    for m in range(D):
        if m == k: continue
        i, j = min(k, m), max(k, m)
        out.append((x, i, j))                         # the plaquette spanned by k, m at x
        out.append((add(x, tuple(-c for c in E[m])), i, j))  # and at x - e_m
    return out
p0 = (tuple([0] * D), 0, 1)
found = Counter(); shapes = {}
nodes = 0
def rec(S, cnt, excl):
    global nodes
    nodes += 1
    defi = [l for l, c in cnt.items() if c == 1]
    if not defi:
        found[len(S)] += 1
        key = len(S)
        shapes.setdefault(key, []).append(frozenset(S)); return
    l = min(defi)
    cands = [p for p in plaqs_through(l) if p not in S and p not in excl]
    room = NMAX - len(S)
    for r in range(1, min(len(cands), room) + 1):
        for T in itertools.combinations(cands, r):
            S2 = S | set(T); cnt2 = Counter(cnt)
            for p in T:
                for ll in links(p): cnt2[ll] += 1
            excl2 = excl | (set(cands) - set(T))
            # every plaquette through l is now decided
            rec(S2, cnt2, frozenset(excl2))
cnt0 = Counter(links(p0))
rec({p0}, cnt0, frozenset())
print(f"D={D}, NMAX={NMAX}: closed sets containing p0 by size: {dict(sorted(found.items()))}; search nodes {nodes}")

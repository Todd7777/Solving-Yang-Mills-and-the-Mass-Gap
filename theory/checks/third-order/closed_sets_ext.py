# Section 16 of Part IV (Remark mg:rem:nextorder): closed sets of plaquettes up to NMAX, also those containing a closed proper subset.
# Usage: python3 closed_sets_ext.py D NMAX
"""Connected sets Y of plaquettes of Z^D containing a fixed plaquette p0, |Y| <= NMAX, in which every link of E(Y) lies in at least two plaquettes of
Y ("closed"). Deficiency-driven search as in closed_sets.py (at the smallest link lying in exactly one plaquette, branch over the further plaquettes
through it), and in addition, at every closed state S with |S| < NMAX, each plaquette sharing a link with S is added and the search continued: a
closed Y of at most 11 plaquettes is link-connected (two closed components would need 12), so it contains a plaquette sharing a link with any closed
proper subset. Closed sets are counted once (a set of found closed sets)."""
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
        out.append((x, i, j)); out.append((add(x, tuple(-c for c in E[m])), i, j))
    return out
p0 = (tuple([0] * D), 0, 1)
found = set(); seen_states = set()
def rec(S, cnt, excl):
    key = frozenset(S)
    if key in seen_states: return
    seen_states.add(key)
    defi = [l for l, c in cnt.items() if c == 1]
    if not defi:
        found.add(key)
        if len(S) < NMAX:
            nb = {p for q in S for l in links(q) for p in plaqs_through(l) if p not in S}
            for p in nb:
                cnt2 = Counter(cnt)
                for ll in links(p): cnt2[ll] += 1
                rec(S | {p}, cnt2, frozenset())
        return
    l = min(defi)
    cands = [p for p in plaqs_through(l) if p not in S and p not in excl]
    room = NMAX - len(S)
    for r in range(1, min(len(cands), room) + 1):
        for T in itertools.combinations(cands, r):
            cnt2 = Counter(cnt)
            for p in T:
                for ll in links(p): cnt2[ll] += 1
            rec(S | set(T), cnt2, excl | frozenset(c for c in cands if c not in T))
rec({p0}, Counter(links(p0)), frozenset())
print(f"D={D}, NMAX={NMAX}: closed sets containing p0 by size: {dict(sorted(Counter(len(Y) for Y in found).items()))}; states {len(seen_states)}")

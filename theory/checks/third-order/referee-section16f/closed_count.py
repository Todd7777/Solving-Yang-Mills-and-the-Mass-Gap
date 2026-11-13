# Referee (Section 16, sixth pass): closed polymers through a plaquette of Z^d, for theta (Lemma mg:lem:theta7, Remark mg:rem:nextorder(ii),(ix)).
# Usage: python3 closed_count.py d NMAX
"""Link-connected sets Y of plaquettes of Z^d containing the plaquette p0 = (0; e_0, e_1), in which every link of E(Y) lies in at least two
plaquettes of Y, counted by |Y| <= NMAX.  Partition search (each set reached once, checked by a set of the results): at the smallest link of
coverage one, branch i = "c_i is the first plaquette through it (in a fixed order) that Y contains", c_1..c_{i-1} excluded below; at a state
without such links, either Y = X or Y contains a plaquette sharing a link with X (Y is link-connected), branch as before.  Prune: the further
plaquettes number at most NMAX - |X| and each covers at most four links of coverage one."""
import sys, itertools
from collections import Counter
d = int(sys.argv[1]); NMAX = int(sys.argv[2])
E = [tuple(1 if k == i else 0 for k in range(d)) for i in range(d)]
add = lambda x, y: tuple(a + b for a, b in zip(x, y))
sub = lambda x, y: tuple(a - b for a, b in zip(x, y))
def plinks(c):
    x, i, j = c
    return ((x, i), (add(x, E[i]), j), (add(x, E[j]), i), (x, j))
def plaqs_through(l):
    x, k = l; out = []
    for m in range(d):
        if m == k: continue
        i, j = min(k, m), max(k, m)
        out.append((x, i, j)); out.append((sub(x, E[m]), i, j))
    return out
p0 = ((0,) * d, 0, 1)
found = Counter(); allY = set(); dups = 0
cov = Counter(plinks(p0)); X = [p0]; inX = {p0}; excl = set()
def rec():
    global dups
    defi = [l for l in set(l for c in X for l in plinks(c)) if cov[l] == 1]
    if -(-len(defi) // 4) > NMAX - len(X): return
    if not defi:
        key = frozenset(X)
        if key in allY: dups += 1
        else: allY.add(key); found[len(X)] += 1
        if len(X) == NMAX: return
        cands = []
        for l in sorted(set(l for c in X for l in plinks(c))):
            for c in plaqs_through(l):
                if c not in inX and c not in excl and c not in cands: cands.append(c)
    else:
        l = min(defi)
        if len(X) == NMAX: return
        cands = [c for c in plaqs_through(l) if c not in inX and c not in excl]
    done = []
    for c in cands:
        X.append(c); inX.add(c)
        for m in plinks(c): cov[m] += 1
        rec()
        for m in plinks(c): cov[m] -= 1
        X.pop(); inX.discard(c)
        excl.add(c); done.append(c)
    for c in done: excl.discard(c)
rec()
print(f"d={d}, n<={NMAX}: closed link-connected sets through a plaquette, by size: {dict(sorted(found.items()))}  (duplicates: {dups})")

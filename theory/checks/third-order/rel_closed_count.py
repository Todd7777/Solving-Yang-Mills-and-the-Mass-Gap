# Part IV, the upper bound at practical beta (Lemma mg:lem:twopoint): the sets Y of its proof, counted exactly.  Usage: python3 rel_closed_count.py D A|Q NMAX
"""For (S, C) = (S_A, {p0}) (S_A the four links of p0) or (S_Q, Q) (S_Q the twelve edges of the unit 3-cube Q spanned by p0 and the direction 2,
and Q its six faces), count the sets Y of plaquettes of Z^D with
  (1) Y and C disjoint, |Y| = n <= NMAX;
  (2) every link of E(Y) not in S lies in at least two plaquettes of Y ("closed relative to S");
  (3) every link-connected component of Y contains a plaquette of P_S minus C (P_S: the plaquettes with a link in S).
These are the sets Y = X minus C of the proof of Lemma mg:lem:twopoint, which bounds their number by binom(J' + Kn, n) <= (1 + 1/K)^J' (eK)^n.
Each component satisfies (2) and (3) itself and has at least five plaquettes (Lemma mg:lem:five), so for n <= 9 the sets are connected.

Search: from each plaquette r of P_S minus C, a deficiency-driven search over connected sets containing r (at the smallest link outside S that lies in
exactly one plaquette, branch over the nonempty sets of further plaquettes through it), extended at every closed state by each plaquette sharing a
link with it (a connected Y contains one with any closed proper subset), as in closed_sets_ext.py; the sets found are kept in one set, so each Y is
counted once.  Sets with two components (n >= 10) are counted from the pairs of five-plaquette components that share no link.
Also printed: the five-plaquette sets that are 'caps' (the five faces other than one face F of a unit cube, with F in C)."""
import sys, itertools
from collections import Counter
D = int(sys.argv[1]); WHICH = sys.argv[2]; NMAX = int(sys.argv[3])
E = [tuple(1 if k == i else 0 for k in range(D)) for i in range(D)]
O = tuple([0] * D)
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
p0 = (O, 0, 1)
if WHICH == 'A':
    C = {p0}
else:                                     # the unit cube spanned by directions 0, 1, 2 at the origin
    C = {(O, 0, 1), (E[2], 0, 1), (O, 0, 2), (E[1], 0, 2), (O, 1, 2), (E[0], 1, 2)}
S = {l for p in C for l in links(p)}
PS = sorted({p for l in S for p in plaqs_through(l)} - C)
print(f"D={D}, (S, C) = {'(S_A, {p0})' if WHICH == 'A' else '(S_Q, Q)'}: |S| = {len(S)}, |C| = {len(C)}, |P_S minus C| = J' = {len(PS)}")
found = set()
def rec(Y, cnt, excl):
    defi = [l for l, c in cnt.items() if c == 1 and l not in S]
    if not defi:
        key = frozenset(Y)
        if key in found: return
        found.add(key)
        if len(Y) < NMAX:
            nb = {p for q in Y for l in links(q) for p in plaqs_through(l) if p not in Y and p not in C}
            for p in nb:
                cnt2 = Counter(cnt)
                for ll in links(p): cnt2[ll] += 1
                rec(Y | {p}, cnt2, frozenset())
        return
    l = min(defi)
    cands = [p for p in plaqs_through(l) if p not in Y and p not in excl and p not in C]
    room = NMAX - len(Y)
    for r in range(1, min(len(cands), room) + 1):
        for T in itertools.combinations(cands, r):
            cnt2 = Counter(cnt)
            for p in T:
                for ll in links(p): cnt2[ll] += 1
            rec(Y | set(T), cnt2, excl | frozenset(c for c in cands if c not in T))
for r in PS:
    rec({r}, Counter(links(r)), frozenset())
conn = Counter(len(Y) for Y in found)
total = Counter(conn)
fives = [Y for Y in found if len(Y) == 5]
if NMAX >= 10:                            # two components of five plaquettes each, sharing no link
    two = 0
    for Y1, Y2 in itertools.combinations(fives, 2):
        if not ({l for p in Y1 for l in links(p)} & {l for p in Y2 for l in links(p)}): two += 1
    total[10] += two
def is_cap(Y):
    for F in C:                            # Y and F are the six faces of a unit cube
        Z = set(Y) | {F}
        cnt = Counter(l for p in Z for l in links(p))
        if all(c == 2 for c in cnt.values()) and len(cnt) == 12: return True
    return False
print(f"connected sets by size: {dict(sorted(conn.items()))}")
print(f"all sets Y by size (N_n): {dict(sorted(total.items()))}")
print(f"five-plaquette sets: {len(fives)}, of which caps (a unit cube minus a face in C): {sum(is_cap(Y) for Y in fives)}")

# Referee check for Part IV, Section 16, Lemma mg:lem:seven (written independently of closed_sets.py).
#
# Enumerates ALL link-connected sets Y of plaquettes containing a fixed plaquette p0, |Y| <= NMAX, in which every link of E(Y)
# lies in at least two plaquettes of Y ("closed"), in Z^D or on a periodic lattice (Z/L)^D.
#
# Search by decisions (each set reached exactly once): a node is (IN, OUT).
#   * If IN has a deficient link (in exactly one plaquette of IN), take the one with fewest undecided plaquettes through it and branch
#     over the non-empty subsets T of those (T joins IN, the rest joins OUT).
#   * If IN is closed, branch on the smallest undecided plaquette sharing a link with IN: include it, or exclude it. IN is recorded at the
#     end of this chain (no undecided plaquette left next to IN, or |IN| = NMAX).
# Unlike a pure deficiency search this also finds closed sets that contain a closed proper subset (e.g. two cubes sharing an edge).
# Completeness: if Y is closed and connected, IN is a subset of Y and OUT is disjoint from Y, then either IN has a deficient link (and Y
# contains one of the undecided plaquettes through it), or IN is closed and, if IN != Y, Y contains an undecided plaquette sharing a link
# with IN (Y is connected). So the branch consistent with Y always exists, and Y is recorded once.
# Usage: python3 ref_closed.py D NMAX [L]      (L: side of the periodic lattice; omit for Z^D)
import sys, itertools
from collections import Counter

D = int(sys.argv[1]); NMAX = int(sys.argv[2]); L = int(sys.argv[3]) if len(sys.argv) > 3 else None
UNIT = [tuple(1 if k == i else 0 for k in range(D)) for i in range(D)]

def norm(x):
    return tuple(c % L for c in x) if L else x
def add(x, y): return norm(tuple(a + b for a, b in zip(x, y)))
def sub(x, y): return norm(tuple(a - b for a, b in zip(x, y)))

def links(p):
    x, i, j = p
    return ((x, i), (add(x, UNIT[i]), j), (add(x, UNIT[j]), i), (x, j))

def through(l):
    x, k = l
    out = []
    for m in range(D):
        if m == k: continue
        i, j = min(k, m), max(k, m)
        out.append((x, i, j)); out.append((sub(x, UNIT[m]), i, j))
    return out

p0 = (tuple([0] * D), 0, 1)
found = Counter(); examples = {}
nodes = 0

def rec(IN, OUT, cnt):
    global nodes
    nodes += 1
    deficient = [l for l, c in cnt.items() if c == 1]
    if len(IN) + (len(deficient) + 3) // 4 > NMAX:
        return
    if deficient:
        best = None
        for l in deficient:
            c = [p for p in through(l) if p not in IN and p not in OUT]
            if best is None or len(c) < len(best[1]):
                best = (l, c)
                if len(c) <= 1: break
        l, cands = best
        room = NMAX - len(IN)
        for r in range(1, min(len(cands), room) + 1):
            for T in itertools.combinations(cands, r):
                cnt2 = Counter(cnt)
                for p in T:
                    for ll in links(p): cnt2[ll] += 1
                rec(IN | frozenset(T), OUT | (frozenset(cands) - frozenset(T)), cnt2)
        return
    # IN is closed
    if len(IN) < NMAX:
        undecided = sorted({p for l in cnt for p in through(l) if p not in IN and p not in OUT})
        if undecided:
            c = undecided[0]
            cnt2 = Counter(cnt)
            for ll in links(c): cnt2[ll] += 1
            rec(IN | {c}, OUT, cnt2)          # include c
            rec(IN, OUT | {c}, cnt)           # exclude c
            return
    found[len(IN)] += 1
    examples.setdefault(len(IN), []).append(IN)

rec(frozenset([p0]), frozenset(), Counter(links(p0)))
lat = f"(Z/{L})^{D}" if L else f"Z^{D}"
print(f"{lat}, NMAX={NMAX}: closed connected sets containing p0, by size: {dict(sorted(found.items()))}; search nodes {nodes}")

def is_cube(Y):
    # six plaquettes = boundary of a unit cube: all vertices within one unit box, 8 distinct vertices, 12 distinct links
    if len(Y) != 6: return False
    E = Counter(l for p in Y for l in links(p))
    V = {v for (v, k) in E} | {add(v, UNIT[k]) for (v, k) in E}
    return len(E) == 12 and len(V) == 8 and all(c == 2 for c in E.values())
for n, ex in sorted(examples.items()):
    if n == 6:
        print(f"  size 6: {len(ex)} sets, all unit cubes: {all(is_cube(Y) for Y in ex)}  (2(D-2) = {2*(D-2)})")
    elif n <= 10:
        shapes = Counter()
        for Y in ex:
            E = Counter(l for p in Y for l in links(p))
            shapes[(len(E), tuple(sorted(Counter(E.values()).items())))] += 1
        print(f"  size {n}: {len(ex)} sets; (number of links, link multiplicities) -> count: {dict(shapes)}")

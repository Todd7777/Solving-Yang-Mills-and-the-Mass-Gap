# Referee check for Part IV, Section 16: two auxiliary facts.
# (1) In Z^D (D = 3, 4, 5) a plaquette that is not a face of a unit cube Gamma contains at most one link of E(Gamma)
#     (the step "two links of E(Gamma) in a common plaquette span a face of Gamma" in the case t = 1 of Lemma mg:lem:seven).
# (2) The 60 closed ten-plaquette sets through a plaquette of Z^4 found by ref_closed.py are the 20 straight and 40 bent "two-cube boxes"
#     (two unit cubes sharing a face, the shared face removed); counted here directly from that description.
import itertools
from collections import Counter

def unit(D, i): return tuple(1 if k == i else 0 for k in range(D))
def add(x, y): return tuple(a + b for a, b in zip(x, y))
def links(D, p):
    x, i, j = p
    return {(x, i), (add(x, unit(D, i)), j), (add(x, unit(D, j)), i), (x, j)}
def cube(D, x, dirs):
    a, b, c = dirs; F = []
    for (i, j, k) in ((a, b, c), (a, c, b), (b, c, a)):
        i2, j2 = min(i, j), max(i, j)
        F.append((x, i2, j2)); F.append((add(x, unit(D, k)), i2, j2))
    return frozenset(F)

for D in (3, 4, 5):
    G = cube(D, tuple([0] * D), (0, 1, 2)); EG = set().union(*(links(D, p) for p in G))
    worst = 0
    for x in itertools.product(range(-2, 3), repeat=D):
        for i in range(D):
            for j in range(i + 1, D):
                p = (x, i, j)
                if p in G: continue
                worst = max(worst, len(links(D, p) & EG))
    print(f"(1) Z^{D}: max number of links of a unit cube in a plaquette that is not one of its faces: {worst}")

# (2) Z^4, p0 = (0, 01). Enumerate pairs of unit cubes sharing a face f, both containing f, union minus f, containing p0.
D = 4
p0 = (tuple([0] * D), 0, 1)
boxes = set()
rng = range(-2, 2)
cubes = []
for x in itertools.product(rng, repeat=D):
    for dirs in itertools.combinations(range(D), 3):
        cubes.append(cube(D, x, dirs))
for C1, C2 in itertools.combinations(cubes, 2):
    common = C1 & C2
    if len(common) != 1: continue
    Y = (C1 | C2) - common
    if p0 in Y: boxes.add(frozenset(Y))
def straight(Y):
    dirs = Counter((i, j) for (_, i, j) in Y)
    return len(dirs) == 3
cnt = Counter('straight' if straight(Y) else 'bent' for Y in boxes)
print(f"(2) Z^4: two-cube boxes through p0: {len(boxes)} ({dict(cnt)})")

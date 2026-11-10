# Exhaustive enumeration of small fronts through a vertex (Lemma mg:lem:smallfronts of Part IV).
#
# A front is a connected set of edges of the nearest-neighbour graph of the spatial torus in which every vertex has
# degree >= 2. N_k^v is the number of fronts with k edges containing a given vertex v (the torus is vertex-transitive).
# The script enumerates, by Redelmeier's algorithm (every connected edge set containing an edge at v exactly once),
# all connected edge sets with at most KMAX edges that contain v, and compares with
#   N_4^v = 2n(n-1) + #{sides = 4},   N_5^v = #{sides = 5},
#   N_6^v = 12 C(n,2) + 96 C(n,3) + 18(n-1) #{sides = 4} + #{sides = 6}     (n = d-1 spatial dimensions),
# and with the crude bound N_k^v <= n * C((4n-2)k, k-1) used for k >= 6; also P_v = 2(d-1)(d-2), the plaquettes at a vertex.
# Sides > k behave like Z^n for sets of k edges, so sides in {4,...,KMAX+1} cover every torus with sides >= 4.
import itertools, math, sys, time
from collections import defaultdict

def torus_graph(shape):
    n = len(shape)
    sites = list(itertools.product(*[range(L) for L in shape]))
    def add(x, i, s=1):
        y = list(x); y[i] = (y[i] + s) % shape[i]; return tuple(y)
    edges = [(x, i) for x in sites for i in range(n)]
    ends = {e: (e[0], add(e[0], e[1])) for e in edges}
    at = defaultdict(list)
    for e in edges:
        for v in ends[e]:
            at[v].append(e)
    nbr = {e: sorted(set(at[ends[e][0]]) | set(at[ends[e][1]]) - {e}) for e in edges}
    return sites, edges, ends, at, nbr

def count_fronts(shape, kmax):
    sites, edges, ends, at, nbr = torus_graph(shape)
    v0 = sites[0]
    counts = defaultdict(int)      # fronts (min degree >= 2) by size
    animals = defaultdict(int)     # all connected edge sets through v0 by size
    deg = defaultdict(int)
    current = []
    reached = set(at[v0])

    def low_degree_count():
        return sum(1 for v, c in deg.items() if c == 1)

    def rec(untried):
        untried = list(untried)
        while untried:
            e = untried.pop()
            current.append(e)
            for v in ends[e]:
                deg[v] += 1
            k = len(current)
            animals[k] += 1
            if all(c != 1 for c in deg.values()):
                counts[k] += 1
            if k < kmax:
                new = [f for f in nbr[e] if f not in reached]
                reached.update(new)
                rec(untried + new)
                reached.difference_update(new)
            for v in ends[e]:
                deg[v] -= 1
                if deg[v] == 0:
                    del deg[v]
            current.pop()

    rec(list(at[v0]))
    return counts, animals

def plaquettes_at_vertex(shape):
    # number of plaquettes of the spatial torus having a given vertex as a corner (P_v = 2(d-1)(d-2))
    n = len(shape); v0 = (0,) * n
    def add(x, i, s=1):
        y = list(x); y[i] = (y[i] + s) % shape[i]; return tuple(y)
    cnt = 0
    for x in itertools.product(*[range(L) for L in shape]):
        for i in range(n):
            for j in range(i + 1, n):
                if v0 in {x, add(x, i), add(x, j), add(add(x, i), j)}:
                    cnt += 1
    return cnt

def predicted(shape, k):
    n = len(shape)
    c = lambda L: sum(1 for x in shape if x == L)
    if k == 4:
        return 2 * n * (n - 1) + c(4)
    if k == 5:
        return c(5)
    if k == 6:
        return 12 * math.comb(n, 2) + 96 * math.comb(n, 3) + 18 * (n - 1) * c(4) + c(6)
    return None

def crude(n, k):
    return n * math.comb((4 * n - 2) * k, k - 1)

if __name__ == "__main__":
    KMAX = 6
    ok = True
    for n in (2, 3):
        sides = range(4, KMAX + 2)
        for shape in itertools.combinations_with_replacement(sides, n):
            t = time.time()
            counts, animals = count_fronts(shape, KMAX)
            line = []
            for k in range(4, KMAX + 1):
                p = predicted(shape, k)
                good = (p is None or counts[k] == p) and counts[k] <= crude(n, k)
                ok &= good
                line.append(f"N_{k}^v={counts[k]}" + ("" if p is None else f" (formula {p})") + ("" if good else "  <-- MISMATCH"))
            Pv = plaquettes_at_vertex(shape)
            ok &= Pv == 2 * n * (n - 1)
            line.append(f"P_v={Pv} (2(d-1)(d-2)={2*n*(n-1)})")
            print(f"d={n+1} spatial {shape}: " + ", ".join(line) + f";  animals through v: {[animals[k] for k in range(1, KMAX+1)]}  [{time.time()-t:.1f}s]")
            sys.stdout.flush()
    # k = 7 on the large torus (Z^n behaviour), for the record, against the crude bound
    for n, shape in ((2, (9, 9)), (3, (8, 8, 8))):
        t = time.time()
        counts, animals = count_fronts(shape, 7)
        print(f"d={n+1} spatial {shape} (as Z^{n}): N_7^v={counts[7]} vs crude bound {crude(n, 7)};  [{time.time()-t:.1f}s]")
        ok &= counts[7] <= crude(n, 7)
    print("all formulas confirmed" if ok else "MISMATCH FOUND")

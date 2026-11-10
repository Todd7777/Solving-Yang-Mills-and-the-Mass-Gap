# Lemma mg:lem:twofront of Part IV (Section 8). Exhaustive classification of small admissible junctions with exactly two front-pairs (entry p = (phi, s-1), exit out),
# on a spatial torus. phi is the boundary of the plaquette q0 (or a straight 4-cycle). Q = set of spatial plaquettes at slice s
# (the markers' plaquettes), |Q| <= NMAX. Admissibility at every edge e of E(Q) u phi u out:
#   n(e) = [e in F_{s-1}] + [e in F_s] + #{q in Q : e in q}  in {0} u [2,inf),
# with F_{s-1}, F_s the edges of the junction's fronts at the two layers (exactly p and out; no other temporal edges nearby).
# Also: fronts at the same layer are distinct components (vertex-disjoint), each front is touched by a plaquette of Q, and the
# junction (clusters of Q, p, out; adjacency = common vertex) is connected.
import itertools, sys
from collections import defaultdict, Counter
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from small_fronts import torus_graph

def setup(shape):
    n = len(shape)
    sites, edges, ends, at, nbr = torus_graph(shape)
    def add(x, i):
        y = list(x); y[i] = (y[i] + 1) % shape[i]; return tuple(y)
    plaqs = {}
    for x in sites:
        for i in range(n):
            for j in range(i + 1, n):
                plaqs[(x, i, j)] = frozenset([(x, i), (add(x, i), j), (add(x, j), i), (x, j)])
    return sites, edges, ends, at, nbr, plaqs

def fronts_of_size(shape, k, sites, edges, ends, at, nbr, near):
    # all fronts with k edges having a vertex in `near` (Redelmeier from each vertex in near)
    out = set()
    for v0 in near:
        reached = set(at[v0]); cur = []; deg = {}
        def rec(untried):
            untried = list(untried)
            while untried:
                e = untried.pop(); cur.append(e)
                for v in ends[e]: deg[v] = deg.get(v, 0) + 1
                if len(cur) == k:
                    if all(c >= 2 for c in deg.values()): out.add(frozenset(cur))
                else:
                    new = [f for f in nbr[e] if f not in reached]; reached.update(new); rec(untried + new); reached.difference_update(new)
                for v in ends[e]:
                    deg[v] -= 1
                    if deg[v] == 0: del deg[v]
                cur.pop()
        rec(list(at[v0]))
    return out

def run(shape, NMAX=2, KOUT=(4, 5, 6, 7, 8), straight=False):
    sites, edges, ends, at, nbr, plaqs = setup(shape)
    V = lambda E: {v for e in E for v in ends[e]}
    o = tuple(0 for _ in shape)
    if straight:
        phi = frozenset(((tuple((j if a == 0 else 0) for a in range(len(shape)))), 0) for j in range(shape[0]))
        assert len(phi) == 4
    else:
        phi = plaqs[(o, 0, 1)]
    Vphi = V(phi)
    # candidate plaquettes: within touching distance 2 of phi
    P1 = [p for p, E in plaqs.items() if V(E) & Vphi]
    V1 = Vphi | set().union(*[V(plaqs[p]) for p in P1])
    P2 = [p for p, E in plaqs.items() if V(E) & V1]
    results = Counter(); examples = {}
    near_out = V1 | set().union(*[V(plaqs[p]) for p in P2])
    outs = {k: fronts_of_size(shape, k, sites, edges, ends, at, nbr, near_out) for k in KOUT}
    for nQ in range(1, NMAX + 1):
        for Q in itertools.combinations(P2, nQ):
            EQ = [plaqs[p] for p in Q]
            EQall = set().union(*EQ)
            # markers' plaquettes must be touching-connected to phi (shown in the text); quick filter
            comp = set(Vphi); rest = list(EQ); changed = True
            while changed:
                changed = False
                for E in list(rest):
                    if V(E) & comp: comp |= V(E); rest.remove(E); changed = True
            if rest: continue
            cntQ = Counter(e for E in EQ for e in E)
            for k, OUTS in outs.items():
                for out in OUTS:
                    for lay in ('s', 's-1'):
                        if lay == 's-1' and (V(out) & Vphi or out == phi): continue      # distinct components at layer s-1
                        if lay == 's' and False: pass
                        Fm = set(phi) | (set(out) if lay == 's-1' else set())
                        Fp = set(out) if lay == 's' else set()
                        ok = True
                        for e in EQall | set(phi) | set(out):
                            n_e = (e in Fm) + (e in Fp) + cntQ[e]
                            if n_e == 1: ok = False; break
                        if not ok: continue
                        # fronts touched by Q; junction connectivity (clusters of Q, phi, out)
                        if not any(V(E) & Vphi for E in EQ): continue
                        if not any(V(E) & V(out) for E in EQ): continue
                        # clusters of Q
                        cl = []
                        for E in EQ:
                            merged = [c for c in cl if any(E & F for F in c)]
                            newc = [E] + [F for c in merged for F in c]
                            cl = [c for c in cl if c not in merged] + [newc]
                        nodes = [V(set().union(*c)) for c in cl] + [Vphi, V(out)]
                        seen = {0}; st = [0]
                        while st:
                            a = st.pop()
                            for b in range(len(nodes)):
                                if b not in seen and nodes[a] & nodes[b]: seen.add(b); st.append(b)
                        if len(seen) != len(nodes): continue
                        # the fronts must be exactly components: out is a front; at layer s-1 phi, out vertex-disjoint (checked)
                        key = (nQ, k, lay, 'phi=out' if out == phi else ('out=bdry' if any(out == plaqs[p] for p in Q) else 'other'))
                        results[key] += 1
                        examples.setdefault(key, (Q, sorted(out)))
    return results, examples, len(P1)

if __name__ == "__main__":
    # Lemma mg:lem:twofront: entry with four edges, exit with 4 (and 6) edges, at most two plaquettes (three for d=3).
    cases = [((6, 6), False, 3, (4,)), ((7, 7), False, 3, (4,)), ((4, 6), True, 3, (4,)),
             ((6, 6, 6), False, 2, (4, 6)), ((4, 6, 6), True, 2, (4, 6))]
    for shape, straight, nmax, kout in cases:
        res, ex, nP1 = run(shape, NMAX=nmax, KOUT=kout, straight=straight)
        print(f"spatial {shape}, entry {'straight 4-cycle' if straight else 'plaquette boundary'}, |Q| <= {nmax}, exits {kout}: "
              f"{nP1} plaquettes touch the entry")
        for key in sorted(res): print(f"   |Q|={key[0]} |exit|={key[1]} exit layer {key[2]:4s} {key[3]:9s}: {res[key]}")
        sys.stdout.flush()

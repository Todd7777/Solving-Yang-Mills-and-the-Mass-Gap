# Written by the independent referee of Section 13 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs. Code written from the definitions alone (no lemma of the paper is used).
# Independent check of Lemma 13.1(ii) and (iii) (one plaquette), from the definitions (see lattice.py, conditions J1-J4).
# p = (bd q, layer -1), q the plaquette at the origin in the plane (0,1); WLOG by lattice symmetry and time reflection.
#  (ii): one plaquette c, exactly two front-pairs p, p' with p' any front of <= KMAX edges (all connected min-degree-2 edge sets
#        through a vertex of c, found by brute force), either layer.  Lemma: |p'| >= 5 forces the other layer, c != q, and
#        p' = bd q xor bd c (6 edges) or bd q u bd c (7 edges), c sharing a link with q, or bd q u bd c (8 edges), c sharing exactly
#        one vertex with q; counts K_S, K_S, 4(P_v-1)-2K_S; none with 5 or >= 9 edges.
#  (iii): one plaquette c, exactly three front-pairs p, x, y with |x|, |y| <= 7: Lemma says impossible.
import sys, time, itertools
from collections import Counter, defaultdict
from lattice import *

def shift(F, v):
    return frozenset(E(add(a, v), add(b, v)) for (a, b) in F)

def run(n, KMAX_ii, KMAX_iii):
    t0 = time.time()
    o = tuple([0] * n)
    F0 = fronts_through(o, KMAX_ii, n)
    print(f"n={n} (d={n+1}): {len(F0)} fronts with <= {KMAX_ii} edges through a vertex  ({time.time()-t0:.1f}s)"); sys.stdout.flush()
    q = plaq(o, 0, 1, n); bq = q; Vq = verts(q)
    Pv = 2 * n * (n - 1); KS = 4 * (2 * n - 3)
    T = sorted(plaqs_touching(q, n), key=lambda p: sorted(p))
    def rel(c):
        k = len(c & q)
        if c == q: return 'c=q'
        if k == 1: return 'c shares a link'
        return 'c shares one vertex' if len(verts(c) & Vq) == 1 else 'c shares %d vertices?' % len(verts(c) & Vq)
    def fronts_touching(c, K):
        out = set()
        for v in verts(c):
            for f in F0:
                if len(f) <= K: out.add(shift(f, v))
        return [f for f in out if verts(f) & verts(c)]
    # ---------- (ii) ----------
    tally = Counter(); bad = []
    for c in T:
        Vc = verts(c)
        for f in fronts_touching(c, KMAX_ii):
            for lay in (-1, 1):
                if f == bq and lay == -1: continue
                # J1: same layer -> distinct and vertex-disjoint
                if lay == -1 and (verts(f) & Vq): continue
                # J4
                cnt = Counter()
                for e in bq: cnt[e] += 1
                for e in f: cnt[e] += 1
                for e in c: cnt[e] += 1
                if min(cnt.values()) < 2: continue
                # J2: both fronts touch c (the only marker); J3 then holds (c adjacent to both)
                if not (Vq & Vc): continue
                key = (len(f), 'same layer' if lay == -1 else 'other layer', rel(c),
                       'bd q ^ bd c' if f == (bq ^ c) else ('bd q | bd c' if f == (bq | c) else ('bd q' if f == bq else 'OTHER')))
                tally[key] += 1
    print(f"  (ii) one plaquette, two front-pairs, |p'| <= {KMAX_ii}:")
    for k, v in sorted(tally.items()): print(f"      {v:4d}  {k}")
    print(f"      expected: cap (4 edges, c=q, bd q, other layer) once; K_S={KS} with 6 and with 7 edges; 4(P_v-1)-2K_S={4*(Pv-1)-2*KS} with 8")
    sys.stdout.flush()
    # ---------- (iii) ----------
    t1 = time.time(); found = Counter(); npairs = 0
    for c in T:
        Vc = verts(c)
        cand = [f for f in fronts_touching(c, KMAX_iii)]
        # bit encoding of edges and vertices of this neighbourhood
        alle = set(bq) | set(c)
        for f in cand: alle |= f
        eidx = {e: i for i, e in enumerate(sorted(alle))}
        allv = verts(alle); vidx = {v: i for i, v in enumerate(sorted(allv))}
        em = lambda F: sum(1 << eidx[e] for e in F)
        vm = lambda F: sum(1 << vidx[v] for v in verts(F))
        Q, VQ, C = em(bq), vm(bq), em(c)
        FP = []
        for f in cand:
            for lay in (-1, 1):
                if f == bq and lay == -1: continue
                FP.append((em(f), vm(f), lay, len(f)))
        for (x, y) in itertools.combinations(FP, 2):
            npairs += 1
            ex, vx, lx, sx = x; ey, vy, ly, sy = y
            # J1 among the three fronts (p has layer -1)
            if lx == -1 and (vx & VQ): continue
            if ly == -1 and (vy & VQ): continue
            if lx == ly and ((vx & vy) or ex == ey): continue
            A = Q; B = 0
            if lx == -1: A |= ex
            else: B |= ex
            if ly == -1: A |= ey
            else: B |= ey
            U = A | B | C
            if U != ((A & B) | (A & C) | (B & C)): continue      # J4: every edge of E(J) in at least two of A, B, bd c
            # J2 holds by construction (every candidate touches c, and so must p): check p
            if not (Vq & Vc): continue
            found[(4, min(sx, sy), max(sx, sy))] += 1
    print(f"  (iii) one plaquette, three front-pairs (4, x, y), |x|,|y| <= {KMAX_iii}: {npairs} pairs tested, junctions found: {dict(found)} (Lemma: none)"
          f"  ({time.time()-t1:.0f}s)")
    sys.stdout.flush()

if __name__ == "__main__":
    n = int(sys.argv[1]); K2 = int(sys.argv[2]); K3 = int(sys.argv[3])
    run(n, K2, K3)

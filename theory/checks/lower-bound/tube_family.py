# The dominant family of connecting sets in Theorem mg:thm:fronts (Remark mg:rem:fronts(iii)), summed exactly.
#
# Family: one tube whose front in every layer is the boundary of a spatial plaquette, with a set Q of at most QMAX spatial plaquettes
# at each slice, subject to "no free link": every link at the slice lies in 0 or >= 2 plaquettes of X (the temporal plaquettes of the
# two adjacent layers count once each). With weights a per temporal plaquette and b per spatial plaquette (the weights W(X) of the proof),
# the sum over the family is a transfer-matrix sum: a^{4t} (M^{t-1})_{p0,.} with M[p,p'] = sum_Q b^{|Q|}.
# The per-layer growth a^4 lambda_max(M) is to be compared with a^4 (1 + delta) of Theorem mg:thm:fronts: lambda_max - 1 is about b,
# the cap plaquette (Q = {p}, p' = p), whereas delta ~ 128 P_v^2 beta. A compass for the constants, not part of any proof.
# (An independent referee's script, 10 November 2026, found the same growth factor a^4(1+b).)
import itertools, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fronts_constants import theorem_fronts, u_su2
from mpmath import mpf

def torus(shape):
    n = len(shape)
    sites = list(itertools.product(*[range(L) for L in shape]))
    def add(x, i):
        y = list(x); y[i] = (y[i] + 1) % shape[i]; return tuple(y)
    eidx = {(x, i): k for k, (x, i) in enumerate((x, i) for x in sites for i in range(n))}
    keys = [(x, i, j) for x in sites for i in range(n) for j in range(i + 1, n)]
    plaqs = [(eidx[(x, i)], eidx[(add(x, i), j)], eidx[(add(x, j), i)], eidx[(x, j)]) for (x, i, j) in keys]
    edge_plaqs = [[] for _ in eidx]
    for q, P in enumerate(plaqs):
        for e in P: edge_plaqs[e].append(q)
    return plaqs, edge_plaqs, len(eidx), keys

def entry(p, pp, plaqs, edge_plaqs, E, b, qmax):
    base = np.zeros(E, dtype=int)
    for e in plaqs[p]: base[e] += 1
    for e in plaqs[pp]: base[e] += 1
    local = sorted({q for e in set(plaqs[p]) | set(plaqs[pp]) for q in edge_plaqs[e]})
    tot = 0.0
    for r in range(qmax + 1):
        for Q in itertools.combinations(local, r):
            c = base.copy()
            for q in Q:
                for e in plaqs[q]: c[e] += 1
            if not np.any(c == 1):
                tot += b ** r
    return tot

def transfer(shape, b, qmax=3):
    """M[p,p'] = sum over admissible Q of b^|Q|; M is translation invariant, so one row per plaquette orientation is computed."""
    plaqs, edge_plaqs, E, keys = torus(shape)
    P = len(plaqs); index = {k: q for q, k in enumerate(keys)}
    origin = tuple(0 for _ in shape)
    rows = {}
    for (x, i, j) in keys:
        if x == origin:
            p0 = index[(x, i, j)]
            rows[(i, j)] = [entry(p0, pp, plaqs, edge_plaqs, E, b, qmax) for pp in range(P)]
    M = np.zeros((P, P))
    for p, (x, i, j) in enumerate(keys):
        for pp, (y, k, l) in enumerate(keys):
            rel = tuple((yc - xc) % L for yc, xc, L in zip(y, x, shape))
            M[p, pp] = rows[(i, j)][index[(rel, k, l)]]
    return M

if __name__ == "__main__":
    for d, shape in ((3, (4, 4)), (3, (5, 5)), (4, (4, 4, 4))):
        for beta in ('1e-7', '1e-6', '1e-5'):
            r = theorem_fronts(mpf(beta), d, u_su2(mpf(beta)), 'finite')
            b, delta = float(r['b']), float(r['delta'])
            M = transfer(shape, b)
            lam = max(abs(np.linalg.eigvals(M)))
            print(f"d={d} spatial {shape} beta={beta}: lambda_max(M)-1 = {lam-1:.4e}, b = {b:.4e}, (lambda_max-1)/b = {(lam-1)/b:.4f}; "
                  f"delta of Theorem mg:thm:fronts = {delta:.4e}")
            sys.stdout.flush()

# Fourth order (Part IV, Section 12): exact transfer-matrix spectra of Z_2 lattice gauge theory on spatial blocks of unit cubes
# (open boundary) and on tori, in dual plaquette variables, against the effective hopping matrix of the plaquette state.
#   States: the plaquette configurations b = (b_p) in Z_2^P reachable from 0 by flipping links (the image of the coboundary map;
#   for a block this is the Bianchi constraint on every cube, for a torus there are global constraints as well), i.e. the
#   gauge-invariant functions of the links that do not depend on the holonomies (the sector of the plaquette state).
#   T = V^{1/2} K V^{1/2}, V = prod_p exp(beta (-1)^{b_p}), K = prod_links (w0 + w1 X_l)/2, X_l = flip of the plaquettes at l.
#   One-plaquette levels: m_i = -log(lambda_i/lambda_0), i = 1..P; prediction m_i = 4 log(1/u) - log(1 - u^2) - log(1 + h_i u^4)
#   + O(u^6), u = tanh(beta), h_1 >= h_2 >= ... the eigenvalues of the effective matrix H on the plaquettes of the geometry:
#     H[q][q]  += #{c : c shares a link with q}        (bulge: the cube c x [t,t+1] attached to the side face of the tube)
#     H[c][q]  += 1 for every c sharing a link with q  (hop through the six-edge front of q u c)
#     H[q'][q] += 1 for the other five faces q' of every unit cube containing q   (four plaquettes of the cube at one slice: a
#                                                        turn to an adjacent face, or a jump to the opposite face; d >= 4)
import itertools, sys, math, time
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

def geometry(dims, periodic=False):
    """links (count), plaquettes (frozensets of link indices), cubes (lists of plaquette indices) of an open block of unit cubes
    [0,L_1] x ... x [0,L_n] or of the torus Z_{L_1} x ... x Z_{L_n}."""
    n = len(dims)
    def wrap(x): return tuple(xi % L for xi, L in zip(x, dims)) if periodic else tuple(x)
    def ok(x): return periodic or all(0 <= xi <= L for xi, L in zip(x, dims))
    rng = [range(L) if periodic else range(L + 1) for L in dims]
    verts = list(itertools.product(*rng))
    def e(i): return tuple(1 if j == i else 0 for j in range(n))
    def add(x, y): return tuple(a + b for a, b in zip(x, y))
    links = {}
    for x in verts:
        for i in range(n):
            if ok(add(x, e(i))): links[(x, i)] = len(links)
    def L(x, i): return links[(wrap(x), i)]
    plaqs, pidx = [], {}
    for x in verts:
        for i, j in itertools.combinations(range(n), 2):
            if not (ok(add(x, e(i))) and ok(add(x, e(j)))): continue
            E = frozenset([L(x, i), L(add(x, e(i)), j), L(add(x, e(j)), i), L(x, j)])
            pidx[(x, i, j)] = len(plaqs); plaqs.append(E)
    cubes = []
    for x in verts:
        for i, j, k in itertools.combinations(range(n), 3):
            if not all(ok(add(x, e(a))) for a in (i, j, k)): continue
            fs = [pidx[(x, i, j)], pidx[(wrap(add(x, e(k))), i, j)], pidx[(x, i, k)], pidx[(wrap(add(x, e(j))), i, k)],
                  pidx[(x, j, k)], pidx[(wrap(add(x, e(i))), j, k)]]
            cubes.append(fs)
    return len(links), plaqs, cubes

def reachable_basis(nl, plaqs):
    """GF(2) basis of the span of the link-flip vectors, and the coordinates of every link flip in that basis."""
    flips = [sum(1 << p for p, E in enumerate(plaqs) if l in E) for l in range(nl)]
    red = []                   # rows (pivot, vector); row k is basis vector k
    def reduce(v):
        c = 0
        for k, (pv, bv) in enumerate(red):
            if v >> pv & 1: v ^= bv; c ^= 1 << k
        return v, c
    for f in flips:
        v, _ = reduce(f)
        if v: red.append((v.bit_length() - 1, v))
    coords = []
    for f in flips:
        v, c = reduce(f)
        assert v == 0
        coords.append(c)
    return [bv for (pv, bv) in red], coords

def spectrum(dims, beta, periodic=False, k=2, tol=1e-13):
    nl, plaqs, cubes = geometry(dims, periodic)
    P = len(plaqs)
    basis, coords = reachable_basis(nl, plaqs)
    r = len(basis); D = 1 << r
    idx = np.arange(D, dtype=np.int64)
    # b_p(s) = xor_k s_k basis[k][p]
    nup = np.zeros(D, dtype=np.int32)
    for p in range(P):
        bp = np.zeros(D, dtype=np.int8)
        for kk, bv in enumerate(basis):
            if bv >> p & 1: bp ^= ((idx >> kk) & 1).astype(np.int8)
        nup += bp
    logV = beta*(P - 2*nup)
    sqrtV = np.exp(0.5*(logV - logV.max()))
    w0, w1 = math.exp(beta), math.exp(-beta)
    masks = coords
    def K(psi):
        for m in masks:
            psi = (w0*psi + w1*psi[idx ^ m])/2 if m else (w0 + w1)/2*psi
        return psi
    op = LinearOperator((D, D), matvec=lambda x: sqrtV*K(sqrtV*x), dtype=float)
    vals = eigsh(op, k=k, which='LA', tol=tol, return_eigenvectors=False, ncv=max(2*k + 1, 20))
    return np.sort(vals)[::-1], P, len(cubes), r

def effective_H(dims, periodic=False, cube_processes=True):
    nl, plaqs, cubes = geometry(dims, periodic)
    P = len(plaqs); H = np.zeros((P, P))
    for q in range(P):
        for c in range(P):
            if c != q and plaqs[q] & plaqs[c]:
                H[q, q] += 1; H[c, q] += 1
    if cube_processes:
        for fs in cubes:
            for q in fs:
                for qq in fs:
                    if qq != q: H[qq, q] += 1
    assert np.allclose(H, H.T)
    return np.sort(np.linalg.eigvalsh(H))[::-1]

def run(dims, periodic, betas=(0.2, 0.1, 0.05), nlev=1):
    h = effective_H(dims, periodic); h0 = effective_H(dims, periodic, cube_processes=False)
    P = len(h)
    print(f"{'torus' if periodic else 'block'} {'x'.join(map(str, dims))}: {P} plaquettes; top eigenvalues of H: "
          f"{np.round(h[:nlev], 6).tolist()} (without cube processes {np.round(h0[:nlev], 6).tolist()})")
    for beta in betas:
        t0 = time.time()
        vals, P, C, r = spectrum(dims, beta, periodic, k=nlev + 1)
        u = math.tanh(beta); base = 4*math.log(1/u) - math.log(1 - u*u)
        res = [(-math.log(vals[i]/vals[0]) - base)/u**4 for i in range(1, nlev + 1)]
        pred = [-math.log(1 + hh*u**4)/u**4 for hh in h[:nlev]]
        print(f"   beta={beta}: dim 2^{r}; (m_i - 4log(1/u) + log(1-u^2))/u^4 = {np.round(res, 4).tolist()}")
        print(f"   {'':{len(str(beta))+6}}  prediction -log(1 + h_i u^4)/u^4 = {np.round(pred, 4).tolist()}   [{time.time()-t0:.1f}s]")
        sys.stdout.flush()

if __name__ == "__main__":
    cases = {'3x3t': ((3, 3), True, 9), '4x4t': ((4, 4), True, 16), '1x1x1': ((1, 1, 1), False, 6),
             '1x1x2': ((1, 1, 2), False, 11), '1x2x2': ((1, 2, 2), False, 20), '1x2x3': ((1, 2, 3), False, 4),
             '2x2t': ((2, 2), True, 4), '3x4t': ((3, 4), True, 12)}
    for name in sys.argv[1:]:
        dims, per, nlev = cases[name]
        run(dims, per, nlev=nlev)

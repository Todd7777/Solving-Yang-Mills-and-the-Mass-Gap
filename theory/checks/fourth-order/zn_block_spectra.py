# Fourth order (Part IV, Section 12): exact transfer-matrix spectra of Z_n lattice gauge theory (n >= 3: complex characters) on
# spatial blocks of unit cubes (open boundary) and on tori, in dual plaquette variables b_p in Z_n, in the C-even and C-odd sectors,
# against the effective hopping matrices H_+ (unsigned) and H_- (signed by the orientations) of the plaquette state.
#   States: the configurations b reachable from 0 by link shifts (b_p -> b_p + [coefficient of l in the boundary of p]);
#   T = V^{1/2} K V^{1/2}, V = prod_p exp(beta cos(2 pi b_p/n)), K = prod_l sum_k exp(beta cos(2 pi k/n)) X_l^k (normalised).
#   C: b -> -b. One-plaquette levels in the sector C: m_i = 4 log(1/u) - log Upsilon_C - log(1 + h^C_i u^4) + O(beta^5),
#   u = c_1/c_0, Upsilon_C = Var(sqrt2 cos) (C = +) or Var(sqrt2 sin) (C = -) under exp(beta cos theta)/c_0,
#   h^C_i the eigenvalues of H_C:  H_C[q][q] += #{c sharing a link with q} (bulges, sign +);  H_C[c][q] += sigma(q,c) (hops);
#   H_C[q'][q] += sigma(q,q') for the other five faces q' of every unit cube containing q (turns and jumps); sigma = 1 for C = +,
#   and for C = - the orientation sign: the tube on q (standard orientation) ends on q' with orientation sigma(q,q').
import itertools, sys, math, time
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

def geometry(dims, periodic=False):
    """oriented plaquette boundaries {link: +-1} and cubes (lists of six plaquettes) of a block or torus."""
    nd = len(dims)
    def wrap(x): return tuple(xi % L for xi, L in zip(x, dims)) if periodic else tuple(x)
    def ok(x): return periodic or all(0 <= xi <= L for xi, L in zip(x, dims))
    rng = [range(L) if periodic else range(L + 1) for L in dims]
    verts = list(itertools.product(*rng))
    def e(i): return tuple(1 if j == i else 0 for j in range(nd))
    def add(x, y): return tuple(a + b for a, b in zip(x, y))
    links = {}
    for x in verts:
        for i in range(nd):
            if ok(add(x, e(i))): links[(x, i)] = len(links)
    def L(x, i): return links[(wrap(x), i)]
    plaqs, pidx = [], {}
    for x in verts:
        for i, j in itertools.combinations(range(nd), 2):
            if not (ok(add(x, e(i))) and ok(add(x, e(j)))): continue
            bd = {}
            for l, s in ((L(x, i), 1), (L(add(x, e(i)), j), 1), (L(add(x, e(j)), i), -1), (L(x, j), -1)):
                bd[l] = bd.get(l, 0) + s
            pidx[(x, i, j)] = len(plaqs); plaqs.append(bd)
    cubes = []
    for x in verts:
        for i, j, k in itertools.combinations(range(nd), 3):
            if not all(ok(add(x, e(a))) for a in (i, j, k)): continue
            cubes.append([pidx[(x, i, j)], pidx[(wrap(add(x, e(k))), i, j)], pidx[(x, i, k)], pidx[(wrap(add(x, e(j))), i, k)],
                          pidx[(x, j, k)], pidx[(wrap(add(x, e(i))), j, k)]])
    return len(links), plaqs, cubes

def chain_add(*terms):
    out = {}
    for s, bd in terms:
        for l, c in bd.items(): out[l] = out.get(l, 0) + s*c
    return {l: c for l, c in out.items() if c}

def sign_to(plaqs, q, S, q2):
    """sigma such that bd(q) + sum_{f in S} s_f bd(f) = sigma bd(q2) for some signs s_f (brute force); None if impossible."""
    for signs in itertools.product((1, -1), repeat=len(S)):
        front = chain_add((1, plaqs[q]), *[(s, plaqs[f]) for s, f in zip(signs, S)])
        for sg in (1, -1):
            if front == chain_add((sg, plaqs[q2])): return sg
    return None

def effective_H(dims, periodic=False, C='+'):
    nl, plaqs, cubes = geometry(dims, periodic)
    P = len(plaqs); H = np.zeros((P, P))
    for q in range(P):
        for c in range(P):
            if c != q and set(plaqs[q]) & set(plaqs[c]):
                H[q, q] += 1
                # hop: insert c at t (front bd(q u c)), then q at t+1: the front is +-bd(c)
                sg = sign_to(plaqs, q, [c, q], c) if False else None
                s_c = None
                for a in (1, -1):
                    fr = chain_add((1, plaqs[q]), (a, plaqs[c]))
                    if len(fr) == 6: s_c = a
                H[c, q] += 1 if C == '+' else s_c
    for fs in cubes:
        for q in fs:
            for q2 in fs:
                if q2 == q: continue
                S = [f for f in fs if f not in (q, q2)]
                sg = sign_to(plaqs, q, S, q2)
                assert sg is not None
                H[q2, q] += 1 if C == '+' else sg
    assert np.allclose(H, H.T)
    return np.sort(np.linalg.eigvalsh(H))[::-1]

def basis_mod_n(vectors, n):
    """basis (row echelon, mod prime n) of the span of integer vectors; returns rows and a function giving coordinates."""
    rows = []                                       # (pivot, vector)
    def reduce(v):
        v = [x % n for x in v]; c = {}
        for k, (pv, r) in enumerate(rows):
            if v[pv]:
                f = v[pv]                             # rows are normalised with r[pv] = 1
                v = [(a - f*b) % n for a, b in zip(v, r)]; c[k] = (c.get(k, 0) + f) % n
        return v, c
    for vec in vectors:
        v, _ = reduce(vec)
        nz = [i for i, x in enumerate(v) if x]
        if nz:
            pv = nz[0]; inv = pow(v[pv], n - 2, n)
            rows.append((pv, [(x*inv) % n for x in v]))
    def coords(vec):
        v, c = reduce(vec); assert not any(v); return [c.get(k, 0) for k in range(len(rows))]
    return [r for (pv, r) in rows], coords

def spectrum(n, dims, beta, periodic=False, k=2, sectors=('+', '-'), tol=1e-13):
    nl, plaqs, cubes = geometry(dims, periodic)
    P = len(plaqs)
    shifts = [[plaqs[p].get(l, 0) for p in range(P)] for l in range(nl)]
    rows, coords = basis_mod_n(shifts, n)
    r = len(rows); D = n**r
    digits = np.array(list(itertools.product(range(n), repeat=r))[::1], dtype=np.int16) if D <= 3**12 else None
    # index <-> digit vector, little-endian base n
    idx = np.arange(D, dtype=np.int64)
    dig = [(idx // n**kk) % n for kk in range(r)]
    b = np.zeros((P, D), dtype=np.int16)
    for p in range(P):
        acc = np.zeros(D, dtype=np.int64)
        for kk in range(r):
            if rows[kk][p]: acc += rows[kk][p]*dig[kk]
        b[p] = acc % n
    logV = beta*np.cos(2*np.pi*b/n).sum(axis=0)
    sqrtV = np.exp(0.5*(logV - logV.max()))
    w = np.array([math.exp(beta*math.cos(2*math.pi*kk/n)) for kk in range(n)]); w /= w.sum()
    shift_idx = []
    for l in range(nl):
        cl = coords(shifts[l])
        if not any(cl): shift_idx.append(None); continue
        perms = []
        for kk in range(1, n):
            newdig = [(dig[j] + kk*cl[j]) % n for j in range(r)]
            perms.append(sum(newdig[j]*n**j for j in range(r)))
        shift_idx.append(perms)
    neg = sum(((-dig[j]) % n)*n**j for j in range(r))          # C: b -> -b
    def K(psi):
        for perms in shift_idx:
            if perms is None: continue
            out = w[0]*psi
            for kk, pm in enumerate(perms): out = out + w[kk + 1]*psi[pm]
            psi = out
        return psi
    res = {}
    for C in sectors:
        sgn = 1 if C == '+' else -1
        proj = lambda x: 0.5*(x + sgn*x[neg])
        op = LinearOperator((D, D), matvec=lambda x: proj(sqrtV*K(sqrtV*proj(x))), dtype=float)
        vals = eigsh(op, k=k, which='LA', tol=tol, return_eigenvectors=False, ncv=max(2*k + 1, 20))
        res[C] = np.sort(vals)[::-1]
    return res, r

def zn_constants(n, beta):
    th = 2*np.pi*np.arange(n)/n; w = np.exp(beta*np.cos(th)); c0 = w.mean(); c1 = (w*np.cos(th)).mean()
    u = c1/c0
    ups = {}
    for C, f in (('+', np.sqrt(2)*np.cos(th)), ('-', np.sqrt(2)*np.sin(th))):
        m1 = (w*f).sum()/w.sum(); m2 = (w*f*f).sum()/w.sum(); ups[C] = m2 - m1**2
    return u, ups

def run(n, dims, periodic, betas=(0.04, 0.02), nlev=None):
    H = {C: effective_H(dims, periodic, C) for C in ('+', '-')}
    P = len(H['+']); nlev = nlev or P
    print(f"Z_{n} {'torus' if periodic else 'block'} {'x'.join(map(str, dims))}: {P} plaquettes")
    for C in ('+', '-'): print(f"   eigenvalues of H_{C}: {np.round(H[C], 6).tolist()}")
    for beta in betas:
        t0 = time.time()
        res, r = spectrum(n, dims, beta, periodic, k=nlev + 1)
        u, ups = zn_constants(n, beta)
        vac = max(res['+'][0], res['-'][0]); assert res['+'][0] >= res['-'][0]
        for C in ('+', '-'):
            lev = res[C][1:] if C == '+' else res[C][:nlev]
            base = 4*math.log(1/u) - math.log(ups[C])
            resid = [(-math.log(l/vac) - base)/u**4 for l in lev[:nlev]]
            pred = [-math.log(1 + hh*u**4)/u**4 for hh in H[C][:nlev]]
            print(f"   beta={beta} C{C}: dim {n}^{r}; (m_i - 4log(1/u) + log Ups_C)/u^4 = {np.round(resid, 3).tolist()}")
            print(f"   {'':{len(str(beta))+9}} prediction -log(1 + h_i u^4)/u^4 = {np.round(pred, 3).tolist()}   [{time.time()-t0:.1f}s]")
        sys.stdout.flush()

if __name__ == "__main__":
    cases = {'3x3t': ((3, 3), True), '1x1x1': ((1, 1, 1), False), '1x1x2': ((1, 1, 2), False)}
    n = int(sys.argv[1])
    for name in sys.argv[2:]:
        dims, per = cases[name]
        run(n, dims, per)

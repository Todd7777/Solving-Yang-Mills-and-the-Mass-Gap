# zn_links_tm.py -- written by the independent referee of Section 12 of Part IV (2026-11-12), added to the repository with credit;
# unchanged apart from this header. A compass, not part of the proofs.
# Referee code (independent of theory/checks): exact transfer matrix of Z_n lattice gauge theory (Wilson action
# w(U) = exp(beta Re chi(U)), chi(k) = exp(2 pi i k/n)) on a small spatial lattice, in LINK variables.
#   Configurations: a in Z_n^{links}.  Gauge orbits are computed by label propagation under the vertex generators
#   (a_l -> a_l + 1 for links leaving x, -1 for links entering x), optionally also under the centre transformations of a torus
#   (shift all direction-i links crossing the cut x_i = 0), which removes the winding (torelon) sectors.
#   T = V^{1/2} K V^{1/2} on gauge-invariant functions; in the basis of normalised orbit indicators
#       T_{OO'} = V^{1/2}(O) V^{1/2}(O') sum_{a' in O'} prod_l w(a0_l - a'_l)/n      (a0 any point of O; orbits of equal size).
#   Then: (i) the plaquette band in each C sector against -log(1 + h_i u^4); (ii) the equal-time and unit-time covariance
#   matrices of the plaquette functions psi_C(U_q) against Steps 2 and 5 of the proof of Theorem 12.8:
#       Cov(A_q, A_q')        = Ups delta + u^4 Ups^2 C_s + O(beta^5)
#       Cov(A_q(0), A_q'(1))  = u^4 Ups^2 (delta + u^4 (K_S delta + A_s + 2 C_s)) + O(beta^9)      [sum over ALL fronts]
import sys, math, time
import numpy as np
from geometry_H import Lattice, build_H, unit, vadd

def orbit_labels(lat, n, center=False):
    nl = len(lat.links); N = n ** nl
    idx = np.arange(N, dtype=np.int64)
    pw = n ** np.arange(nl, dtype=np.int64)
    dig = ((idx[None, :] // pw[:, None]) % n).astype(np.int8)          # (nl, N)
    gens = []
    for x in lat.verts:
        g = np.zeros(nl, dtype=np.int64)
        for k, (y, i) in enumerate(lat.links):
            if y == x: g[k] += 1
            if vadd(y, unit(lat.n, i), lat.Lw) == x: g[k] -= 1
        gens.append(g)
    if center:
        for i in range(lat.n):
            g = np.zeros(nl, dtype=np.int64)
            for k, (y, j) in enumerate(lat.links):
                if j == i and y[i] == 0: g[k] = 1
            gens.append(g)
    perms = []
    for g in gens:
        new = np.zeros(N, dtype=np.int64)
        for l in range(nl):
            new += ((dig[l].astype(np.int64) + g[l]) % n) * pw[l]
        perms.append(new)
    lab = idx.copy()
    while True:
        old = lab.copy()
        for p in perms:
            np.minimum(lab, lab[p], out=lab)
        if np.array_equal(lab, old): break
    reps, inv = np.unique(lab, return_inverse=True)
    return dig, reps, inv

def build(lat, n, center=False):
    t0 = time.time()
    dig, reps, inv = orbit_labels(lat, n, center)
    nl, N = dig.shape; norb = len(reps)
    sizes = np.bincount(inv, minlength=norb)
    assert len(set(sizes)) == 1, 'orbits of unequal size'
    # mismatch counts N_{O O'}(m), m = number of links with a'_l != a0_l  (K depends only on m for n = 2, 3)
    assert n in (2, 3)
    Ncnt = np.zeros((norb, norb, nl + 1), dtype=np.int64)
    for o, r in enumerate(reps):
        mism = (dig != dig[:, r][:, None]).sum(0)
        Ncnt[o] = np.bincount(inv * (nl + 1) + mism, minlength=norb * (nl + 1)).reshape(norb, nl + 1)
    # plaquette values b_p (standard orientation) at the representatives
    P = lat.P()
    B = np.zeros((P, norb), dtype=np.int64)
    for p, d in enumerate(lat.bd):
        for l, c in d.items():
            B[p] += c * dig[l, reps].astype(np.int64)
    B %= n
    # charge conjugation a -> -a on orbits
    neg_idx = np.zeros(norb, dtype=np.int64)
    pw = n ** np.arange(nl, dtype=np.int64)
    negrep = (((-dig[:, reps].astype(np.int64)) % n) * pw[:, None]).sum(0)
    conj = inv[negrep]
    print(f'   built: {N} configurations, {norb} orbits of size {sizes[0]}, {time.time()-t0:.1f}s'); sys.stdout.flush()
    return dict(Ncnt=Ncnt, B=B, conj=conj, nl=nl, norb=norb, n=n)

def T_matrix(data, beta):
    n, nl = data['n'], data['nl']
    th = 2 * np.pi * np.arange(n) / n
    w = np.exp(beta * np.cos(th))
    w0, w1 = w[0] / n, w[1] / n            # n = 2, 3: w(1) = w(n-1)
    m = np.arange(nl + 1)
    kern = (w0 ** (nl - m)) * (w1 ** m)
    S = np.cos(2 * np.pi * data['B'] / n).sum(0)
    sv = np.exp(0.5 * beta * (S - S.max()))
    T = (data['Ncnt'] * kern[None, None, :]).sum(2)
    T = sv[:, None] * T * sv[None, :]
    assert np.allclose(T, T.T, rtol=1e-12, atol=0)
    return T

def zn_consts(n, beta):
    th = 2 * np.pi * np.arange(n) / n
    w = np.exp(beta * np.cos(th)); c0 = w.mean(); pr = w / w.sum()
    u = (w * np.cos(th)).mean() / c0
    def var(f): return (pr * f * f).sum() - (pr * f).sum() ** 2
    if n == 2:
        return u, {'+': var(np.cos(th))}
    return u, {'+': var(np.sqrt(2) * np.cos(th)), '-': var(np.sqrt(2) * np.sin(th))}

def analyse(lat, data, beta, KS_block=True, verbose=True):
    n = data['n']; P = lat.P()
    T = T_matrix(data, beta)
    u, Ups = zn_consts(n, beta)
    Hp, Hm, parts, _ = build_H(lat)
    ev, evec = np.linalg.eigh(T)
    lam0 = ev[-1]; Om = evec[:, -1]
    if Om.sum() < 0: Om = -Om
    conj = data['conj']
    out = {}
    for C in (['+'] if n == 2 else ['+', '-']):
        # sector projector in orbit basis
        Pc = np.eye(data['norb']) + (1 if C == '+' else -1) * np.eye(data['norb'])[conj]
        Q, _ = np.linalg.qr(Pc)                     # columns span the sector (rank = dim of sector)
        rank = np.linalg.matrix_rank(Pc)
        U, s, _ = np.linalg.svd(Pc); Q = U[:, :rank]
        evC = np.sort(np.linalg.eigvalsh(Q.T @ T @ Q))[::-1]
        levels = evC[1:P + 1] if C == '+' else evC[:P]
        H = Hp if C == '+' else Hm
        h = np.sort(np.linalg.eigvalsh(H))[::-1]
        base = 4 * math.log(1 / u) - math.log(Ups[C])
        resid = (-np.log(levels / lam0) - base) / u ** 4
        pred = -np.log(1 + h * u ** 4) / u ** 4
        # covariance matrices
        f = np.sqrt(2) * (np.cos if C == '+' else np.sin)(2 * np.pi * data['B'] / n)
        if n == 2: f = np.cos(2 * np.pi * data['B'] / n)
        V = f * Om[None, :]                                  # rows: A_q Omega
        mean = V @ Om
        Vt = V - mean[:, None] * Om[None, :]
        G = Vt @ Vt.T
        T1 = Vt @ T @ Vt.T / lam0
        Cs = parts['Cp'] if C == '+' else parts['Cm']
        As = parts['Ap'] if C == '+' else parts['Am']
        KS = np.diag(parts['KS'])
        Y = Ups[C]
        off = ~np.eye(P, dtype=bool)
        gdev = np.abs(G[off] / (u ** 4 * Y ** 2) - Cs[off]).max()
        gdiag = np.abs(np.diag(G) - Y).max() / u ** 4
        t1off = np.abs(T1[off] / (u ** 8 * Y ** 2) - (As + 2 * Cs)[off]).max()
        t1diag = np.abs((np.diag(T1) / (u ** 4 * Y ** 2) - 1) / u ** 4 - parts['KS']).max()
        # effective band from the 2x2 data: eigenvalues of G^{-1/2} T1 G^{-1/2} / (u^4 Ups)
        gw, gv = np.linalg.eigh(G); Gm = gv @ np.diag(gw ** -0.5) @ gv.T
        eff = np.sort(np.linalg.eigvalsh(Gm @ T1 @ Gm))[::-1]
        heff = (eff / (u ** 4 * Y) - 1) / u ** 4
        out[C] = dict(resid=resid, pred=pred, h=h, gdev=gdev, gdiag=gdiag, t1off=t1off, t1diag=t1diag, heff=heff,
                      T1=T1, G=G, u=u, Y=Y)
        if verbose:
            print(f'  beta={beta} C{C}: u={u:.5f}, Ups={Y:.6f}')
            print(f'     levels (m_i - 4log(1/u) + log Ups)/u^4 = {np.round(resid, 3).tolist()}')
            print(f'     prediction -log(1+h_i u^4)/u^4         = {np.round(pred, 3).tolist()}   (h = {np.round(h, 6).tolist()})')
            print(f'     Step 2: max|G_qq\'/(u^4 Ups^2) - C| = {gdev:.2e} (off-diag), max|G_qq - Ups|/u^4 = {gdiag:.2e}')
            print(f'     Step 5: max|T1_qq\'/(u^8 Ups^2) - (A+2C)| = {t1off:.2e} (off-diag), max|(T1_qq/(u^4 Ups^2) - 1)/u^4 - K_S| = {t1diag:.2e}')
            print(f'     effective h from (G, T1): {np.round(heff, 3).tolist()}')
            sys.stdout.flush()
    return out

if __name__ == '__main__':
    which = sys.argv[1:] or ['z2cube', 'z3cube', 'z2t33']
    betas = (0.2, 0.1, 0.05)
    if 'z2cube' in which:
        lat = Lattice((1, 1, 1), False); print('Z_2, one cube (d=4), link variables')
        data = build(lat, 2)
        for b in betas: analyse(lat, data, b)
    if 'z3cube' in which:
        lat = Lattice((1, 1, 1), False); print('Z_3, one cube (d=4), link variables')
        data = build(lat, 3)
        for b in betas: analyse(lat, data, b)
    if 'z2t33' in which:
        lat = Lattice((3, 3), True); print('Z_2, 3x3 torus (d=3), link variables, gauge + centre orbits (winding zero)')
        data = build(lat, 2, center=True)
        for b in betas: analyse(lat, data, b)
    if 'z3t33' in which:
        lat = Lattice((3, 3), True); print('Z_3, 3x3 torus (d=3), link variables, gauge + centre orbits (winding zero)')
        data = build(lat, 3, center=True)
        for b in betas: analyse(lat, data, b)

# Written by the independent referee of Section 15 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs.
# Referee of Section 15 (2026-11-12). A compass, not part of the proofs.
# Lemma 15.4 (the hop and the turn) for G = Z_3, rho = faithful character, from the definitions:
#  - the two-copy space of the seven links of dq u dc (3^14 = 4.8e6 configurations), real arithmetic;
#  - G_g as the 9x9 kernel g(UV^-1,U'V'^-1) and G_g^{1/2} by matrix square root; Pi_0 x Pi_0 as the average;
#  - the gauge projection P^Gamma as the product of the averages over the gauge transformations at the six vertices, in each copy;
#  - the four faces of the turn: their g's integrated over the five far links of Gamma (on which every K is Pi_0 x Pi_0 on both sides);
#  - C acts on link variables (U -> U^-1 on every link), S exchanges the copies.
# No character theory is used for the operators: they are represented in the basis of normalised indicators of the gauge-invariant
# variables (U_q, U_c, U'_q, U'_c), which span the range of P^Gamma (checked).
# Checked: P commutes with the block; Y0_turn = Y_hop C; Q_-(Y_hop+Y_turn)Q_- has norm <= gamma_g u^12; the norm of Y_hop is u^12;
# the norms on the four joint sectors of (S, C); the twist read off the kernel (a function of U_q - sigma U_c), with sigma of
# Definition 12.1 and sigma_Gamma computed from d(dGamma) = 0.
import sys, time
import numpy as np
from geom import add, unit, boundary, sigma, cube_faces, eps_gamma

n = 3
beta = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
a_ = np.arange(n)
w = np.exp(beta * np.cos(2 * np.pi * a_ / n)); c0 = w.mean(); f = w / c0 - 1          # f_beta on Z_3
u = (f * np.exp(-2j * np.pi * a_ / n)).mean().real                                   # u = int f conj(chi_rho)
g2 = (1 + f[:, None]) * (1 + f[None, :]) - 1                                         # g(a, a')
gamma_g = np.abs(g2).max()
# single link, two copies: kernel (G h)(U,U') = (1/9) sum_{V,V'} g(U-V, U'-V') h(V,V')
G = np.zeros((n, n, n, n))
for U1 in range(n):
    for U2 in range(n):
        for V1 in range(n):
            for V2 in range(n):
                G[U1, U2, V1, V2] = g2[(U1 - V1) % n, (U2 - V2) % n] / n**2
Gm = G.reshape(n * n, n * n)
ev, evec = np.linalg.eigh(Gm)
assert ev.min() > -1e-12, ev.min()
Ghalf = (evec * np.sqrt(np.clip(ev, 0, None))) @ evec.T
Pi0 = np.full((n * n, n * n), 1.0 / n**2)
print(f'Z3, beta={beta}: u={u:.6f}, |G_g|={ev.max():.6f} (=u), gamma_g={gamma_g:.4f}')


def run(q, c):
    # the unit cube having both q and c as faces
    cands = []
    for dx in itertools_product():
        x0 = add(q[0], dx)
        F = cube_faces(x0)
        if q in F and c in F:
            cands.append(F)
    assert len(cands) == 1
    Gam = cands[0]
    eps = eps_gamma(Gam)
    sg = sigma(q, c); sgG = -eps[q] * eps[c]
    links = sorted(set(boundary(q)) | set(boundary(c)))
    far = sorted({l for p in Gam for l in boundary(p)} - set(links))
    assert len(links) == 7 and len(far) == 5
    verts = sorted({v for (b, d) in links for v in (b, add(b, unit(d)))})
    L = len(links)
    # configurations of one copy: 3^7; arrays of shape (3,)*7
    grids = np.indices((n,) * L)
    Uc = [grids[k] for k in range(L)]
    def hol(p, Ul, ls):
        h = 0
        for l, s in boundary(p).items():
            h = h + s * Ul[ls.index(l)]
        return h % n
    hq1, hc1 = hol(q, Uc, links), hol(c, Uc, links)
    # four faces integrated over the far links, one copy: F_S(U_7) = E_far prod_{p in S} f(U_p), S subset of the four faces
    rest = [p for p in Gam if p not in (q, c)]
    allL = links + far
    gr = np.indices((n,) * (L + 5))
    Ua = [gr[k] for k in range(L + 5)]
    fr = [f[hol(p, Ua, allL)] for p in rest]
    FS = {}
    for S in range(16):
        prod = np.ones((n,) * (L + 5))
        for i in range(4):
            if S >> i & 1:
                prod = prod * fr[i]
        FS[S] = prod.reshape((n,) * L + (n**5,)).mean(axis=-1)
    Mrest = np.zeros((n,) * (2 * L))
    for S1 in range(16):
        for S2 in range(16):
            if S1 | S2 == 15:
                Mrest += np.multiply.outer(FS[S1], FS[S2])
    gq = g2[hq1.reshape((n,) * L + (1,) * L), hq1.reshape((1,) * L + (n,) * L)]
    gc = g2[hc1.reshape((n,) * L + (1,) * L), hc1.reshape((1,) * L + (n,) * L)]
    idx_q = [links.index(l) for l in boundary(q)]
    idx_c = [links.index(l) for l in boundary(c)]
    hexl = [links.index(l) for l in links if (l in boundary(q)) != (l in boundary(c))]
    def K(psi, F, half):
        op_on = Ghalf if half else Gm
        for k in range(L):
            M = (op_on if k in F else Pi0).reshape(n, n, n, n)
            psi = np.moveaxis(np.tensordot(M, psi, axes=([2, 3], [k, L + k])), [0, 1], [k, L + k])
        return psi
    def Pgauge(psi):
        for copy in (0, L):
            for v in verts:
                acc = 0
                for h in range(n):
                    phi = psi
                    for k, (b, d) in enumerate(links):
                        sh = h * ((b == v) - (add(b, unit(d)) == v))
                        if sh % n:
                            phi = np.roll(phi, -sh, axis=copy + k)
                    acc = acc + phi
                psi = acc / n
        return psi
    # basis of range P: normalised indicators of (U_q, U_c, U'_q, U'_c)
    cls = (hq1.reshape((n,) * L + (1,) * L) * 27 + hc1.reshape((n,) * L + (1,) * L) * 9
           + hq1.reshape((1,) * L + (n,) * L) * 3 + hc1.reshape((1,) * L + (n,) * L))
    clsf = cls.ravel()
    def vec(j):
        return (n**2) * (cls == j).astype(float)
    # sanity: P fixes the basis vectors, and the basis spans range P (random test)
    rng = np.random.default_rng(0)
    r = rng.normal(size=(n,) * (2 * L)); Pr = Pgauge(r)
    co = np.bincount(clsf, weights=Pr.ravel(), minlength=81) * 81 / n**14
    recon = sum(co[j] * (cls == j) for j in range(81))
    err_span = np.abs(Pr - recon).max()
    err_fix = np.abs(Pgauge(vec(5)) - vec(5)).max()
    def blockop(kind, j, cap=True):
        psi = vec(j)
        if kind == 'hop':
            psi = K(psi, idx_c, True); psi = gq * psi; psi = K(psi, hexl, False); psi = gc * psi; psi = K(psi, idx_q, True)
        else:
            psi = K(psi, idx_c, True); psi = Mrest * psi; psi = K(psi, idx_q, False)
            if cap:
                psi = (1 + gq) * psi
            psi = K(psi, idx_q, True)
        return psi
    def matrix(kind, cap=True):
        Y = np.zeros((81, 81))
        for j in range(81):
            out = blockop(kind, j, cap)
            Y[:, j] = np.bincount(clsf, weights=out.ravel(), minlength=81) * 81 / n**14 / n**2
        return Y
    t0 = time.time()
    Yh = matrix('hop'); Yt = matrix('turn'); Y0 = matrix('turn', cap=False)
    # commutation with P (one column each)
    o = blockop('turn', 7); err_P = np.abs(Pgauge(o) - o).max() / np.abs(o).max()
    o = blockop('hop', 7); err_P = max(err_P, np.abs(Pgauge(o) - o).max() / np.abs(o).max())
    # C and S on the basis: (a,b,a',b') -> (-a,-b,-a',-b') and -> (a',b',a,b)
    Cm = np.zeros((81, 81)); Sm = np.zeros((81, 81))
    for j in range(81):
        A, B, A2, B2 = j // 27, j // 9 % 3, j // 3 % 3, j % 3
        Cm[((-A) % 3) * 27 + ((-B) % 3) * 9 + ((-A2) % 3) * 3 + (-B2) % 3, j] = 1
        Sm[A2 * 27 + B2 * 9 + A * 3 + B, j] = 1
    I = np.eye(81)
    Qm, Qp, Pm, Pp = (I - Cm) / 2, (I + Cm) / 2, (I - Sm) / 2, (I + Sm) / 2
    nrm = lambda X: np.linalg.norm(X, 2)
    print(f'  q={q}, c={c}: sigma(q,c)={sg:+d}, sigma_Gamma(q,c)={sgG:+d}  (Lemma 12.2(i): sigma_Gamma=-sigma: {sgG == -sg})')
    print(f'    P^Gamma: basis fixed {err_fix:.1e}, spans range {err_span:.1e}; block commutes with P: {err_P:.1e}   [{time.time()-t0:.0f}s]')
    print(f'    |Y_hop| = {nrm(Yh):.6e}, u^12 = {u**12:.6e}, ratio {nrm(Yh)/u**12:.10f}')
    print(f'    |Y0_turn - Y_hop C| / |Y_hop| = {nrm(Y0 - Yh @ Cm)/nrm(Yh):.1e}')
    print(f'    |Q-(Y_hop+Y0_turn)Q-| / |Y_hop| = {nrm(Qm @ (Yh + Y0) @ Qm)/nrm(Yh):.1e}')
    r1 = nrm(Qm @ (Yh + Yt) @ Qm)
    print(f'    |Q-(Y_hop+Y_turn)Q-| = {r1:.4e} <= |Y_turn-Y0_turn| = {nrm(Yt-Y0):.4e} <= gamma_g u^12 = {gamma_g*u**12:.4e}: {r1 <= nrm(Yt-Y0)*(1+1e-9) <= gamma_g*u**12}')
    for (sn, PS) in (('S-', Pm), ('S+', Pp)):
        for (cn, PC) in (('C-', Qm), ('C+', Qp)):
            X = PS @ PC @ (Yh + Yt) @ PC @ PS
            print(f'    sector ({sn},{cn}): |P(Y_hop+Y_turn)P| / u^12 = {nrm(X)/u**12:.6f}')
    # twist: the kernel of Y_hop, as a function of (U_q, U_c) in copy 1 (copy 2 fixed at the trivial class), depends on U_q - sigma U_c
    Kh = Yh.reshape(3, 3, 3, 3, 3, 3, 3, 3)[:, 0, 0, 0, 0, :, 0, 0]   # out (a_q, a_c=0, a'_q=0, a'_c=0) <- in (0, b_c, 0, 0)
    Kt = Y0.reshape(3, 3, 3, 3, 3, 3, 3, 3)[:, 0, 0, 0, 0, :, 0, 0]
    def dep(Kmat, s):
        return max(abs(Kmat[A, B] - Kmat[(A - s * B) % 3, 0]) for A in range(3) for B in range(3))
    print(f'    hop kernel a function of U_q - sigma U_c: {dep(Kh, sg):.1e} (of U_q + sigma U_c: {dep(Kh, -sg):.1e});'
          f'  turn kernel a function of U_q - sigma_Gamma U_c: {dep(Kt, sgG):.1e}')
    return nrm(Y0 - Yh @ Cm) / nrm(Yh), r1 / (gamma_g * u**12)


def itertools_product():
    import itertools
    return [d for d in itertools.product((-1, 0), repeat=3)]


if __name__ == '__main__':
    q = ((0, 0, 0), 0, 1)
    worst = 0
    for c in (((0, 0, 0), 0, 2), ((0, 0, -1), 0, 2), ((0, 0, 0), 1, 2), ((1, 0, 0), 1, 2)):
        e1, e2 = run(q, c)
        worst = max(worst, e1)
    print(f'worst |Y0_turn - Y_hop C|/|Y_hop| = {worst:.1e}')

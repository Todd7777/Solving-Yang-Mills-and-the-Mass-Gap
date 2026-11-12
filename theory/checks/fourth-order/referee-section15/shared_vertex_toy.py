# Written by the independent referee of Section 15 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs.
# Referee of Section 15 (2026-11-12). A compass, not part of the proofs.
# The step in the proof of Lemma 15.7(ii): at a vertex v shared by the local operator L and the rest R, the gauge projection factorises,
#     P_v (L x R) P_v = (P^L_v L P^L_v) x (P^R_v R P^R_v),
# if one side of L is Pi_0 x Pi_0 on all of L's links at v (L and R commuting with the gauge transformations). For abelian G the test is
# degenerate (a single-vertex gauge transformation acts trivially once the other vertices are projected), so we use G = S3 (two copies).
# Toy graph: two digons sharing the vertex 0 (L on a: 0->1, b: 1->0; R on c: 0->2, d: 2->0); operators built from G_g^{1/2}, Pi_0 and
# multiplication by g of the digon holonomies, with the explicit gauge projection (average over all of S3 at each vertex, in each copy).
# Negative control: G_g^{1/2} instead of Pi_0 on the left side of L.
import sys
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 2)[0] + '/referee-section14')
from groups import S3
Gr = S3(); n = Gr.n; irho = int(np.argmax(Gr.dims))          # the two-dimensional irrep
f, u = Gr.setup(irho, 0.9)
g2 = (1 + f[:, None]) * (1 + f[None, :]) - 1
UVi = Gr.mul[:, Gr.inv]                                          # UVi[U, V] = U V^-1
G = np.array([[g2[UVi[U1, V1], UVi[U2, V2]] / n**2 for V1 in range(n) for V2 in range(n)] for U1 in range(n) for U2 in range(n)])
ev, Vv = np.linalg.eigh(G); Gh = (Vv * np.sqrt(np.clip(ev, 0, None))) @ Vv.T
Pi0 = np.full((n * n, n * n), 1 / n**2)
links = [(0, 1), (1, 0), (0, 2), (2, 0)]                          # a, b (L) ; c, d (R)
NL = len(links)
idx = np.indices((n,) * (2 * NL))
hL = [Gr.mul[idx[c + 0], idx[c + 1]] for c in (0, NL)]
hR = [Gr.mul[idx[c + 2], idx[c + 3]] for c in (0, NL)]
gL, gR = g2[hL[0], hL[1]], g2[hR[0], hR[1]]

def op1(M, psi, k):
    return np.moveaxis(np.tensordot(M.reshape(n, n, n, n), psi, axes=([2, 3], [k, NL + k])), [0, 1], [k, NL + k])

def gauge(psi, vert, which):
    """Average over g in S3 of psi(U^g), U^g_e = g_{s(e)} U_e g_{t(e)}^-1, at one vertex, in each copy, on the links `which`."""
    for copy in (0, NL):
        acc = 0
        for g in range(n):
            phi = psi
            for k, (s, t) in enumerate(links):
                if k in which and (s == vert or t == vert):
                    perm = np.arange(n)
                    if s == vert:
                        perm = Gr.mul[g, perm]
                    if t == vert:
                        perm = Gr.mul[perm, Gr.inv[g]]
                    phi = np.take(phi, perm, axis=copy + k)
            acc = acc + phi
        psi = acc / n
    return psi

def P(psi, spec):
    for vert, which in spec:
        psi = gauge(psi, vert, which)
    return psi

Pfull = [(0, [0, 1, 2, 3]), (1, [0, 1]), (2, [2, 3])]
PL, PR = [(0, [0, 1]), (1, [0, 1])], [(0, [2, 3]), (2, [2, 3])]

def Lop(psi, pi_left):
    for k in (0, 1):
        psi = op1(Gh, psi, k)
    psi = gL * psi
    for k in (0, 1):
        psi = op1(Pi0 if pi_left else Gh, psi, k)
    return psi

def Rop(psi):
    for k in (2, 3):
        psi = op1(Gh, psi, k)
    psi = gR * psi
    for k in (2, 3):
        psi = op1(Gh, psi, k)
    return psi

rng = np.random.default_rng(1)
for pi_left in (True, False):
    errs, sizes = [], []
    for _ in range(3):
        psi = rng.normal(size=(n,) * (2 * NL))
        lhs = P(Lop(Rop(P(psi, Pfull)), pi_left), Pfull)
        rhs = P(Lop(P(P(Rop(P(psi, PR)), PR), PL), pi_left), PL)
        errs.append(np.abs(lhs - rhs).max() / np.abs(lhs).max()); sizes.append(np.abs(lhs).max())
    print(f'S3, left side of L is Pi_0 at the shared vertex: {pi_left}:  |P(LxR)P psi| ~ {min(sizes):.2e};'
          f'  max |P(LxR)P - (PLP)x(PRP)| / |P(LxR)P| = {max(errs):.1e}')

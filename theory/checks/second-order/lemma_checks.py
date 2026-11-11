# Brute-force checks of three lemmas of Section 10 of Part IV ("The second-order term"). A compass, not part of the proofs.
#  (A) eq. (mg:eq:cubereal): for a real class function a that need not be invariant under inversion,
#      int a(U_p0) a(U_p1) prod_{4 sides} f(U_s) dU = sum_r |a^_r|^2 u_r^4   (one unit cube, maximal-tree gauge, 5 free links)
#      groups Z_3, Z_5 (faithful character) and A_4 (three-dimensional representation; a separates the two classes of 3-cycles)
#  (B) Lemma mg:lem:dressed: on L^2(G^4 x G^4) (four links of a square, two copies), D_q = (G_g^{1/2})^{(x)4} M_{1+g_q} (G_g^{1/2})^{(x)4};
#      the largest eigenvalue of P_- P^q D_q P^q P_- (swap-odd sector, gauge projection at the corners), matrix-free Lanczos on the full
#      space, against u^4 Upsilon and the bound u^4 Upsilon-bar; and ||P^q D_q P^q|| <= u^4 (1+gamma_g). Z_3, Z_4, S_3.
#  (C) Lemma mg:lem:twoplaq: for a shift (two squares sharing a link, 7 links) and a jump (sharing a vertex, 8 links), two copies,
#      P^J O_J P^J has rank one, singular value 2 sigma_6 + sigma_6^2, and vanishes on the swap-odd sector. Z_2, Z_3 (shift), Z_2 (jump).
# Usage: python3 lemma_checks.py [A|B|C]   (default: all; B with S_3 and C with Z_3 take a few minutes)
import numpy as np, itertools, math, sys
from scipy.sparse.linalg import LinearOperator, eigsh, svds

def group_Zn(n):
    els = list(range(n)); om = np.exp(2j*np.pi/n)
    irreps = {str(k): (lambda k: (lambda g: np.array([[om**(k*g)]])))(k) for k in range(n)}
    return dict(els=els, mul=lambda a, b: (a+b) % n, inv=lambda a: (-a) % n, irreps=irreps, rho='1', name=f'Z{n}')

def group_S3():
    els = list(itertools.permutations(range(3)))
    mul = lambda a, b: tuple(a[b[i]] for i in range(3)); inv = lambda a: tuple(int(i) for i in np.argsort(a))
    B = np.array([[1, -1, 0], [1, 1, -2]], float); B[0] /= np.sqrt(2); B[1] /= np.sqrt(6)
    def perm(p):
        P = np.zeros((3, 3))
        for i in range(3): P[p[i], i] = 1
        return P
    sgn = lambda p: round(np.linalg.det(perm(p)))
    irreps = {'0': lambda g: np.array([[1.]]), 'sgn': lambda g: np.array([[float(sgn(g))]]), '2': lambda g: B@perm(g)@B.T}
    return dict(els=els, mul=mul, inv=inv, irreps=irreps, rho='2', name='S3')

def group_A4():
    def parity(p): return sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j]) % 2
    els = [p for p in itertools.permutations(range(4)) if parity(p) == 0]
    mul = lambda a, b: tuple(a[b[i]] for i in range(4)); inv = lambda a: tuple(int(i) for i in np.argsort(a))
    V4 = [p for p in els if all(p[p[i]] == i for i in range(4))]      # identity and the double transpositions
    c = (1, 2, 0, 3); c2 = mul(c, c); e = (0, 1, 2, 3)
    def k(g): return next(j for j, cj in enumerate((e, c, c2)) if mul(g, inv(cj)) in V4)
    om = np.exp(2j*np.pi/3)
    B = np.array([[1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]], float); B /= np.linalg.norm(B, axis=1)[:, None]
    def perm(p):
        P = np.zeros((4, 4))
        for i in range(4): P[p[i], i] = 1
        return P
    irreps = {'0': lambda g: np.array([[1.]]), 'w': lambda g: np.array([[om**k(g)]]), 'w2': lambda g: np.array([[om**(2*k(g))]]),
              '3': lambda g: B@perm(g)@B.T}
    return dict(els=els, mul=mul, inv=inv, irreps=irreps, rho='3', name='A4')

def setup(G, beta):
    els = G['els']; n = len(els); idx = {g: i for i, g in enumerate(els)}
    chi = {r: np.array([np.trace(R(g)) for g in els]) for r, R in G['irreps'].items()}
    d = {r: G['irreps'][r](els[0]).shape[0] for r in chi}
    rho = G['rho']; N = d[rho]; iota = 1 if np.allclose(chi[rho], np.conj(chi[rho])) else 0
    w = np.exp(beta/N*np.real(chi[rho])); c0 = w.mean(); f = w/c0-1
    u = {r: np.real(np.mean(np.conj(chi[r])*w))/(d[r]*c0) for r in chi}
    mul, inv = G['mul'], G['inv']
    MUL = np.array([[idx[mul(a, b)] for b in els] for a in els]); INV = np.array([idx[inv(a)] for a in els])
    return dict(n=n, idx=idx, chi=chi, d=d, N=N, iota=iota, f=f, u=u, MUL=MUL, INV=INV, ur=u[rho], rho=rho)

# ---------------------------------------------------------------- (A) the cube identity for real class functions
def check_A():
    rng = np.random.default_rng(1)
    for G in (group_Zn(3), group_Zn(5), group_A4()):
        for beta in (0.3, 0.05):
            S = setup(G, beta); n, MUL, INV, f, chi = S['n'], S['MUL'], S['INV'], S['f'], S['chi']
            # a real class function, not invariant under inversion: random real combination of Re/Im of the characters
            a = np.zeros(n)
            for r in chi: a += rng.normal()*np.real(chi[r]) + rng.normal()*np.imag(chi[r])
            assert not np.allclose(a, a[INV])
            ahat = {r: np.mean(a*np.conj(chi[r])) for r in chi}
            rhs = sum(abs(ahat[r])**2*S['u'][r]**4 for r in chi if r != '0')
            other = sum((ahat[r]**2).real*S['u'][r]**4 for r in chi if r != '0')        # what a wrong orientation would give
            # free links x1 = U_DC, x2 = U_A'B', x3 = U_B'C', x4 = U_D'C', x5 = U_A'D' (tree links = identity)
            x = np.indices((n,)*5).reshape(5, -1)
            Up0 = INV[x[0]]
            Up1 = MUL[MUL[MUL[x[1], x[2]], INV[x[3]]], INV[x[4]]]
            lhs = np.mean(a[Up0]*a[Up1]*f[x[1]]*f[x[2]]*f[MUL[x[3], INV[x[0]]]]*f[x[4]])
            print(f"(A) {G['name']:3s} beta={beta}: lhs={lhs:+.12e}  sum_(r!=0)|a_r|^2u_r^4={rhs:+.12e}  (sum Re(a_r^2)u_r^4={other:+.4e})  rel.diff={abs(lhs-rhs)/abs(rhs):.1e}")

# ---------------------------------------------------------------- shared: link operators and gauge projections
def link_ops(S):
    n, MUL, INV, f = S['n'], S['MUL'], S['INV'], S['f']
    F = f[MUL[np.arange(n)[:, None], INV[np.arange(n)[None, :]]]]/n       # kernel f(U V^{-1}) / |G|
    Pi0 = np.ones((n, n))/n
    Gg = np.kron(Pi0+F, Pi0+F)-np.kron(Pi0, Pi0)
    ev, V = np.linalg.eigh(Gg); Gh = (V*np.sqrt(np.clip(ev, 0, None)))@V.T
    return Gh, np.kron(Pi0, Pi0)

class Graph:
    """links (start, end) on vertices; configurations of both copies as an array of shape (n,)*L x (n,)*L."""
    def __init__(self, S, links):
        self.S, self.links, self.L = S, links, len(links); n = S['n']
        self.verts = sorted({v for l in links for v in l})
        self.grid = np.indices((n,)*self.L).reshape(self.L, -1)
    def hol(self, path):          # path: list of (link index, +1/-1); returns holonomy index per configuration
        MUL, INV = self.S['MUL'], self.S['INV']; h = None
        for l, s in path:
            x = self.grid[l] if s > 0 else INV[self.grid[l]]
            h = x if h is None else MUL[h, x]
        return h
    def gauge_perms(self, v):     # permutations of the one-copy configurations for the gauge transformations at vertex v
        MUL, INV, n = self.S['MUL'], self.S['INV'], self.S['n']; perms = []
        for g in range(n):
            y = self.grid.copy()
            for l, (a, b) in enumerate(self.links):
                if a == v: y[l] = MUL[g, y[l]]
                if b == v: y[l] = MUL[y[l], INV[g]]
            perms.append(np.ravel_multi_index(tuple(y), (n,)*self.L))
        return perms
    def projector(self):
        P = [self.gauge_perms(v) for v in self.verts]; n = self.S['n']
        def apply(M):            # M shape (D, D): rows copy 1, columns copy 2
            for perms in P:
                M = sum(M[p, :] for p in perms)/n
                M = sum(M[:, p] for p in perms)/n
            return M
        return apply
    def apply_links(self, M, ops):   # ops[l]: operator on the pair (U_l, U'_l) (n^2 x n^2) or None
        n, L = self.S['n'], self.L
        v = M.reshape((n,)*(2*L))
        for l, O in enumerate(ops):
            if O is None: continue
            v = np.moveaxis(v, (l, L+l), (0, 1)); sh = v.shape
            v = (O@v.reshape(n*n, -1)).reshape(sh)
            v = np.moveaxis(v, (0, 1), (l, L+l))
        return v.reshape(M.shape)

def upsilon(S):
    """Upsilon and the bound Upsilon-bar of Lemma mg:lem:dressed, with the exact u_2 and sigma_6 of the finite group."""
    chi, rho, f, u, d, N, iota = S['chi'], S['rho'], S['f'], S['u'], S['d'], S['N'], S['iota']
    w1 = 1+f
    basis = [np.sqrt(2)*np.real(chi[rho]), np.sqrt(2)*np.imag(chi[rho])] if iota == 0 else [np.real(chi[rho])]
    V = np.array([[np.mean(p*q*w1)-np.mean(p*w1)*np.mean(q*w1) for q in basis] for p in basis])
    Ups = np.linalg.eigvalsh(V).max()
    others = [r for r in chi if r != '0' and not np.allclose(chi[r], chi[rho]) and not np.allclose(chi[r], np.conj(chi[rho]))]
    u2 = max([u[r] for r in others], default=0.0); ur = S['ur']; beta_eps = np.max(np.abs(f))
    s6 = sum(d[r]**2*u[r]**6 for r in chi if r != '0'); gam = np.max(np.abs(np.add.outer(f, f)+np.multiply.outer(f, f)))
    v2 = sum(d[r]**2*u[r]**6 for r in others)
    nB = (u2/ur)**2*beta_eps + math.sqrt(2-iota)*N/ur*math.sqrt(v2)
    ups1 = Ups + nB**2/(Ups-(u2/ur)**4*np.max(w1)) if nB > 0 else Ups
    Ubar = ups1 + 2*np.max(w1)**2*s6/(ups1-ur**4*(1+gam))
    return Ups, Ubar, gam, s6

# ---------------------------------------------------------------- (B) the dressed tube on the swap-odd sector
def check_B(groups):
    for G in groups:
        for beta in (0.3, 0.1):
            S = setup(G, beta); n = S['n']
            g = Graph(S, [(0, 1), (1, 2), (3, 2), (0, 3)])          # square a=0,b=1,c=2,d=3: U_q = U1 U2 U3^-1 U4^-1
            fq = S['f'][g.hol([(0, 1), (1, 1), (2, -1), (3, -1)])]
            one_plus_g = (1+fq)[:, None]*(1+fq)[None, :]
            Gh, _ = link_ops(S); P = g.projector(); D1 = n**4
            def Dq(x, odd):
                M = x.reshape(D1, D1)
                if odd: M = (M-M.T)/2
                M = P(M); M = g.apply_links(M, [Gh]*4); M = one_plus_g*M; M = g.apply_links(M, [Gh]*4); M = P(M)
                if odd: M = (M-M.T)/2
                return M.reshape(-1)
            ur = S['ur']; Ups, Ubar, gam, s6 = upsilon(S)
            top_odd = eigsh(LinearOperator((D1*D1,)*2, matvec=lambda x: Dq(x, True), dtype=float), k=1, which='LA', tol=1e-12,
                            return_eigenvectors=False)[0]
            top_all = eigsh(LinearOperator((D1*D1,)*2, matvec=lambda x: Dq(x, False), dtype=float), k=1, which='LA', tol=1e-12,
                            return_eigenvectors=False)[0]
            print(f"(B) {G['name']} beta={beta}: u={ur:.5e}  odd sector: lambda_max/u^4={top_odd/ur**4:.10f}  Upsilon={Ups:.10f}  "
                  f"Upsilon-bar={Ubar:.10f}  |  ||P D P||/u^4={top_all/ur**4:.6f} <= 1+gamma_g={1+gam:.6f}")
            sys.stdout.flush()

# ---------------------------------------------------------------- (C) shifts and jumps: rank one, swap-even
def check_C(cases):
    for G, kind in cases:
        for beta in (0.3, 0.1):
            S = setup(G, beta); n = S['n']
            if kind == 'shift':   # squares [0,1]x[0,1] and [1,2]x[0,1]; vertices (x,y) -> 3y+x
                links = [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]       # h00 h10 h01 h11 v00 v10 v20
                q = [(0, 1), (5, 1), (2, -1), (4, -1)]; qp = [(1, 1), (6, 1), (3, -1), (5, -1)]
            else:                 # jump: squares [0,1]x[0,1] and [1,2]x[1,2] (common vertex (1,1)); vertices (x,y) -> 3y+x
                links = [(0, 1), (3, 4), (0, 3), (1, 4), (4, 5), (7, 8), (4, 7), (5, 8)]
                q = [(0, 1), (3, 1), (1, -1), (2, -1)]; qp = [(4, 1), (7, 1), (5, -1), (6, -1)]
            g = Graph(S, links); L = g.L; D1 = n**L
            Eq = {l for l, _ in q}; Eqp = {l for l, _ in qp}
            gq = S['f'][g.hol(q)]; gqp = S['f'][g.hol(qp)]
            Gtwo = lambda h: (1+h)[:, None]*(1+h)[None, :]-1
            M = Gtwo(gq)*Gtwo(gqp)
            Gh, Pi = link_ops(S); P = g.projector()
            Kin = [Gh if l in Eq else Pi for l in range(L)]           # K~^{1/2}_{J,s-1}: G_g^{1/2} on dq, Pi_0 x Pi_0 elsewhere
            Kout = [Gh if l in Eqp else Pi for l in range(L)]         # K~^{1/2}_{J,s}
            def mv(x): X = P(x.reshape(D1, D1)); X = g.apply_links(X, Kout); X = M*X; X = g.apply_links(X, Kin); return P(X).reshape(-1)
            def rmv(x): X = P(x.reshape(D1, D1)); X = g.apply_links(X, Kin); X = M*X; X = g.apply_links(X, Kout); return P(X).reshape(-1)
            op = LinearOperator((D1*D1,)*2, matvec=mv, rmatvec=rmv, dtype=float)
            sv = np.sort(svds(op, k=3, tol=1e-12, return_singular_vectors=False))[::-1]
            def odd(x): X = x.reshape(D1, D1); return ((X-X.T)/2).reshape(-1)
            opo = LinearOperator((D1*D1,)*2, matvec=lambda x: odd(mv(odd(x))), rmatvec=lambda x: odd(rmv(odd(x))), dtype=float)
            svo = svds(opo, k=1, tol=1e-12, return_singular_vectors=False)[0]
            s6 = sum(S['d'][r]**2*S['u'][r]**6 for r in S['chi'] if r != '0')
            print(f"(C) {G['name']} {kind} beta={beta}: singular values {sv[0]:.10e} {sv[1]:.2e} {sv[2]:.2e};  2 sigma6+sigma6^2={2*s6+s6**2:.10e};  "
                  f"odd sector: {svo:.2e}")
            sys.stdout.flush()

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else 'ABC'
    if 'A' in which: check_A()
    if 'B' in which: check_B([group_Zn(3), group_Zn(4), group_S3()])
    if 'C' in which: check_C([(group_Zn(2), 'shift'), (group_Zn(2), 'jump'), (group_Zn(3), 'shift')])

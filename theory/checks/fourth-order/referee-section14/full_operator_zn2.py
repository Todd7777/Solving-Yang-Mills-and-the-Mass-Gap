# Written by the independent referee of Section 14 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs.
# ||P^J O_J P^J|| from the DEFINITIONS of Lemma 9.2 (pieces), for Z_n gauge theory, without the class-function reduction of Lemma 14.1:
#   O_J = K_{J,s-1}^{1/2} M_J K_{J,s}^{1/2} on L^2(G^{E(J)} x G^{E(J)}),  K^{1/2} = tensor product over links of G_g^{1/2} (front links) or Pi0 x Pi0,
#   G_g = integral operator with kernel g(UV^-1, U'V'^-1) (matrix square root computed numerically), M_J = prod_p g_p, P^J = average over
#   the gauge transformations at all vertices, in each copy. Top singular value by power iteration; compared with u_*^k.
import numpy as np, itertools, sys

def setup(n, beta):
    z = np.arange(n); chi = np.exp(2j * np.pi * z / n)
    w = np.exp(beta * chi.real); c0 = w.mean(); f = w / c0 - 1
    u = np.array([(f * np.exp(-2j * np.pi * r * z / n)).mean().real for r in range(n)])
    # one-link two-copy operator G_g:  (G psi)(x,x') = 1/n^2 sum_{y,y'} g(x-y, x'-y') psi(y,y')
    g = lambda a, b: (1 + f[a % n]) * (1 + f[b % n]) - 1
    Gm = np.zeros((n * n, n * n))
    for x, xp, y, yp in itertools.product(range(n), repeat=4):
        Gm[x * n + xp, y * n + yp] = g(x - y, xp - yp) / n**2
    ev, V = np.linalg.eigh(Gm)
    assert ev.min() > -1e-13, ev
    Gh = (V * np.sqrt(np.clip(ev, 0, None))) @ V.T
    Pi = np.ones((n * n, n * n)) / n**2
    return f, u, Gh, Pi

class Junction:
    def __init__(self, n, links, faces, front_minus, front_plus, P):
        """links: list of (v,w) oriented v->w; faces: dict name -> list of (link index, sign); fronts: sets of link indices; P: face names."""
        self.n, self.links, self.L = n, links, len(links)
        self.verts = sorted({v for l in links for v in l})
        self.faces, self.fm, self.fp, self.P = faces, front_minus, front_plus, P

    def holo(self, X, face):     # X: array (n,)*L of link values for one copy -> holonomy of the face, as an array broadcast over configs
        return sum(s * X[i] for i, s in self.faces[face]) % self.n

    def build(self, beta):
        n, L = self.n, self.L
        self.f, self.u, Gh, Pi = setup(n, beta)
        grids = np.meshgrid(*[np.arange(n)] * (2 * L), indexing='ij')    # axis 2e = copy 1 of link e, 2e+1 = copy 2
        X = [grids[2 * e] for e in range(L)]; Xp = [grids[2 * e + 1] for e in range(L)]
        M = np.ones([n] * (2 * L))
        for p in self.P:
            M = M * ((1 + self.f[self.holo(X, p)]) * (1 + self.f[self.holo(Xp, p)]) - 1)
        self.M = M
        self.Km = [Gh if e in self.fm else Pi for e in range(L)]
        self.Kp = [Gh if e in self.fp else Pi for e in range(L)]
        del grids, X, Xp

    def applyK(self, K, v):
        n, L = self.n, self.L
        t = v.reshape([n * n] * L)
        for e in range(L):
            t = np.moveaxis(np.tensordot(K[e], t, axes=([1], [e])), 0, e)
        return t.reshape([n] * (2 * L))

    def gauge(self, v):
        n = self.n
        for vert in self.verts:
            for copy in (0, 1):
                acc = np.zeros_like(v)
                for h in range(n):
                    t = v
                    for e, (a, b) in enumerate(self.links):
                        ax = 2 * e + copy
                        if a == vert: t = np.roll(t, h, axis=ax)
                        if b == vert: t = np.roll(t, -h, axis=ax)
                    acc += t
                v = acc / n
        return v

    def A(self, v):      # P O P
        v = self.gauge(v); v = self.applyK(self.Kp, v); v = self.M * v; v = self.applyK(self.Km, v); return self.gauge(v)

    def AT(self, v):
        v = self.gauge(v); v = self.applyK(self.Km, v); v = self.M * v; v = self.applyK(self.Kp, v); return self.gauge(v)

    def O(self, v):
        v = self.applyK(self.Kp, v); v = self.M * v; return self.applyK(self.Km, v)

    def OT(self, v):
        v = self.applyK(self.Km, v); v = self.M * v; return self.applyK(self.Kp, v)

    def commute_test(self, seed=1):
        rng = np.random.default_rng(seed)
        v = rng.normal(size=[self.n] * (2 * self.L))
        a = self.gauge(self.O(v)); b = self.O(self.gauge(v))
        return np.linalg.norm(a - b) / np.linalg.norm(a)

    def top_singular(self, iters=14, seed=0):
        # P commutes with O (checked numerically by commute_test), so (POP)^T(POP) = P O^T O P; keep v in range(P)
        rng = np.random.default_rng(seed)
        v = self.gauge(rng.normal(size=[self.n] * (2 * self.L))); v /= np.linalg.norm(v)
        s = 0
        for it in range(iters):
            w = self.O(v); s_new = np.linalg.norm(w)          # = ||P O P v|| since v = P v and P commutes with O
            v = self.gauge(self.OT(w)); v /= np.linalg.norm(v)
            if it > 3 and abs(s_new - s) <= 1e-12 * s_new: s = s_new; break
            s = s_new
        return s

def domino_junction(n, bent):
    # vertices as tuples; links oriented in the positive coordinate direction
    if not bent:
        qv = [(0, 0), (1, 0), (1, 1), (0, 1)]; cv = [(1, 0), (2, 0), (2, 1), (1, 1)]
    else:
        qv = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]; cv = [(0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)]
    links = sorted({(min(c[t], c[(t + 1) % 4]), max(c[t], c[(t + 1) % 4])) for c in (qv, cv) for t in range(4)})
    def face(cyc):
        out = []
        for t in range(4):
            a, b = cyc[t], cyc[(t + 1) % 4]
            if (a, b) in links: out.append((links.index((a, b)), 1))
            else: out.append((links.index((b, a)), -1))
        return out
    faces = {'q': face(qv), 'c': face(cv)}
    lq = {i for i, s in faces['q']}; lc = {i for i, s in faces['c']}
    return links, faces, lq, lq ^ lc

def cube_junction(n, kind):
    V = list(itertools.product((0, 1), repeat=3))
    links = sorted((a, b) for a in V for b in V if a < b and sum(abs(a[k] - b[k]) for k in range(3)) == 1)
    faces = {}
    for k in range(3):
        for s in (0, 1):
            i, j = [m for m in range(3) if m != k]
            def P_(x, y): return tuple(s if m == k else (x if m == i else y) for m in range(3))
            cyc = [P_(0, 0), P_(1, 0), P_(1, 1), P_(0, 1)]
            out = []
            for t in range(4):
                a, b = cyc[t], cyc[(t + 1) % 4]
                out.append((links.index((a, b)), 1) if (a, b) in links else (links.index((b, a)), -1))
            faces[(k, s)] = out
    q = (2, 0); qp = (0, 0) if kind == 'turn' else (2, 1)
    P = [F for F in faces if F not in (q, qp)]
    return links, faces, {i for i, s in faces[q]}, {i for i, s in faces[qp]}, P

if __name__ == '__main__':
    for n in (2, 3):
        for beta in ((0.3, 0.9) if n == 2 else (0.9,)):
            for bent in (False, True):
                links, faces, fm, fp = domino_junction(n, bent)
                for name, (a, b) in (('widening', (fm, fp)), ('narrowing', (fp, fm))):
                    J = Junction(n, links, faces, a, b, ['c']); J.build(beta)
                    ct = J.commute_test(); s = J.top_singular(); us = max(J.u[1:])
                    print(f'Z{n} beta={beta} {"bent" if bent else "flat"} domino {name:9s}: ||P O P|| = {s:.12e}, u_*^6 = {us**6:.12e}, ratio-1 = {s/us**6-1:.1e}, [P,O] rel. {ct:.1e}', flush=True)
    n = 2
    for beta in (0.9,):
        for kind in ('turn', 'jump'):
            links, faces, fm, fp, P = cube_junction(n, kind)
            J = Junction(n, links, faces, fm, fp, P); J.build(beta)
            s = J.top_singular(iters=10); us = max(J.u[1:])
            print(f'Z{n} beta={beta} cube {kind}: ||P O P|| = {s:.12e}, u_*^8 = {us**8:.12e}, ratio-1 = {s/us**8-1:.1e}', flush=True)

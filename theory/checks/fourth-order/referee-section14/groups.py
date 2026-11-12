# Written by the independent referee of Section 14 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs.
# Finite groups with numerically computed character tables (Dixon-style: simultaneous eigenvectors of class-multiplication matrices).
# Independent of the repository's scripts. Elements are indices 0..n-1 with a multiplication table.
import numpy as np, itertools

def closure(gens, key):
    els = [np.eye(gens[0].shape[0], dtype=complex)]
    seen = {key(els[0]): 0}
    i = 0
    while i < len(els):
        for g in gens:
            h = els[i] @ g
            k = key(h)
            if k not in seen:
                seen[k] = len(els); els.append(h)
        i += 1
    n = len(els)
    mul = np.zeros((n, n), dtype=np.int64)
    for a in range(n):
        for b in range(n):
            mul[a, b] = seen[key(els[a] @ els[b])]
    return els, mul

def mkkey(M):
    return tuple(np.round(M.real, 6).ravel()) + tuple(np.round(M.imag, 6).ravel())

class Group:
    def __init__(self, name, gens):
        self.name = name
        self.mats, self.mul = closure(gens, mkkey)
        n = self.n = len(self.mats)
        self.e = 0
        self.inv = np.array([[b for b in range(n) if self.mul[a, b] == 0][0] for a in range(n)])
        # conjugacy classes
        cls = -np.ones(n, dtype=int); classes = []
        for a in range(n):
            if cls[a] < 0:
                orb = sorted({self.mul[self.mul[g, a], self.inv[g]] for g in range(n)})
                for x in orb: cls[x] = len(classes)
                classes.append(orb)
        self.cls, self.classes = cls, classes
        h = len(classes)
        # class multiplication: C_i C_j = sum_k a_ijk C_k ; a_ijk = #{(x,y) in K_i x K_j: xy = z_k}
        A = np.zeros((h, h, h))
        for i in range(h):
            for j in range(h):
                for x in classes[i]:
                    for y in classes[j]:
                        A[i, j, cls[self.mul[x, y]]] += 1
                for k in range(h):
                    A[i, j, k] /= len(classes[k])
        rng = np.random.default_rng(1)
        c = rng.normal(size=h)
        Mt = sum(c[i] * A[i].T for i in range(h))     # (M_i)_{kj} = a_ij^k ; left eigenvectors -> eigenvectors of M^T
        # M_i acting: (M_i)_{kj}=A[i,j,k] -> M_i = A[i].T ; left eigvec w: w M_i = omega_i w  <=> M_i^T w^T = omega w^T, M_i^T = A[i]
        Mt = sum(c[i] * A[i] for i in range(h))
        w, V = np.linalg.eig(Mt)
        chars = []
        for t in range(h):
            om = V[:, t] / V[0, t]            # omega(C_1)=1 (class 0 is the identity)
            # chi(g_k) = om_k * chi(1) / |K_k| ; normalise sum_k |K_k| |chi_k|^2 = n
            vals = np.array([om[k] / len(classes[k]) for k in range(h)])
            nrm = sum(len(classes[k]) * abs(vals[k])**2 for k in range(h))
            d = np.sqrt(n / nrm)
            chik = vals * d
            chars.append(chik[cls])          # as a function on elements
        # sort: trivial first, then by dimension
        chars.sort(key=lambda ch: (round(ch[0].real), -round(np.mean(ch).real, 6), round(np.angle(ch[1]) if abs(ch[1]) > 1e-9 else 0, 4)))
        self.chars = np.array(chars)
        self.dims = np.array([int(round(ch[0].real)) for ch in chars])
        # checks: orthonormality and sum d^2 = n
        G = (self.chars @ self.chars.conj().T) / n
        assert np.allclose(G, np.eye(h), atol=1e-9), (name, G)
        assert sum(self.dims**2) == n
        self.R = h
        self.conjperm = np.array([[s for s in range(h) if np.allclose(self.chars[s], self.chars[r].conj())][0] for r in range(h)])

    def setup(self, irho, beta):
        N = self.dims[irho]
        w = np.exp(beta / N * self.chars[irho].real)
        c0 = w.mean(); f = w / c0 - 1
        u = np.array([(f * self.chars[r].conj()).mean().real / self.dims[r] for r in range(self.R)])
        ui = np.array([(f * self.chars[r].conj()).mean().imag for r in range(self.R)])
        assert np.allclose(ui, 0, atol=1e-12)
        return f, u

def Q8():
    i = np.array([[1j, 0], [0, -1j]]); j = np.array([[0, 1], [-1, 0]], dtype=complex)
    return Group('Q8', [i, j])

def A4():
    def P(p):
        M = np.zeros((4, 4), dtype=complex)
        for a, b in enumerate(p): M[b, a] = 1
        return M
    return Group('A4', [P((1, 2, 0, 3)), P((1, 0, 3, 2))])

def F21():
    z = np.exp(2j * np.pi / 7)
    a = np.diag([z, z**2, z**4]); b = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=complex)
    return Group('F21', [a, b])

def S3():
    def P(p):
        M = np.zeros((3, 3), dtype=complex)
        for a, b in enumerate(p): M[b, a] = 1
        return M
    return Group('S3', [P((1, 2, 0)), P((1, 0, 2))])

def Zn(n):
    return Group(f'Z{n}', [np.array([[np.exp(2j * np.pi / n)]])])

if __name__ == '__main__':
    for G in (Q8(), A4(), F21(), S3(), Zn(4)):
        print(G.name, G.n, 'dims', G.dims, 'conj', G.conjperm)

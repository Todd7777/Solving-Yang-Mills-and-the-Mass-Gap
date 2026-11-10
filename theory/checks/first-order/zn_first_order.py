# Corollary mg:cor:firstorder2 (Section 9 of Part IV), exact finite-volume check with abelian groups, d = 3:
#  * Z_3 on the 3x3 spatial torus, zero-winding sector (3^8 states): kappa(beta) = u(beta) exactly; the lowest C-even state has
#    m - 4log(1/u) = -u(1+O(u)) and the lowest C-odd state +u(1+O(u)) (the cap mixes the plaquette states chi and chi-bar);
#  * Z_2 on the 4x4 spatial torus, zero-winding sector (2^15 states): kappa = 0, and m - 4log(1/u) = +u^2(1+O(u)).
# The 3x3 torus has girth 3 (the theorem needs spatial sides >= 4 for the lower bound); in the zero-winding sector the torelons
# are absent and the plaquette glueball is the lightest state. A compass, not part of the proofs.
import numpy as np, math
from scipy.sparse.linalg import LinearOperator, eigsh

def z3(beta):
    L = 3; nl = 2*L*L; shape = (3,)*8; om = np.exp(2j*np.pi/3)
    def lid(x, y, i): return ((x % L)*L + (y % L))*2 + i
    def bd(x, y):
        v = np.zeros(nl, dtype=np.int64)
        v[lid(x, y, 0)] += 1; v[lid(x+1, y, 1)] += 1; v[lid(x, y+1, 0)] -= 1; v[lid(x, y, 1)] -= 1
        return v
    P = [bd(x, y) for x in range(L) for y in range(L)]
    labels = (np.indices(shape).reshape(8, -1).T @ np.array(P[:8])) % 3      # electric charges of the 3^8 states
    def coeffs(fun): return np.array([np.real(sum(np.conj(om**(k*m))*fun(om**m) for m in range(3))/3) for k in range(3)])
    c = coeffs(lambda U: math.exp(beta*U.real)); b = coeffs(lambda U: math.exp(beta*U.real/2))
    Kd = np.prod(c[labels]/c[0], axis=1)
    def shift(v, p, k):
        if p < 8: return np.roll(v, k, axis=p)
        return np.roll(v, (-k,)*8, axis=tuple(range(8)))     # the ninth plaquette boundary is minus the sum of the others
    def Vop(v):
        v = v.reshape(shape)
        for p in range(9): v = (b[0]*v + b[1]*shift(v, p, 1) + b[2]*shift(v, p, -1))/b[0]
        return v.reshape(-1)
    idxC = np.ix_(*[np.mod(-np.arange(3), 3)]*8)
    def Cop(v): return np.ascontiguousarray(v.reshape(shape)[idxC]).reshape(-1)
    n = 3**8; top = {}
    for sgn in (+1, -1):
        Pc = lambda v: 0.5*(v + sgn*Cop(v))
        op = LinearOperator((n, n), matvec=lambda v: Pc(Vop(Kd*Vop(Pc(v)))), dtype=float)
        top[sgn] = np.sort(eigsh(op, k=3, which='LA', tol=1e-14, return_eigenvectors=False))[::-1]
    u = c[1]/c[0]; kappa = 2*np.mean(np.cos(2*np.pi*np.arange(3)/3)**2*(np.exp(beta*np.cos(2*np.pi*np.arange(3)/3))/c[0]-1))
    return u, kappa, math.log(top[1][0]/top[1][1]), math.log(top[1][0]/top[-1][0])

def z2(beta):
    L = 4; nl = 2*L*L
    def lid(x, y, i): return ((x % L)*L + (y % L))*2 + i
    def bd(x, y):
        v = np.zeros(nl, dtype=np.int64)
        for l in (lid(x, y, 0), lid(x+1, y, 1), lid(x, y+1, 0), lid(x, y, 1)): v[l] ^= 1
        return v
    P = [bd(x, y) for x in range(L) for y in range(L)]
    n = 2**15; bits = ((np.arange(n)[:, None] >> np.arange(15)[None, :]) & 1)
    popc = ((bits @ np.array(P[:15])) % 2).sum(1); idx = np.arange(n)
    flip = [idx ^ (1 << k) for k in range(15)] + [idx ^ ((1 << 15)-1)]
    th, tb = math.tanh(beta), math.tanh(beta/2); Kd = th**popc
    def Vop(v):
        for fl in flip: v = v + tb*v[fl]
        return v
    op = LinearOperator((n, n), matvec=lambda v: Vop(Kd*Vop(v)), dtype=float)
    w = np.sort(eigsh(op, k=3, which='LA', tol=1e-14, return_eigenvectors=False))[::-1]
    return th, math.log(w[0]/w[1])

if __name__ == "__main__":
    print("Z_3, 3x3 torus, zero winding:  kappa = u;  (m - 4log(1/u))/u for the lowest C-even and C-odd states")
    for beta in (0.03, 0.01, 0.003, 0.001):
        u, kappa, me, mo = z3(beta)
        print(f"  beta={beta:<6g} u={u:.5e} kappa/u={kappa/u:.12f}   C-even {(me-4*math.log(1/u))/u:+.5f}   C-odd {(mo-4*math.log(1/u))/u:+.5f}")
    print("Z_2, 4x4 torus, zero winding:  kappa = 0;  (m - 4log(1/u))/u and /u^2")
    for beta in (0.03, 0.01, 0.003, 0.001):
        u, m = z2(beta)
        print(f"  beta={beta:<6g} u={u:.5e}   (m-4log(1/u))/u = {(m-4*math.log(1/u))/u:+.3e}   /u^2 = {(m-4*math.log(1/u))/u**2:+.4f}")

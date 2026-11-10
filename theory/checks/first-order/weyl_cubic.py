# Lemma mg:lem:capelement(iii) (Section 9 of Part IV): int (Re chi)^3 dU over SU(N) for the fundamental representation, by Weyl
# integration on the maximal torus (trapezoidal rule, exact for trigonometric polynomials of low degree): 0 for N = 2, 4, 5 and
# 1/4 for N = 3 (int chi^3 = 1, the determinant); also int chi^2 chi-bar = 0. And kappa(beta) for SU(3) against beta/6 = 3u.
import numpy as np, math
def torus(N, M=48):
    th = np.linspace(0, 2*np.pi, M, endpoint=False)
    grids = np.meshgrid(*([th]*(N-1)), indexing='ij')
    ang = [g.ravel() for g in grids]; ang.append(-sum(ang))
    z = np.exp(1j*np.array(ang)); vd = np.ones(z.shape[1])
    for i in range(N):
        for j in range(i+1, N): vd *= np.abs(z[i]-z[j])**2
    return z, vd/(math.factorial(N)*M**(N-1))
if __name__ == "__main__":
    for N in (2, 3, 4, 5):
        z, w = torus(N); X = z.sum(0)
        print(f"SU({N}): int 1 = {np.sum(w):.12f}; int (Re chi)^3 = {np.sum(np.real(X)**3*w):+.12f}; "
              f"int chi^3 = {np.sum(X**3*w).real:+.10f}; |int chi^2 chibar| = {abs(np.sum(X**2*np.conj(X)*w)):.1e}")
    z, w = torus(3, 64); X = z.sum(0)
    for beta in (1e-2, 1e-3, 1e-4, 1.538e-4, 4.026e-4):
        W = np.exp(beta/3*np.real(X)); c0 = np.sum(W*w); u = np.sum(np.conj(X)*W*w).real/(3*c0)
        kappa = np.sum(np.real(X)**2*(W/c0-1)*w)/0.5
        print(f"  SU(3) beta={beta:.4e}: u={u:.8e}  kappa={kappa:.8e}  kappa/(beta/6)={kappa/(beta/6):.6f}  kappa/(3u)={kappa/(3*u):.6f}")

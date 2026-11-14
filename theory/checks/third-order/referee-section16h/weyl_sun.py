# Referee (Section 16h) -- independent group data of SU(N) by Weyl integration on the maximal torus (numpy, double precision).
# Usage: python3 weyl_sun.py [N] [M]     (default N=5, M=28: a 28^4 grid; about 1 minute)
"""Weyl integration: for a class function F on SU(N),
      int F dU = (1/N!) * mean over the torus grid of F(diag z) |Delta(z)|^2,
z_j = e^{i theta_j} (j < N), z_N = e^{-i(theta_1+...+theta_{N-1})}, theta_j on a uniform M-point grid (trapezoidal rule, exact for
trigonometric polynomials whose frequencies in each theta_j are < M in absolute value).  Characters are Schur polynomials, by Jacobi-Trudi
from the complete homogeneous polynomials h_k (Newton's identities from the power sums), so nothing is divided by Delta.
The weight is w_beta = exp((beta/N) Re chi_rho) (the normalisation of Section mg:sec:upper: y = (beta/N) Re chi in the proof of Lemma mg:lem:char(iv),
t = beta/(2N) in its (iii)).  Computed: c_0, u = c_rho/(N c_0), U_+ = Var(sqrt2 Re chi), U_- = Var(sqrt2 Im chi) under w/c_0, and the character
coefficients u_r = c_r/(d_r c_0) of the irreps of low order, against Lemma mg:lem:u2sun.  Nothing here is imported from the authors' code."""
import sys, itertools
import numpy as np
from math import comb, factorial, expm1

def torus(N, M):
    th = 2*np.pi*np.arange(M)/M
    grids = np.meshgrid(*([th]*(N - 1)), indexing='ij')
    ang = [g.ravel() for g in grids]
    ang.append(-sum(ang))
    z = np.exp(1j*np.array(ang))                       # shape (N, M^(N-1))
    D = np.ones(z.shape[1], dtype=complex)
    for i in range(N):
        for j in range(i + 1, N):
            D *= z[i] - z[j]
    return z, (D*D.conj()).real/factorial(N)          # Weyl weight, mean = 1

def h_polys(z, kmax):
    """complete homogeneous symmetric polynomials h_0..h_kmax of the eigenvalues"""
    p = [None] + [np.sum(z**k, axis=0) for k in range(1, kmax + 1)]
    h = [np.ones(z.shape[1], dtype=complex)]
    for k in range(1, kmax + 1):
        h.append(sum(p[i]*h[k - i] for i in range(1, k + 1))/k)
    return h

def schur(lam, h):
    """s_lambda = det(h_{lam_i - i + j}) (Jacobi-Trudi)"""
    l = len(lam)
    if l == 0: return h[0]
    n = h[0].shape[0]
    Mx = np.zeros((n, l, l), dtype=complex)
    for i in range(l):
        for j in range(l):
            k = lam[i] - i + j
            if 0 <= k < len(h): Mx[:, i, j] = h[k]
    return np.linalg.det(Mx)

def dim(lam, N):
    lam = list(lam) + [0]*(N - len(lam))
    num = den = 1
    for i in range(N):
        for j in range(i + 1, N):
            num *= lam[i] - lam[j] + j - i; den *= j - i
    return num//den

def conj(lam, N):
    lam = list(lam) + [0]*(N - len(lam))
    c = [lam[0] - lam[N - 1 - i] for i in range(N)]
    return tuple(x for x in c if x > 0)

def partitions_box(N, lmax):
    """partitions with at most N-1 rows and parts <= lmax (the irreps of SU(N) up to that size)"""
    out = []
    for lam in itertools.product(range(lmax + 1), repeat=N - 1):
        if all(lam[i] >= lam[i + 1] for i in range(N - 2)):
            out.append(tuple(x for x in lam if x > 0))
    return out

def lemma_bound(N, beta):
    """Lemma mg:lem:u2sun as stated (max of the four entries; Lambda^2 entry for N >= 4)"""
    R = expm1(beta) - beta - beta**2/2
    dS, dA, dL = N*(N + 1)/2, N*N - 1, N*(N - 1)/2
    ent = [(beta**2/(8*N*N) + R/dS)/dS, (beta**2/(4*N*N) + R/dA)/dA, R/N**2]
    if N >= 4: ent.append((beta**2/(8*N*N) + R/dL)/dL)
    return max(ent)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 28
    z, W = torus(N, M)
    chi = z.sum(axis=0); Re, Im = chi.real, chi.imag
    print(f"SU({N}), grid {M}^{N-1} = {z.shape[1]} points: mean Weyl weight = {W.mean():.15f};  <|chi|^2> = {(W*abs(chi)**2).mean():.12f},"
          f"  <chi^N> = {(W*chi**N).mean().real:.12f}")
    lmax = 4
    h = h_polys(z, lmax + N)
    irreps = partitions_box(N, lmax)
    S = {lam: schur(lam, h) for lam in irreps}
    # orthonormality check of the characters on the grid
    worst = 0.0
    for a in irreps:
        for b in irreps:
            v = (W*S[a].conj()*S[b]).mean()
            worst = max(worst, abs(v - (1 if a == b else 0)))
    print(f"{len(irreps)} irreps with parts <= {lmax}: max |<chi_a, chi_b> - delta_ab| = {worst:.2e}")
    # multiplicities at order 2: m_r(a,b) = int conj(chi_r) chi^a chibar^b, a + b = 2
    print("order-2 multiplicities m_r(2,0), m_r(1,1), m_r(0,2) (nonzero only):")
    for lam in irreps:
        ms = [(W*S[lam].conj()*chi**a*chi.conj()**b).mean().real for a, b in ((2, 0), (1, 1), (0, 2))]
        if max(abs(m) for m in ms) > 0.5:
            print(f"   r = {lam!s:14s} (dim {dim(lam, N)}, conj {conj(lam, N)}):  {[int(round(m)) for m in ms]}")
    rho, rhob = (1,), conj((1,), N)
    for beta in (1e-3, 1e-2, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0):
        w = np.exp(beta/N*Re); c0 = (W*w).mean()
        mR = (W*w*Re).mean()/c0; u = mR/N
        Up = 2*((W*w*Re**2).mean()/c0 - mR**2); Um = 2*(W*w*Im**2).mean()/c0
        cr = {lam: (W*S[lam].conj()*w).mean() for lam in irreps}
        ur = {lam: cr[lam].real/(dim(lam, N)*c0) for lam in irreps if lam}
        others = {lam: v for lam, v in ur.items() if lam not in (rho, rhob)}
        top = max(others, key=others.get)
        maxim = max(abs(cr[lam].imag) for lam in irreps)
        bnd = lemma_bound(N, beta)
        print(f"beta={beta:g}: c0={c0:.15f}  u={u:.15e} (u_rho from c_rho: {ur[rho]:.15e})  U+ - 1 = {Up-1:+.6e}  U- - 1 = {Um-1:+.6e}"
              f"  [max Im c_r {maxim:.1e}]")
        print(f"      u_2 = max over r not in {{0,rho,rhobar}} (parts <= {lmax}) = {others[top]:.6e} at r = {top} (dim {dim(top, N)});"
              f"  Lemma 16.11 bound = {bnd:.6e};  ratio bound/u_2 = {bnd/others[top]:.6f};  e^b-1-b = {expm1(beta)-beta:.3e}")
        if beta in (1e-2, 0.5):
            for lam in sorted(others, key=others.get, reverse=True)[:8]:
                print(f"         u_r = {others[lam]:.6e}  r = {lam} (dim {dim(lam, N)})")

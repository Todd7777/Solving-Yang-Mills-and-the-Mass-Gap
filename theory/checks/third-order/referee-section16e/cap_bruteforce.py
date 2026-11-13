# Referee check (Part IV, Lemma mg:lem:uppertails): the five-face cap by brute force.
#
# Claim: integrating the eight inner links of the five faces other than F of a unit cube, with f_beta = sum_{r != 0} d_r u_r chi_r on each
# face, gives a class function g of the holonomy V of dF with g = sum_{r != 0} d_r u_r^5 chi_r, hence |g| <= sigma_5 = sum d_r^2 u_r^5.
#
# Basis-free form used here: for class functions the convolution (a*b)(x) = int a(y) b(y^-1 x) dy has Fourier coefficients
# hat(a*b)_r = hat a_r hat b_r / d_r, so the claim is g = f*f*f*f*f (five-fold convolution), and sigma_5 = g(1).  We compute g by summing
# over ALL configurations of the eight inner links (no gauge fixing) for Z_3, Z_4 and S_3, with the bottom links arbitrary (not gauge
# fixed, several random choices per class), and with the four vertical links gauge fixed (an exact change of variables) for the
# Frobenius group F_21 (complex three-dimensional irreducible representations with trivial determinant, like SU(3)).  For Z_n and S_3 the
# result is also compared with sum_r d_r u_r^5 chi_r from the character tables, and sigma_5 with sum d_r^2 u_r^5.
import itertools, sys
import numpy as np

def group_Zn(n):
    M = np.array([[(a + b) % n for b in range(n)] for a in range(n)])
    return M

def group_S3():
    els = list(itertools.permutations(range(3)))
    idx = {p: i for i, p in enumerate(els)}
    comp = lambda p, q: tuple(p[q[i]] for i in range(3))   # (p q)(i) = p(q(i))
    M = np.array([[idx[comp(p, q)] for q in els] for p in els])
    return M, els

def group_F21():
    els = [(a, b) for b in range(3) for a in range(7)]
    idx = {g: i for i, g in enumerate(els)}
    mul = lambda g, h: ((g[0] + pow(2, g[1], 7) * h[0]) % 7, (g[1] + h[1]) % 3)
    M = np.array([[idx[mul(g, h)] for h in els] for g in els])
    return M, els

def inverse_table(M):
    n = len(M); e = [i for i in range(n) if all(M[i, j] == j for j in range(n))][0]
    inv = np.array([[j for j in range(n) if M[i, j] == e][0] for i in range(n)])
    return inv, e

def convolve(a, b, M, inv):
    n = len(M)
    return np.array([sum(a[y] * b[M[inv[y], x]] for y in range(n)) for x in range(n)]) / n

def cap_g(f, M, inv, bottom, gauge_fix_vertical=False):
    """g = average over the 8 inner links of prod_{5 faces} f(U_face); bottom = (b0, b1, b2, b3) the links of F around the square:
    vertices v00, v10, v11, v01 at height 0 and w.. at height 1; b0: v00->v10, b1: v10->v11, b2: v01->v11, b3: v00->v01."""
    n = len(M)
    mul = lambda x, y: M[x, y]
    I = lambda x: inv[x]
    ax = np.arange(n)
    if gauge_fix_vertical:
        e = [i for i in range(n) if all(M[i, j] == j for j in range(n))][0]
        shape = (n,) * 4
        t0, t1, t2, t3 = np.meshgrid(ax, ax, ax, ax, indexing='ij')
        z = np.full(shape, e)
        c00 = c10 = c11 = c01 = z
    else:
        grids = np.meshgrid(*([ax] * 8), indexing='ij')
        t0, t1, t2, t3, c00, c10, c11, c01 = grids
    b0, b1, b2, b3 = bottom
    # top face: w00->w10 (t0), w10->w11 (t1), w01->w11 (t2), w00->w01 (t3); holonomy t0 t1 t2^-1 t3^-1
    top = mul(mul(mul(t0, t1), I(t2)), I(t3))
    # vertical links c_v: v -> w.  side faces (holonomy around each, any orientation: f is a class function invariant under inversion)
    s0 = mul(mul(mul(b0, c10), I(t0)), I(c00))        # v00 -> v10 -> w10 -> w00 -> v00
    s1 = mul(mul(mul(b1, c11), I(t1)), I(c10))        # v10 -> v11 -> w11 -> w10
    s2 = mul(mul(mul(b2, c11), I(t2)), I(c01))        # v01 -> v11 -> w11 -> w01
    s3 = mul(mul(mul(b3, c01), I(t3)), I(c00))        # v00 -> v01 -> w01 -> w00
    val = f[top] * f[s0] * f[s1] * f[s2] * f[s3]
    return val.mean()

def holonomy_F(bottom, M, inv):
    b0, b1, b2, b3 = bottom
    return M[M[M[b0, b1], inv[b2]], inv[b3]]

def check(name, M, chi_rho, N, beta, char_table=None, gauge_fix_vertical=False, ntries=3, rng=np.random.default_rng(1)):
    n = len(M); inv, e = inverse_table(M)
    w = np.exp(beta / N * np.real(chi_rho))
    c0 = w.mean(); f = w / c0 - 1
    assert np.allclose(f, f[inv])
    g5 = f
    for _ in range(4): g5 = convolve(g5, f, M, inv)
    maxdev = 0.0; maxratio = 0.0
    for V in range(n):
        for t in range(ntries):
            # random bottom links with holonomy V: choose b0, b1, b2 at random, solve b3
            b0, b1, b2 = rng.integers(n, size=3)
            # V = b0 b1 b2^-1 b3^-1  =>  b3 = (b0 b1 b2^-1)^-1 V ... b3^-1 = (b0 b1 b2^-1)^-1 V
            x = M[M[b0, b1], inv[b2]]
            b3 = inv[M[inv[x], V]]
            assert holonomy_F((b0, b1, b2, b3), M, inv) == V
            g = cap_g(f, M, inv, (b0, b1, b2, b3), gauge_fix_vertical)
            maxdev = max(maxdev, abs(g - g5[V]))
            maxratio = max(maxratio, abs(g) / g5[e])
    line = f"{name:5s} beta={beta}: max|g_brute - f*f*f*f*f| = {maxdev:.2e}; sigma_5 = g(1) = {g5[e]:.6e}; max|g|/sigma_5 = {maxratio:.6f}"
    if char_table is not None:
        dims, chis = char_table
        cr = [np.mean(np.conj(ch) * w) for ch in chis]
        u = [np.real(c / (dr * c0)) for c, dr in zip(cr, dims)]
        gch = sum(dr * ur**5 * ch for dr, ur, ch in zip(dims, u, chis) if not np.allclose(ch, 1))
        s5 = sum(dr**2 * ur**5 for dr, ur, ch in zip(dims, u, chis) if not np.allclose(ch, 1))
        line += f"; max|f^*5 - sum d_r u_r^5 chi_r| = {np.max(np.abs(g5 - gch)):.2e}; |sigma_5 - sum d_r^2 u_r^5| = {abs(g5[e] - s5):.2e}"
        line += f"; u_r = {[round(x, 6) for x in u]}"
    print(line)
    return maxdev

if __name__ == "__main__":
    tol = 1e-12
    worst = 0.0
    for n in (3, 4):
        M = group_Zn(n)
        chis = [np.exp(2j * np.pi * k * np.arange(n) / n) for k in range(n)]
        for beta in (0.3, 1.5):
            worst = max(worst, check(f"Z_{n}", M, chis[1], 1, beta, ([1] * n, chis)))
    M, els = group_S3()
    sgn = lambda p: np.linalg.det(np.eye(3)[list(p)])
    fix = lambda p: sum(1 for i in range(3) if p[i] == i)
    chi_triv = np.ones(6); chi_sgn = np.array([sgn(p) for p in els]); chi_2 = np.array([fix(p) - 1.0 for p in els])
    for beta in (0.3, 1.5):
        worst = max(worst, check("S_3", M, chi_2, 2, beta, ([1, 1, 2], [chi_triv, chi_sgn, chi_2])))
    M, els = group_F21()
    om = np.exp(2j * np.pi / 7)
    chi3 = np.array([(om**a + om**(2 * a) + om**(4 * a)) if b == 0 else 0 for (a, b) in els])
    for beta in (0.3, 1.5):
        worst = max(worst, check("F_21", M, chi3, 3, beta, None, gauge_fix_vertical=True, ntries=2))
    print("all agree to", tol, ":", worst < tol)

# Written by the second independent referee of Section 15 of Part IV (2026-11-12) and added to the repository with credit; unchanged
# apart from this header. A compass, not part of the proofs.
"""Referee check of Lemma mg:lem:split at finite beta (a compass; independent of upsilon_split.py and haar_moments.py).

(a) The SU(N) generating function Z(alpha, gamma) = <exp(alpha chi + gamma conj(chi))> is evaluated numerically
    (mpmath, 40 digits) from the Toeplitz sum over the centre, sum_q det[fhat_{j-k+q}], with
    fhat_n = alpha^n g_n(alpha gamma), fhat_{-n} = gamma^n g_n(alpha gamma), g_n(w) = sum_m w^m/(m!(m+n)!).
    Derivatives in gamma come from Cauchy integrals (trapezoid rule on a circle, exponentially accurate), so
        Upsilon_+ - Upsilon_- = 2(<cb^2>_y - <cb>_y^2) = 2(Z_gg/Z - (Z_g/Z)^2) at alpha = gamma = beta/(2N).
    Upsilon_+ and Upsilon_- are also computed separately, as 2 Var_y(Re chi) and 2 Var_y(Im chi), via
    Var_y(Re chi) = (1/4)(Z_aa + 2Z_ag + Z_gg)/Z - ((Z_a+Z_g)/(2Z))^2 etc.
(b) The Toeplitz formula itself is validated against direct Weyl integration on the maximal torus for SU(2) (1-d) and SU(3) (2-d).
(c) The finite-beta values are compared with the exact series of split_toeplitz_exact.py: the relative difference must be
    of the order of the first omitted term.
(d) Monte Carlo compass: Haar sampling of SU(6) (QR of complex Ginibre matrices, phases fixed, determinant removed by an
    N-th root) at beta = 12, where Upsilon_+ - Upsilon_- is visible, against the Toeplitz value.
Run: python3 split_numeric_mp.py   (about two minutes)
"""
import mpmath as mp
import numpy as np
from math import factorial
import split_toeplitz_exact as ex

mp.mp.dps = 40


def g(n, w, tol=mp.mpf(10) ** -45):
    s = mp.mpf(0)
    term = mp.mpf(1) / mp.factorial(n)
    m = 0
    while True:
        s += term
        m += 1
        term = term * w / (m * (m + n))
        if abs(term) < tol * (1 + abs(s)) and m > 5:
            return s + term


def fhat(n, a, c):
    w = a * c
    return a ** n * g(n, w) if n >= 0 else c ** (-n) * g(-n, w)


def Z(N, a, c, qmax=None):
    if qmax is None:
        qmax = 16 if N <= 3 else 8
    tot = mp.mpc(0)
    cache = {}
    for q in range(-qmax, qmax + 1):
        M = mp.matrix(N, N)
        for j in range(N):
            for k in range(N):
                n = j - k + q
                if n not in cache:
                    cache[n] = fhat(n, a, c)
                M[j, k] = cache[n]
        tot += mp.det(M)
    return tot


def derivs(N, t, npts=16, r=None):
    """Z and its first and second derivatives in (alpha, gamma) at alpha = gamma = t, by 2-d Cauchy integrals."""
    if r is None:
        r = mp.mpf('0.5') * max(t, mp.mpf('0.2'))
    th = [2 * mp.pi * k / npts for k in range(npts)]
    vals = {}
    for i, ti in enumerate(th):
        for j, tj in enumerate(th):
            vals[i, j] = Z(N, t + r * mp.expj(ti), t + r * mp.expj(tj))
    def coef(p, q):  # coefficient of (da)^p (dg)^q in the Taylor expansion, times p! q!
        s = mp.mpc(0)
        for i, ti in enumerate(th):
            for j, tj in enumerate(th):
                s += vals[i, j] * mp.expj(-p * ti - q * tj)
        return (s / npts ** 2 / r ** (p + q) * mp.factorial(p) * mp.factorial(q))
    return {(p, q): coef(p, q) for (p, q) in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]}


def upsilons(N, beta):
    t = mp.mpf(beta) / (2 * N)
    D = derivs(N, t)
    Z0 = D[0, 0]
    # chi <-> d/dalpha, cb <-> d/dgamma ; Re chi = (chi+cb)/2, Im chi = (chi-cb)/(2i)
    m_re = (D[1, 0] + D[0, 1]) / (2 * Z0)
    m_re2 = (D[2, 0] + 2 * D[1, 1] + D[0, 2]) / (4 * Z0)
    m_im = (D[1, 0] - D[0, 1]) / (2j * Z0)
    m_im2 = -(D[2, 0] - 2 * D[1, 1] + D[0, 2]) / (4 * Z0)
    Up = 2 * (m_re2 - m_re ** 2)
    Um = 2 * (m_im2 - m_im ** 2)
    split = 2 * (D[0, 2] / Z0 - (D[0, 1] / Z0) ** 2)
    return mp.re(Up), mp.re(Um), mp.re(split), mp.im(Up), mp.im(Um)


def weyl_check():
    """Toeplitz vs Weyl integration: [<e^y>, <cb>_y, <cb^2>_y] for SU(2), SU(3), y = (beta/N) Re chi.
    Weyl: int_{SU(N)} f = (1/N!) int_T |Delta(z)|^2 f(z), T the maximal torus {prod z = 1} with normalised Haar measure;
    the integrands are smooth and periodic, so the trapezoid rule on a uniform grid is spectrally accurate."""
    out = []
    for N, beta in ((2, 3.0), (3, 3.0), (3, 0.7)):
        t = mp.mpf(beta) / (2 * N)
        D = derivs(N, t)
        tz = [mp.re(D[0, 0]), mp.re(D[0, 1] / D[0, 0]), mp.re(D[0, 2] / D[0, 0])]
        n = 256
        th = 2 * np.pi * np.arange(n) / n
        if N == 2:
            z = np.stack([np.exp(1j * th), np.exp(-1j * th)])
            V = np.abs(z[0] - z[1]) ** 2 / 2
            ch = z.sum(axis=0)
        else:
            A, B = np.meshgrid(th, th, indexing='ij')
            z = np.stack([np.exp(1j * A), np.exp(1j * B), np.exp(-1j * (A + B))])
            V = np.abs((z[0] - z[1]) * (z[0] - z[2]) * (z[1] - z[2])) ** 2 / 6
            ch = z.sum(axis=0)
        w = np.exp(beta / N * ch.real)
        I0 = (V * w).mean(); I1 = (V * np.conj(ch) * w).mean(); I2 = (V * np.conj(ch) ** 2 * w).mean()
        wy = [I0, (I1 / I0).real, (I2 / I0).real]
        out.append((N, beta, tz, wy, abs((I1 / I0).imag) + abs((I2 / I0).imag)))
    return out


def mc_su(N, beta, nsamp=10_000_000, batch=250_000, seed=12345):
    rng = np.random.default_rng(seed)
    t = beta / N
    sums = []
    for _ in range(nsamp // batch):
        Zg = (rng.standard_normal((batch, N, N)) + 1j * rng.standard_normal((batch, N, N))) / np.sqrt(2)
        Q, R = np.linalg.qr(Zg)
        d = np.diagonal(R, axis1=1, axis2=2)
        Q = Q * (d / np.abs(d))[:, None, :]
        det = np.linalg.det(Q)
        # a uniformly chosen N-th root of det^{-1}; any choice gives Haar on SU(N) (central ambiguity)
        k = rng.integers(0, N, size=batch)
        root = np.exp(-1j * (np.angle(det) + 2 * np.pi * k) / N)
        ch = np.trace(Q, axis1=1, axis2=2) * root
        w = np.exp(t * ch.real)
        sums.append((w.mean(), (np.conj(ch) * w).mean(), (np.conj(ch) ** 2 * w).mean(), (np.abs(ch) ** 6).mean(), (ch ** N).mean()))
    sums = np.array(sums)
    Ebatch = (2 * (sums[:, 2] / sums[:, 0] - (sums[:, 1] / sums[:, 0]) ** 2)).real
    W, C1, C2, A6, CN = sums.mean(axis=0)
    est = 2 * (C2 / W - (C1 / W) ** 2)
    err = Ebatch.std(ddof=1) / np.sqrt(len(Ebatch))
    return est.real, err, A6.real, CN


if __name__ == '__main__':
    print('(b) Toeplitz sum vs Weyl integration, [<e^y>, <cb>_y, <cb^2>_y]:')
    for N, beta, tz, wy, im in weyl_check():
        print('  SU(%d), beta=%g: Toeplitz %s' % (N, beta, [mp.nstr(x, 15) for x in tz]))
        print('                 Weyl     %s  (imaginary parts %.1e)' % (['%.14e' % x.real if isinstance(x, complex) else '%.14e' % x for x in wy], im))
        assert max(abs(float(a) - float(np.real(b))) / abs(float(a)) for a, b in zip(tz, wy) if abs(float(a)) > 1e-12) < 1e-11
    print()
    print('(a),(c) Upsilon_+ - Upsilon_- at finite beta vs the exact series (DEG=14):')
    ex.DEG = 14
    for N in (6, 7, 8):
        mom = ex.moments(N, 14)
        Up_s, Um_s, Voff_s = ex.ups_series(mom, N, 12)
        for beta in (2.0, 1.0, 0.5):
            Up, Um, split, iUp, iUm = upsilons(N, beta)
            ser = sum(2 * float(c) * beta ** n for n, c in enumerate(Voff_s))
            lead = 2 * beta ** (N - 2) / ((2 * N) ** (N - 2) * factorial(N - 2))
            print('  SU(%d) beta=%4.1f: Ups+ - Ups- = %s (from Ups+, Ups- separately: %s); series %.15e; leading %.6e; ratio to leading %.10f'
                  % (N, beta, mp.nstr(split, 16), mp.nstr(Up - Um, 16), ser, lead, float(split) / lead))
    print()
    print('(d) Monte Carlo compass, SU(6) at beta = 12 (alpha = gamma = 1):')
    Up, Um, split, _, _ = upsilons(6, 12.0)
    est, err, A6, CN = mc_su(6, 12.0)
    print('  Toeplitz: Ups+ - Ups- = %s (Ups+ = %s, Ups- = %s)' % (mp.nstr(split, 10), mp.nstr(Up, 10), mp.nstr(Um, 10)))
    print('  Monte Carlo (1e7 Haar samples): %.5f +- %.5f ; also <|chi|^6> = %.3f (exact 6), <chi^6> = %.3f%+.3fi (exact 1)'
          % (est, err, A6, CN.real, CN.imag))

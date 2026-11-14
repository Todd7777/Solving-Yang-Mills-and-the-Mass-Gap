# Referee (Section 16h): the Haar moments E[chi^a conj(chi)^b] of SU(N) by alias-free Weyl integration (exact integers after rounding), against the
# authors' Schur-Weyl formula (sun_group.moment, imported read-only), and the exact Taylor series of u, U_+, U_- and -log U_+- from these moments.
# Usage: python3 moments_check.py [N] [DEG]     (default N=5, DEG=19; the torus grid has M = DEG + 2N - 1 points per angle; ~2 min for N=5)
"""Alias-free: on the grid with M points per angle theta_1..theta_{N-1} (theta_N = -sum), a monomial prod z_j^{k_j} has frequency k_j - k_N in theta_j.
In chi^a chibar^b |Delta|^2 these frequencies are at most a + b + 2(N-1) in absolute value, so the trapezoidal rule is exact when M > a + b + 2(N-1).
The series: with t = beta/N, E_w[F] = E[F e^{t Re chi}]/E[e^{t Re chi}] and e^{t Re chi} = sum_{a,b} (t/2)^{a+b} chi^a chibar^b/(a! b!).  Written here
from the definitions (Section mg:sec:second: U_+ = Var_w(sqrt2 Re chi), U_- = Var_w(sqrt2 Im chi), u = E_w[Re chi]/N), with exact fractions."""
import os, sys
import numpy as np
from fractions import Fraction
from math import factorial
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..'))
from weyl_sun import torus

def weyl_moments(N, DEG):
    M = DEG + 2*(N - 1) + 1
    z, W = torus(N, M)
    chi = z.sum(axis=0); chib = chi.conj()
    P = [np.ones_like(chi)]
    for k in range(DEG): P.append(P[-1]*chi)
    Pb = [p.conj() for p in P]
    out, worst = {}, 0.0
    for a in range(DEG + 1):
        for b in range(DEG + 1 - a):
            v = (W*P[a]*Pb[b]).mean()
            r = int(round(v.real))
            worst = max(worst, abs(v.real - r), abs(v.imag))
            out[(a, b)] = r
    return out, worst, M

def series(N, mom, order):
    """exact Taylor coefficients (in beta) of c0, E[e^y Re chi], ... and of u, U_+, U_-, to the given order"""
    def G(p, q):                       # E[chi^p chibar^q e^{t Re chi}] as a series in t
        c = [Fraction(0)]*(order + 1)
        for a in range(order + 1):
            for b in range(order + 1 - a):
                m = mom.get((a + p, b + q))
                if m is None: raise ValueError('moment degree too small')
                if m: c[a + b] += Fraction(m, 2**(a + b)*factorial(a)*factorial(b))
        return [c[k]/Fraction(N)**k for k in range(order + 1)]          # t^k -> beta^k
    def mul(x, y): return [sum(x[i]*y[k - i] for i in range(k + 1)) for k in range(order + 1)]
    def inv(x):
        y = [Fraction(0)]*(order + 1); y[0] = 1/x[0]
        for k in range(1, order + 1): y[k] = -sum(x[i]*y[k - i] for i in range(1, k + 1))/x[0]
        return y
    def lin(*terms): return [sum(c*s[k] for c, s in terms) for k in range(order + 1)]
    G00, G10, G01, G20, G11, G02 = (G(*pq) for pq in ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)))
    iZ = inv(G00)
    mR = mul(lin((Fraction(1, 2), G10), (Fraction(1, 2), G01)), iZ)
    ER2 = mul(lin((Fraction(1, 4), G20), (Fraction(1, 2), G11), (Fraction(1, 4), G02)), iZ)
    EI2 = mul(lin((Fraction(-1, 4), G20), (Fraction(1, 2), G11), (Fraction(-1, 4), G02)), iZ)
    Up = lin((2, ER2), (-2, mul(mR, mR)))
    Um = lin((2, EI2))
    u = [c/N for c in mR]
    return u, Up, Um

def minus_log(x, order):
    """-log of a series with x[0] = 1"""
    d = [Fraction(0)] + x[1:]                       # x - 1
    out = [Fraction(0)]*(order + 1); powk = [Fraction(1)] + [Fraction(0)]*order
    for k in range(1, order + 1):
        powk = [sum(powk[i]*d[j - i] for i in range(j + 1)) for j in range(order + 1)]
        for j in range(order + 1): out[j] += Fraction((-1)**k, k)*powk[j]
    return out

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    DEG = int(sys.argv[2]) if len(sys.argv) > 2 else 19
    mom, worst, M = weyl_moments(N, DEG)
    from sun_group import moment                     # the authors' Schur-Weyl moments (read-only)
    bad = [(a, b, mom[(a, b)], moment(N, a, b)) for (a, b) in mom if mom[(a, b)] != moment(N, a, b)]
    nz = sum(1 for v in mom.values() if v)
    print(f"SU({N}): {len(mom)} moments with a + b <= {DEG} by Weyl integration on a {M}^{N-1} grid (max distance to an integer {worst:.1e});"
          f" {nz} non-zero; disagreements with sun_group.moment: {bad if bad else 'none'}")
    print("  E|chi|^(2k), k = 0..%d:" % (DEG//2), [mom[(k, k)] for k in range(DEG//2 + 1)], "; E chi^N =", mom.get((N, 0)))
    order = DEG - 2
    u, Up, Um = series(N, mom, order)
    print(f"  u = " + ' + '.join(f'({c}) b^{k}' for k, c in enumerate(u) if c and k <= 7))
    print(f"  U_+ - 1 = " + ' + '.join(f'({c}) b^{k}' for k, c in enumerate(Up) if c and k > 0 and k <= 9))
    print(f"  U_- - 1 = " + ' + '.join(f'({c}) b^{k}' for k, c in enumerate(Um) if c and k > 0 and k <= 9))
    for name, s in (('U_+', Up), ('U_-', Um)):
        ml = minus_log(s, min(order, 10))
        print(f"  -log {name} = " + ' + '.join(f'({c}) b^{k} [{float(c):.6e}]' for k, c in enumerate(ml) if c))
    # the alternating test: U_+ - U_- and the first coefficients
    diff = [a - b for a, b in zip(Up, Um)]
    print("  U_+ - U_- = " + ' + '.join(f'({c}) b^{k}' for k, c in enumerate(diff) if c and k <= 9))

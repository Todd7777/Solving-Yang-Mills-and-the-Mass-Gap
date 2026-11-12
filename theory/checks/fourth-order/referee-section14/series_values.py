# Written by the independent referee of Section 14 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs.
# Independent check of the values quoted in Corollary 14.3 / Conjecture 12.4: -log Upsilon_C as a series in u, and the u^4 coefficient
# 4log(1/u) - log Upsilon_C - h_C u^4 written as m + 4 log u. Moments of R = Re chi_rho from Weyl integration (rounded to exact integers),
# cumulant generating function K(t) = log <exp(t R)>; u = K'(t)/N, Upsilon_+ = 2K''(t) (iota=0) or K''(t) (iota=1).
import numpy as np, sympy as sp
from math import comb

def su3_moment(j, k, M=96):
    th = 2 * np.pi * np.arange(M) / M
    t1, t2 = np.meshgrid(th, th, indexing='ij'); t3 = -t1 - t2
    z = [np.exp(1j * t1), np.exp(1j * t2), np.exp(1j * t3)]
    D2 = abs((z[0] - z[1]) * (z[0] - z[2]) * (z[1] - z[2]))**2
    chi = z[0] + z[1] + z[2]
    val = (chi**j * np.conj(chi)**k * D2).mean() / 6
    r = round(val.real)
    assert abs(val - r) < 1e-8, (j, k, val)
    return r

def su2_moment(n, M=200):
    th = np.pi * (np.arange(M) + 0.5) / M          # Weyl: (2/pi) int_0^pi sin^2 th (2cos th)^n dth
    return round(((2 * np.cos(th))**n * np.sin(th)**2).mean() * 2)

t, uu = sp.symbols('t u')
ORDER = 9

def series_from_moments(mom):      # mom[n] = <R^n>
    Mt = sum(sp.Rational(mom[n]) * t**n / sp.factorial(n) for n in range(ORDER + 1))
    K = sp.series(sp.log(Mt), t, 0, ORDER + 1).removeO()
    return K

def invert(useries, order):
    # useries = a1 t + a2 t^2 + ... ; return t as series in u
    a = [sp.expand(useries).coeff(t, k) for k in range(order + 2)]
    tt = uu / a[1]
    for _ in range(order + 1):
        tt = sp.series(tt - (useries.subs(t, tt) - uu) / a[1], uu, 0, order + 2).removeO()
    return sp.expand(tt)

def report(name, K, N, iota, sector, h):
    Kp = sp.diff(K, t); Kpp = sp.diff(K, t, 2)
    useries = sp.series(Kp / N, t, 0, ORDER).removeO()
    tt = invert(useries, 6)
    if sector == '+': Ups = (2 if iota == 0 else 1) * Kpp
    else: Ups = None
    return useries, tt

# SU(3): R = Re chi, t = beta/3, u = <chi e^{tR}>/(3<e^{tR}>) = K'(t)/3 ; Upsilon_+ = 2 K''(t)
mom3 = []
for n in range(ORDER + 1):
    s = sum(comb(n, j) * su3_moment(j, n - j) for j in range(n + 1))
    mom3.append(sp.Rational(s, 2**n))
K3 = series_from_moments(mom3)
u3 = sp.series(sp.diff(K3, t) / 3, t, 0, ORDER).removeO()
t3 = invert(u3, 6)
Ups3 = sp.series(2 * sp.diff(K3, t, 2), t, 0, 7).removeO()
mlu3 = sp.series(-sp.log(Ups3.subs(t, t3)), uu, 0, 6).removeO()
print('SU(3): moments <R^n> =', mom3)
print('SU(3): -log Upsilon_+ in u :', sp.expand(mlu3))
for d in (3, 4):
    print(f'   d={d}: m + 4 log u = ', sp.expand(mlu3 - (26 * d - 70) * uu**4))

# SU(2): R = chi (real), w = exp(beta/2 chi), t = beta/2, u = K'(t)/2, Upsilon = K''(t)
mom2 = [sp.Rational(su2_moment(n)) for n in range(ORDER + 1)]
K2 = series_from_moments(mom2)
u2 = sp.series(sp.diff(K2, t) / 2, t, 0, ORDER).removeO()
t2 = invert(u2, 6)
Ups2 = sp.series(sp.diff(K2, t, 2), t, 0, 7).removeO()
mlu2 = sp.series(-sp.log(Ups2.subs(t, t2)), uu, 0, 6).removeO()
print('SU(2): -log Upsilon in u :', sp.expand(mlu2))
for d in (3, 4):
    print(f'   d={d}: m + 4 log u = ', sp.expand(mlu2 - (26 * d - 70) * uu**4))
# in beta for d=4: m - 4 log(4/beta) = -4 log(4u/beta) - log Upsilon - 34 u^4, beta = 2t
b = sp.symbols('beta')
u2b = sp.expand(u2.subs(t, b / 2))
expr = sp.series(-4 * sp.log(4 * u2b / b) - sp.log(Ups2.subs(t, b / 2)) - 34 * u2b**4, b, 0, 6).removeO()
print('SU(2), d=4: m - 4log(4/beta) =', sp.expand(expr))

# U(1): R = cos(theta), <R^n> = binom(n, n/2)/2^n for even n; t = beta; u = K'(t); Upsilon_- = 2(1 - M''/M) = 2(1 - K'' - K'^2)
mom1 = [sp.Rational(comb(n, n // 2), 2**n) if n % 2 == 0 else 0 for n in range(ORDER + 1)]
K1 = series_from_moments(mom1)
u1 = sp.series(sp.diff(K1, t), t, 0, ORDER).removeO()
t1 = invert(u1, 6)
Upsm = sp.series(2 * (1 - sp.diff(K1, t, 2) - sp.diff(K1, t)**2), t, 0, 7).removeO()
mlu1 = sp.series(-sp.log(Upsm.subs(t, t1)), uu, 0, 6).removeO()
print('U(1): -log Upsilon_- in u :', sp.expand(mlu1))
for d in (3, 4):
    print(f'   d={d}: m + 4 log u (C-, h_-={10*d-22}) = ', sp.expand(mlu1 - (10 * d - 22) * uu**4),
          f';  with h_+={26*d-70}:', sp.expand(mlu1 - (26 * d - 70) * uu**4))

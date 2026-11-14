# Referee (Section 16h): sun_group.sun(5) (the authors' group data, imported read-only) against (a) direct Weyl integration of the weighted
# integrals in double precision (weyl_sun.torus), and (b) the referee's exact Taylor series from alias-free Weyl moments (moments_check.series),
# evaluated at 60 digits; and the u_2 of sun_group.u2_bound against Lemma 16.11 as stated and against the true u_Lambda^2 (Weyl, double).
# Usage: python3 compare_group_data.py        (about 1 minute)
import os, sys
import numpy as np
from mpmath import mp, mpf, fsum
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..')); sys.path.insert(0, os.path.join(here, '..', '..', 'upper-bound'))
from weyl_sun import torus, lemma_bound
from moments_check import weyl_moments, series
from sun_group import sun, u2_bound
from constants_upper import beta0
mp.dps = 60
N = 5
mom, worst, M = weyl_moments(N, 19)
u_s, Up_s, Um_s = series(N, mom, 17)                       # exact to beta^17
ev = lambda c, b: fsum(mpf(x.numerator)/x.denominator*b**k for k, x in enumerate(c))
g = sun(N)
z, W = torus(N, 28); chi = z.sum(axis=0); Re, Im = chi.real, chi.imag
print("(a) double-precision Weyl integration vs sun_group.sun(5):")
for beta in (1e-3, 1e-2, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0):
    w = np.exp(beta/N*Re); c0 = (W*w).mean(); mR = (W*w*Re).mean()/c0
    Up = 2*((W*w*Re**2).mean()/c0 - mR**2); Um = 2*(W*w*Im**2).mean()/c0
    d = g(mpf(beta))
    print(f"  beta={beta:g}: |u/u_auth - 1| = {abs(mR/N/float(d['u']) - 1):.1e};  |U+ - U+_auth| = {abs(Up - float(d['Upm'][0])):.1e};"
          f"  |U- - U-_auth| = {abs(Um - float(d['Upm'][1])):.1e};  U+ - 1 = {Up - 1:.6e} (auth {float(d['Upm'][0] - 1):.6e});  Ups = U+: {d['Ups'] == d['Upm'][0]}")
print("(b) exact series (degree 17, referee) vs sun_group.sun(5) (degree 40), 60 digits:")
B3, B4 = beta0(3), beta0(4)
for beta in (mpf('1e-12'), mpf('1e-8'), mpf('1e-6'), mpf('1e-5'), mpf('1e-4'), B4, mpf('2e-4'), B3, mpf('1e-3'), mpf('1e-2')):
    d = g(beta)
    us, Ups, Ums = ev(u_s, beta), ev(Up_s, beta), ev(Um_s, beta)
    print(f"  beta={mp.nstr(beta, 6)}: rel.err u {mp.nstr(abs(us/d['u'] - 1), 2)};  (U+ - 1) series {mp.nstr(Ups - 1, 12)}, auth {mp.nstr(d['Upm'][0] - 1, 12)},"
          f" rel.diff {mp.nstr(abs((Ups - 1)/(d['Upm'][0] - 1) - 1), 2)};  rel.diff (U- - 1) {mp.nstr(abs((Ums - 1)/(d['Upm'][1] - 1) - 1), 2)};"
          f"  (U+ - 1)/(beta^3/6000) = {mp.nstr((Ups - 1)/(beta**3/6000), 12)}")
print("(c) sun_group.u2_bound(5, beta) vs Lemma 16.11 as stated (referee's transcription) and vs beta^2/2000 + R/100:")
for beta in (mpf('1e-6'), B4, B3, mpf('1e-3'), mpf('0.05'), mpf('0.5'), mpf('2')):
    a = u2_bound(5, beta); b = lemma_bound(5, float(beta)); R = mp.expm1(beta) - beta - beta**2/2
    print(f"  beta={mp.nstr(beta, 6)}: u2_bound = {mp.nstr(a, 12)}, lemma = {b:.12e}, beta^2/2000 + R/100 = {mp.nstr(beta**2/2000 + R/100, 12)}")

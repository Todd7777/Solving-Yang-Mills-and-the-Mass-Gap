# lmp_fractions.py -- written by the independent referee of Section 12 of Part IV (2026-11-12), added to the repository with credit;
# unchanged apart from this header. A compass, not part of the proofs.
# Coefficients printed in LMP07 (arXiv:0710.0512v1 = v2), eqs. (2.14)-(2.18), transcribed from pdftotext -raw output.
from math import gcd
from fractions import Fraction as Fr
eqs = {
 '2.14 SU2 A1': [(2,1),(98,3),(20984,405),(151496,243)],
 '2.15 SU2 E' : [(2,1),(26,3),(13036,405),(28052,243)],
 '2.16 SU3 A1': [(3,1),(9,1),(27,2),(7,1),(297,2),(858827,10240),(47641149,71680)],
 '2.17 SU3 E' : [(3,1),(9,1),(27,2),(17,1),(153,2),(1104587,10240),(29577789,71680)],
 '2.18 SU3 T1': [(3,1),(9,2),(98,4),(33,4),(36771,1280),(117897,448)],
}
tot = 0
for k, cs in eqs.items():
    for a, b in cs:
        tot += 1
        if gcd(a, b) != 1: print('unreduced:', k, f'{a}/{b}')
print('number of printed coefficients:', tot)
# hopping part = (u^4 coefficient) - (u^4 coefficient of -log Upsilon_C), with -log Ups: SU2 4/3, SU3 C+ 27, SU3 C- -27/4
for k, c4, ml in (('SU2 A1', Fr(-98,3), Fr(4,3)), ('SU2 E', Fr(-26,3), Fr(4,3)), ('SU3 A1', Fr(-7), Fr(27)), ('SU3 E', Fr(17), Fr(27)),
                  ('SU3 T1 printed', Fr(-98,4), Fr(-27,4)), ('SU3 T1 if -99/4', Fr(-99,4), Fr(-27,4))):
    print(f'{k:16s}: u^4 coefficient minus that of -log Upsilon_C = {c4 - ml}')
# sensitivity of the u^4 coefficients of -log Ups_-+ to v4, w4 (v = u_6 = 3/4 u^2 + v3 u^3 + v4 u^4, w = u_8 = ...):
import itertools
def u4coef(sign, v, w):   # sign=-1: b = 1-3u-6v+8w (C-); sign=+1: c = 1+3u+6v+8w-18u^2 (C+); v,w = [v2,v3,v4]
    # x = series of (b-1) or (c-1) up to u^4
    x = [Fr(0)]*5
    x[1] = sign*3
    x[2] = sign*6*v[0] + 8*w[0] - (18 if sign > 0 else 0)
    x[3] = sign*6*v[1] + 8*w[1]
    x[4] = sign*6*v[2] + 8*w[2]
    # -log(1+x) up to u^4
    def mul(p, q): 
        r = [Fr(0)]*5
        for i in range(5):
            for j in range(5-i): r[i+j] += p[i]*q[j]
        return r
    x2 = mul(x, x); x3 = mul(x2, x); x4 = mul(x3, x)
    return -(x[4] - x2[4]/2 + x3[4]/3 - x4[4]/4)
v = [Fr(3,4), Fr(0), Fr(-9,16)]; w = [Fr(9,8), Fr(-9,8), Fr(81,32)]
print('u^4 coefficient of -log Ups_+:', u4coef(+1, v, w), '   of -log Ups_-:', u4coef(-1, v, w))
dv = [Fr(0), Fr(0), Fr(1)]; dw = dv
print('d/dv4: C+', u4coef(+1, [a+b for a,b in zip(v,dv)], w) - u4coef(+1, v, w), ' C-', u4coef(-1, [a+b for a,b in zip(v,dv)], w) - u4coef(-1, v, w))
print('d/dw4: C+', u4coef(+1, v, [a+b for a,b in zip(w,dw)]) - u4coef(+1, v, w), ' C-', u4coef(-1, v, [a+b for a,b in zip(w,dw)]) - u4coef(-1, v, w))

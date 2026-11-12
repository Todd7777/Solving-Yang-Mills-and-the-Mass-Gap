"""Upsilon_+ - Upsilon_- for SU(N), fundamental representation, as a power series in beta (a compass for Lemma mg:lem:split).
Exact moments <chi^a conj(chi)^b>: zero unless a = b (mod N), k! on the diagonal for k <= 3, <chi^N> = 1; these suffice up to the
orders computed (degree <= 7). Output: N=6: beta^4/248832; N=7: beta^5/32269440; N=8: 0 through beta^5.
The moments are checked by Haar sampling in haar_moments.py."""
import sympy as sp
from math import factorial
b = sp.symbols('beta')
def m(a, c, N):
    if (a - c) % N: return 0
    if a == c: return factorial(a) if a <= N else None
    if {a, c} == {N, 0}: return 1
    return None  # not needed below the order we use
def avg(expr, N):
    x, xb = sp.symbols('x xb')
    p = sp.Poly(sp.expand(expr), x, xb)
    tot = 0
    for (i, j), co in p.terms():
        v = m(i, j, N)
        assert v is not None, (i, j)
        tot += co * v
    return sp.nsimplify(tot)
def split(N, order):
    x, xb = sp.symbols('x xb')
    re = (x + xb) / 2; im = (x - xb) / (2 * sp.I)
    y = b / N * re
    ey = sum(y**j / factorial(j) for j in range(order + 1))
    C = avg(ey, N)
    out = {}
    for name, psi in (('+', sp.sqrt(2) * re), ('-', sp.sqrt(2) * im)):
        A = avg(sp.expand(psi**2 * ey), N); B = avg(sp.expand(psi * ey), N)
        U = sp.series(A / C - B**2 / C**2, b, 0, order + 1).removeO()
        out[name] = sp.expand(U)
    return out
for N in (6, 7, 8):
    o = split(N, 5 if N >= 7 else 4)
    print(N, 'Ups+ - Ups- =', sp.simplify(o['+'] - o['-']), '  Ups- =', o['-'])

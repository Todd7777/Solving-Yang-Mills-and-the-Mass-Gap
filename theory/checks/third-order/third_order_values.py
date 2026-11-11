# Section 11 of Part IV: the third-order coefficients of -log Upsilon(beta) for U(1) and SU(N) (fundamental representation).
#  (1) The moments <chi^j chibar^k> (j+k <= 5) by Weyl integration (trapezoidal rule on the maximal torus, exact for these
#      trigonometric polynomials), rounded to integers, and r_n = <(Re chi)^n>, m_n = <(Re chi)^n (Im chi)^2>.
#  (2) The expansion of Upsilon_+ = Var(sqrt2 Re chi) and Upsilon_- = Var(sqrt2 Im chi) under w = exp(t Re chi)/c0, t = beta/N, to order t^3,
#      in exact rational arithmetic:  Upsilon_+ = 1 + 2 r3 t + (r4 - 3/4) t^2 + (r5 - 5 r3) t^3/3,  Upsilon_- = 1 - 2 r3 t + (3/4 - r4) t^2 + (m3 + r3) t^3/3.
#  (3) A numerical check: -log Upsilon at small beta by Weyl integration against the series.
import numpy as np, math, itertools
from fractions import Fraction as Fr

def torus(N, M):
    th = np.linspace(0, 2*np.pi, M, endpoint=False)
    g = np.meshgrid(*([th]*(N-1)), indexing='ij'); ang = [x.ravel() for x in g]; ang.append(-sum(ang))
    z = np.exp(1j*np.array(ang)); vd = np.ones(z.shape[1])
    for i in range(N):
        for j in range(i+1, N): vd *= np.abs(z[i]-z[j])**2
    return z.sum(0), vd/(math.factorial(N)*M**(N-1))

GRID = {2: 64, 3: 48, 4: 24, 5: 16, 6: 15}      # > 2N+2 points per angle: no aliasing for degree <= 5

def moments(N):
    X, w = torus(N, GRID[N])
    mom = {}
    for j in range(6):
        for k in range(6 - j):
            v = np.sum(w*X**j*np.conj(X)**k)
            assert abs(v.imag) < 1e-9 and abs(v.real - round(v.real)) < 1e-9, (N, j, k, v)
            mom[(j, k)] = int(round(v.real))
    return mom

def rn(mom, n):        # <(Re chi)^n> = 2^-n sum_j C(n,j) <chi^j chibar^(n-j)>
    return Fr(sum(math.comb(n, j)*mom[(j, n-j)] for j in range(n+1)), 2**n)

def mn(mom, n):        # <(Re chi)^n (Im chi)^2> = <(Re chi)^n |chi|^2> - r_{n+2}
    s = Fr(sum(math.comb(n, j)*mom[(j+1, n-j+1)] for j in range(n+1)), 2**n)
    return s - rn(mom, n+2)

def series(N):
    mom = moments(N); r = {n: rn(mom, n) for n in range(6)}; m3 = mn(mom, 3)
    assert r[1] == 0 and r[2] == Fr(1, 2)
    t = Fr(1, N)
    up = [Fr(1), 2*r[3]*t, (r[4]-Fr(3, 4))*t**2, (r[5]-5*r[3])*t**3/3]          # coefficients of beta^0..3
    um = [Fr(1), -2*r[3]*t, (Fr(3, 4)-r[4])*t**2, (m3+r[3])*t**3/3]
    return mom, r, m3, up, um

def mlog(c):          # -log(c0 + c1 b + c2 b^2 + c3 b^3) to order 3, c0 = 1
    _, c1, c2, c3 = c
    return [Fr(0), -c1, -c2 + c1**2/2, -c3 + c1*c2 - c1**3/3]

def numeric(N, beta, sector):
    X, w = torus(N, GRID[N])
    W = np.exp(beta/N*X.real); p = w*W/np.sum(w*W)
    phi = np.sqrt(2)*(X.real if sector == '+' else X.imag)
    return -np.log(np.sum(p*phi**2) - np.sum(p*phi)**2)

if __name__ == "__main__":
    for N in (3, 4, 5, 6):
        mom, r, m3, up, um = series(N)
        print(f"SU({N}): <chi^j chibar^k> (j+k<=5, nonzero, j!=k): {[(jk, v) for jk, v in mom.items() if v and jk[0] != jk[1]]};  r3..r5 = {r[3]}, {r[4]}, {r[5]};  m3 = {m3}")
        for name, c in (('+', up), ('-', um)):
            L = mlog(c)
            # numerical check of the beta^3 coefficient: (-log Ups - L1 b - L2 b^2)/b^3, Richardson in b
            vals = [(float(numeric(N, b, name)) - float(L[1])*b - float(L[2])*b**2)/b**3 for b in (0.02, 0.01)]
            rich = 2*vals[1] - vals[0]
            print(f"   C{name}: Upsilon = 1 + ({c[1]}) b + ({c[2]}) b^2 + ({c[3]}) b^3;  -log Upsilon = ({L[1]}) b + ({L[2]}) b^2 + ({L[3]}) b^3"
                  f"   [numerical beta^3 coefficient {rich:+.6f}, exact {float(L[3]):+.6f}]")
    # U(1) and SU(2): Upsilon is even in beta (rho(-1) = -1), so no beta^3 term; numerical confirmation
    from mpmath import besseli, log, mpf
    for b in (mpf('0.02'), mpf('0.01')):
        um = 1 - besseli(2, b)/besseli(0, b); su2 = 1 + 3*besseli(3, b)/besseli(1, b) - 4*(besseli(2, b)/besseli(1, b))**2
        print(f"U(1) C-: (-log Ups - b^2/8)/b^3 = {float((-log(um) - b**2/8)/b**3):+.2e};   SU(2): (-log Ups - b^2/8)/b^3 = {float((-log(su2) - b**2/8)/b**3):+.2e}")

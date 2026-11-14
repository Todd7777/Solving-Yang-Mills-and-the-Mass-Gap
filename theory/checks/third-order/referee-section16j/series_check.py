# Referee check (tenth pass on Section 16, Corollary mg:cor:sectorscap), item 4: the exact series of -log U_+-, log(U_+/U_-) and of Q, from Haar moments of SU(N)
# computed by a method independent of the authors' (sun_group.py uses Schur-Weyl with the hook-length formula):
#   int_{SU(N)} exp(x chi + y conj(chi)) dU = sum_{k in Z} det[ c^(k)_{i-j}(x,y) ]_{i,j=1..N},   c^(k)_n(x,y) = sum_{a-b=n-k} x^a y^b/(a! b!)
# (Heine's identity for the Toeplitz determinant of the symbol z^k exp(xz + y/z) on U(N), and int_{SU(N)} F = sum_k int_{U(N)} F det^k for functions F
# on U(N)). The coefficient of x^a y^b is m(a,b)/(a! b!), m(a,b) = int chi^a conj(chi)^b dU. Everything in exact rational arithmetic, truncated at
# total degree DEG (exact below it).
# Then, with t = beta/N and the weight e^{t Re chi}: Z, Z_1, Z_2, Z_im as in the corollary's proof, U_+ = 2(Z_2 Z - Z_1^2)/Z^2, U_- = 2 Z_im/Z, Q = Z_2 Z - Z_1^2 - Z_im Z.
# Checked: (i) the moments against a few known values and against the authors' sun_group.moment; (ii) the series of (i), (ii) of the corollary against the
# authors' sector_series; (iii) the leading coefficient of log(U_+/U_-) against 2 (beta/2N)^(N-2)/(N-2)! for N = 3..8; (iv) the leading coefficient
# q_N of Q (1/2, 1/8, 1/48 for N = 3, 4, 5) and the exact coefficients of Q against the authors' (split_lower), to degree DEG - 2.
# Usage: python3 series_check.py
import os, sys, itertools
from fractions import Fraction as Fr
from math import factorial, comb
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..'))

def poly_mul(A, B, D):
    out = {}
    for (a1, b1), c1 in A.items():
        for (a2, b2), c2 in B.items():
            if a1 + a2 + b1 + b2 <= D:
                k = (a1 + a2, b1 + b2); out[k] = out.get(k, 0) + c1*c2
    return {k: v for k, v in out.items() if v}

def entry(n, D):
    """c_n(x,y) = sum_{a-b=n} x^a y^b/(a! b!), truncated at total degree D"""
    out = {}
    for b in range(D + 1):
        a = b + n
        if a < 0 or a + b > D: continue
        out[(a, b)] = Fr(1, factorial(a)*factorial(b))
    return out

def det_poly(M, D):
    """determinant of a matrix of truncated bivariate polynomials (Leibniz; N <= 8)"""
    N = len(M); total = {}
    for perm in itertools.permutations(range(N)):
        sign = 1
        for i in range(N):
            for j in range(i + 1, N):
                if perm[i] > perm[j]: sign = -sign
        prod = {(0, 0): Fr(1)}
        for i in range(N):
            prod = poly_mul(prod, M[i][perm[i]], D)
            if not prod: break
        for k, v in prod.items(): total[k] = total.get(k, 0) + sign*v
    return {k: v for k, v in total.items() if v}

def haar_moments(N, D):
    """m(a,b) = int_{SU(N)} chi^a conj(chi)^b dU for a + b <= D, from the Toeplitz determinants"""
    G = {}
    for k in range(-(D//N) - 1, D//N + 2):
        M = [[entry(k + i - j, D) for j in range(N)] for i in range(N)]
        for key, v in det_poly(M, D).items(): G[key] = G.get(key, 0) + v
    return {(a, b): v*factorial(a)*factorial(b) for (a, b), v in G.items()}

def series_t(m, D, g):
    """coefficients in t of int g(chi, chibar) e^{t Re chi} dU, g given as {(p,q): coefficient} (a polynomial in chi, chibar), to degree D - deg g"""
    dg = max(p + q for p, q in g)
    out = []
    for k in range(D - dg + 1):
        s = Fr(0)
        for j in range(k + 1):                 # (Re chi)^k = 2^-k sum_j C(k,j) chi^j chibar^(k-j)
            for (p, q), c in g.items():
                s += c*comb(k, j)*m.get((j + p, k - j + q), 0)
        out.append(s/(2**k*factorial(k)))
    return out

def mul(A, B, K): return [sum(A[i]*B[k - i] for i in range(k + 1)) for k in range(K + 1)]
def inv(A, K):
    out = [Fr(1)/A[0]]
    for k in range(1, K + 1): out.append(-sum(A[i]*out[k - i] for i in range(1, k + 1))/A[0])
    return out
def log1(A, K):          # log A for A(0) = 1, via log A = int A'/A
    dA = [(k + 1)*A[k + 1] for k in range(K)]
    q = mul(dA, inv(A, K), K - 1)
    return [Fr(0)] + [q[k - 1]/k for k in range(1, K + 1)]

def sector(N, D):
    m = haar_moments(N, D)
    K = D - 2
    Z = series_t(m, D, {(0, 0): Fr(1)})[:K + 1]
    Z1 = series_t(m, D, {(1, 0): Fr(1, 2), (0, 1): Fr(1, 2)})[:K + 1]
    Z2 = series_t(m, D, {(2, 0): Fr(1, 4), (1, 1): Fr(1, 2), (0, 2): Fr(1, 4)})[:K + 1]
    Zi = series_t(m, D, {(2, 0): Fr(-1, 4), (1, 1): Fr(1, 2), (0, 2): Fr(-1, 4)})[:K + 1]
    iZ = inv(Z, K)
    Up = [2*c for c in mul([a - b for a, b in zip(mul(Z2, Z, K), mul(Z1, Z1, K))], mul(iZ, iZ, K), K)]
    Um = [2*c for c in mul(Zi, iZ, K)]
    Q = [a - b - c for a, b, c in zip(mul(Z2, Z, K), mul(Z1, Z1, K), mul(Zi, Z, K))]
    tob = lambda L: [c/Fr(N)**k for k, c in enumerate(L)]
    lp, lm = log1(Up, K), log1(Um, K)
    return dict(m=m, K=K, Z=Z, Z1=Z1, Z2=Z2, Zi=Zi, Up=tob(Up), Um=tob(Um), Q=Q, neglog_plus=tob([-c for c in lp]), neglog_minus=tob([-c for c in lm]),
                log_ratio=tob([a - b for a, b in zip(lp, lm)]))

def fmt(L, kmax):
    return ' + '.join(f"({c}) b^{k}" for k, c in enumerate(L[:kmax + 1]) if c != 0)

if __name__ == "__main__":
    import rigorous_bounds as RB, sun_group as SG
    print("Haar moments of SU(N) from Toeplitz determinants (exact):")
    for N in (3, 4, 5):
        m = haar_moments(N, 10)
        known = {(1, 1): 1, (2, 2): 2, (3, 3): 6, (N, 0): 1, (0, N): 1}
        okk = all(m.get(k, 0) == v for k, v in known.items())
        agree = all(m.get((a, b), 0) == SG.moment(N, a, b) for a in range(11) for b in range(11 - a))
        nz = sorted((a, b) for (a, b), v in m.items() if v and a + b <= 8)
        print(f"  SU({N}): m(1,1)=1, m(2,2)=2, m(3,3)=6, m(N,0)=m(0,N)=1: {okk};  all moments with a+b <= 10 equal to sun_group.moment: {agree};"
              f"  e.g. m(4,4) = {m[(4, 4)]}, m(5,5) = {m.get((5, 5))}, m({N + 1},1) = {m.get((N + 1, 1))}, m({2*N},0) = {m.get((2*N, 0))}")
    print()
    DEG = 12
    for N in (3, 4, 5):
        s = sector(N, DEG); A = RB.sector_series(N, K=6)
        same = all(s[key][k] == A[key][k] for key in ('neglog_plus', 'neglog_minus', 'log_ratio') for k in range(7))
        print(f"SU({N}) (moments to total degree {DEG}; series exact to beta^{s['K']}):")
        print(f"  -log U_+     = {fmt(s['neglog_plus'], 7)} + ...")
        print(f"  -log U_-     = {fmt(s['neglog_minus'], 7)} + ...")
        print(f"  log(U_+/U_-) = {fmt(s['log_ratio'], 7)} + ...")
        print(f"  the authors' sector_series agrees to beta^6: {same}")
        # Q exactly, against the authors' split_lower coefficients (recomputed here from sun_group with their formulas)
        C = SG.coeffs(N, s['K'])
        Za = C[(0, 0)]; Z1a = [(x + y)/2 for x, y in zip(C[(1, 0)], C[(0, 1)])]
        Z2a = [(x + 2*y + z)/4 for x, y, z in zip(C[(2, 0)], C[(1, 1)], C[(0, 2)])]; Zia = [-(x - 2*y + z)/4 for x, y, z in zip(C[(2, 0)], C[(1, 1)], C[(0, 2)])]
        Qa = [a - b - c for a, b, c in zip(mul(Z2a, Za, s['K']), mul(Z1a, Z1a, s['K']), mul(Zia, Za, s['K']))]
        k0 = next(k for k, q in enumerate(s['Q']) if q)
        print(f"  Q = Z_2 Z - Z_1^2 - Z_im Z: first non-zero coefficient q_{k0} = {s['Q'][k0]} (t^{k0}, t = beta/N); coefficients equal to the authors' to t^{s['K']}: {Qa == s['Q']}")
        print(f"  Q coefficients (t^k): " + ', '.join(f"{c}" for c in s['Q'][:s['K'] + 1]))
        # signs of Z_im's coefficients (the proof says they may have either sign)
        print(f"  Z_im coefficients (t^k): " + ', '.join(f"{c}" for c in s['Zi'][:8]) + "  (signs: " + ''.join('+' if c > 0 else ('-' if c < 0 else '0') for c in s['Zi'][:s['K'] + 1]) + ")")
        print()
    print("Leading term of log(U_+/U_-) against 2 (beta/2N)^(N-2)/(N-2)!:")
    for N in range(3, 9):
        D = N + 2
        s = sector(N, D)
        lr = s['log_ratio']; k0 = next(k for k, c in enumerate(lr) if c)
        pred = 2*Fr(1, (2*N)**(N - 2)*factorial(N - 2))
        print(f"  SU({N}): first non-zero term ({lr[k0]}) beta^{k0}; predicted ({pred}) beta^{N - 2}: {k0 == N - 2 and lr[k0] == pred};"
              f"  q_N = {s['Q'][N - 2]} = 2^-(N-2)/(N-2)! = {Fr(1, 2**(N - 2)*factorial(N - 2))}: {s['Q'][N - 2] == Fr(1, 2**(N - 2)*factorial(N - 2))}")

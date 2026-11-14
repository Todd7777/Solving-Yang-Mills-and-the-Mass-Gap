# Referee check (ninth pass on Section 16) of Lemma mg:lem:upsbounds(ii): U_+ >= 1 for SU(5) on (0, beta_max].
# Independent of the authors' code: the Taylor coefficients of Z(t) = int_{SU(N)} exp(t Re tr U) dU come from the Bessel-determinant formula
#     Z(t) = sum_{k in Z} det[ I_{k+i-j}(t) ]_{i,j=1..N}          (U(N) Haar measure restricted to det U = 1, by Fourier in the phase of det U),
# not from Schur-Weyl moments. Z_1 = Z', Z_2 = Z'' (differentiate under the integral), so P := 2 Z_2 Z - 2 Z_1^2 - Z^2 = (U_+ - 1) Z^2.
# Everything below is exact rational arithmetic (fractions.Fraction): no floating point enters the proof of P > 0.
# Usage: python3 ups_ge1_check.py [N [ORDER [BETA_MAX_NUM BETA_MAX_DEN]]]     (default N = 5, ORDER = 34, beta_max = 403/1000000)
import sys
from fractions import Fraction as Fr
from math import factorial
from itertools import permutations

N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
DEG = int(sys.argv[2]) if len(sys.argv) > 2 else 34
BMAX = Fr(int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else Fr(403, 1000000)

def bessel_I(n, deg):
    """Taylor coefficients of I_n(t) = sum_m (t/2)^(2m+|n|)/(m! (m+|n|)!), truncated at degree deg"""
    n = abs(n); c = [Fr(0)]*(deg + 1)
    m = 0
    while 2*m + n <= deg:
        c[2*m + n] = Fr(1, 2**(2*m + n)*factorial(m)*factorial(m + n)); m += 1
    return c

def mul(A, B, deg):
    out = [Fr(0)]*(deg + 1)
    for i, a in enumerate(A):
        if a == 0: continue
        for j in range(0, deg + 1 - i):
            if B[j]: out[i + j] += a*B[j]
    return out

def sign(p):
    s = 1; p = list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s

def det_series(k, deg):
    """det[I_{k+i-j}]_{i,j} as a truncated power series (Leibniz expansion)"""
    I = {n: bessel_I(n, deg) for n in range(k - N, k + N + 1)}
    out = [Fr(0)]*(deg + 1)
    for p in permutations(range(N)):
        term = [Fr(1)] + [Fr(0)]*deg
        for i in range(N):
            term = mul(term, I[k + i - p[i]], deg)
            if not any(term): break
        s = sign(p)
        for t in range(deg + 1): out[t] += s*term[t]
    return out

# the order of det[I_{k+i-j}] is at least N|k| (each entry I_m has order |m|, and sum_i |k + i - p(i)| >= N|k|)
Z = [Fr(0)]*(DEG + 3)
K = 0
while N*K <= DEG + 2:
    for k in ({K, -K} if K else {0}):
        dk = det_series(k, DEG + 2)
        Z = [x + y for x, y in zip(Z, dk)]
    K += 1
Z1 = [(n + 1)*Z[n + 1] for n in range(DEG + 1)]           # Z'
Z2 = [(n + 1)*(n + 2)*Z[n + 2] for n in range(DEG)]      # Z''
Zt = Z[:DEG + 1]
D = DEG - 1   # P is exact to degree DEG - 1 (Z'' known to degree DEG - 1)
P = [2*a - 2*b - c for a, b, c in zip(mul(Z2 + [Fr(0)], Zt, D), mul(Z1, Z1, D), mul(Zt, Zt, D))]
print(f"SU({N}): Z(t) = " + ' + '.join(f'{Z[n]} t^{n}' for n in range(8) if Z[n]) + ' + ...')
print("  Haar moments m_k = k! z_k:", [Z[n]*factorial(n) for n in range(11)])
print("  P(t) = 2Z''Z - 2Z'^2 - Z^2, first coefficients:", [str(P[n]) for n in range(9)])
k0 = next(n for n, p in enumerate(P) if p != 0)
print(f"  first non-zero coefficient: p_{k0} = {P[k0]}")

# compare with the authors' P (rigorous_bounds.upsilon_plus_minus_one, from Schur-Weyl moments) to their order 30
try:
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    from rigorous_bounds import upsilon_plus_minus_one
    Pa = upsilon_plus_minus_one(N, 30)[0]
    same = all(Pa[n] == P[n] for n in range(31))
    print(f"  equal to the authors' P (Schur-Weyl) for every degree <= 30: {same}")
except Exception as err:
    print("  (authors' P not compared:", err, ")")

# rigorous lower bound on (0, t_max], t = beta/N: P(t)/t^k0 >= p_k0 - sum_{k0<n<=D} |p_n| t^(n-k0) - tail(t)/t^k0, every subtracted term increasing in t.
# Tail: the coefficient of t^n in Z, Z', Z'' is at most N^n/n!, N^(n+1)/n!, N^(n+2)/n! (|Re tr U| <= N, moments >= 0), so
#   |p_n| <= sum_{i+j=n} (2 N^(i+2) N^j + 2 N^(i+1) N^(j+1) + N^i N^j)/(i! j!) = (4N^2 + 1) (2N)^n/n!,
# and sum_{n>D} |p_n| t^n <= (4N^2+1) (2Nt)^(D+1)/(D+1)! / (1 - 2Nt/(D+2)).
tm = BMAX/N
lead = P[k0]
rest = sum(abs(P[n])*tm**(n - k0) for n in range(k0 + 1, D + 1))
x = 2*N*tm
tail = Fr(4*N*N + 1)*x**(D + 1)/factorial(D + 1)/(1 - x/(D + 2))
margin = lead - rest - tail/tm**k0
print(f"  on (0, beta_max], beta_max = {BMAX} = {float(BMAX):.6g} (t_max = beta_max/{N}):")
print(f"    p_{k0} = {float(lead):.12g}, sum_(n>{k0}) |p_n| t^(n-{k0}) = {float(rest):.6g}, tail/t^{k0} = {float(tail/tm**k0):.3g}")
print(f"    margin (exact rational) = {float(margin):.12g} > 0: {margin > 0}   =>  U_+ - 1 = P/Z^2 >= 0 on (0, beta_max]")

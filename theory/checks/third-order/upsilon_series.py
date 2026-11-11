# upsilon_series.py -- written by the independent referee of Section 11 of Part IV (2026-11-11), added to the repository with credit;
# unchanged apart from this header. A compass, not part of the proofs.
# Referee check (independent code) for Lemma 11.1 of Part IV, Section 11.
#
# Method (different from the manuscript and from theory/checks/third-order/third_order_values.py):
#   the one-plaquette integral of SU(N) is a sum of Toeplitz determinants of modified Bessel functions,
#       c0(a,b) := int_{SU(N)} exp(a Re chi + b Im chi) dU = sum_{q in Z} e^{i N q theta} D_q(r),   r e^{i theta} = a + i b,
#       D_q(r) := det[ I_{i-j+q}(r) ]_{i,j=1..N}          (Heine/Andreief identity + projection U(N) -> SU(N)).
#   Hence, with F := sum_q D_q and G := sum_q q^2 D_q, at b = 0 and a = x = beta/N:
#       Upsilon_+ = Var(sqrt2 Re chi) = 2 (log F)''(x),
#       Upsilon_- = Var(sqrt2 Im chi) = 2 [F'(x)/x - N^2 G(x)/x^2] / F(x).
#   Everything is done in exact rational arithmetic on truncated power series in x.
#   For U(1): F = I_0(x), G = 0 (no determinant constraint).
#   Also the moment <chi^4 chibar> for SU(3) is computed from <(Re chi)^m (Im chi)^n> read off the 2-variable series.
from fractions import Fraction as Fr
from math import factorial
import itertools, sys

DEG = 12   # truncation degree in x

def mul(p, q, deg=DEG):
    r = [Fr(0)]*(deg+1)
    for i, a in enumerate(p):
        if a == 0: continue
        for j, b in enumerate(q):
            if i+j > deg: break
            if b: r[i+j] += a*b
    return r

def add(p, q):
    return [a+b for a, b in zip(p, q)]

def scal(c, p):
    return [c*a for a in p]

def bessel(n, deg=DEG):
    n = abs(n); p = [Fr(0)]*(deg+1)
    k = 0
    while 2*k+n <= deg:
        p[2*k+n] = Fr(1, 2**(2*k+n)*factorial(k)*factorial(k+n)); k += 1
    return p

def det_series(M, deg=DEG):
    N = len(M); tot = [Fr(0)]*(deg+1)
    for perm in itertools.permutations(range(N)):
        sgn = 1
        for i in range(N):
            for j in range(i+1, N):
                if perm[i] > perm[j]: sgn = -sgn
        t = [Fr(1)] + [Fr(0)]*deg
        for i in range(N):
            t = mul(t, M[i][perm[i]], deg)
            if not any(t): break
        tot = add(tot, scal(sgn, t))
    return tot

def F_G(N, deg=DEG):
    if N == 1:
        return bessel(0, deg), [Fr(0)]*(deg+1)
    F = [Fr(0)]*(deg+1); G = [Fr(0)]*(deg+1)
    qmax = deg//N + 1
    for q in range(-qmax, qmax+1):
        M = [[bessel(i-j+q, deg) for j in range(N)] for i in range(N)]
        D = det_series(M, deg)
        F = add(F, D); G = add(G, scal(q*q, D))
    return F, G

def deriv(p):
    return [i*p[i] for i in range(1, len(p))] + [Fr(0)]

def inv(p, deg):            # 1/p, p[0] != 0
    r = [Fr(0)]*(deg+1); r[0] = 1/p[0]
    for n in range(1, deg+1):
        s = sum(p[k]*r[n-k] for k in range(1, n+1) if k < len(p))
        r[n] = -s/p[0]
    return r

def log_series(p, deg):     # log p, p[0] = 1
    assert p[0] == 1
    dp = deriv(p); q = mul(dp, inv(p, deg), deg)
    out = [Fr(0)]*(deg+1)
    for n in range(1, deg+1): out[n] = q[n-1]/n
    return out

def shift_div(p, k):        # p / x^k (p must vanish to order k)
    assert all(c == 0 for c in p[:k]), p[:k]
    return p[k:] + [Fr(0)]*k

def upsilons(N, deg=DEG):
    F, G = F_G(N, deg)
    assert F[0] == 1
    L = log_series(F, deg)
    Up = scal(2, deriv(deriv(L)))                                  # 2 (log F)''
    Fp = deriv(F)
    num = add(shift_div(Fp, 1), scal(-N*N, shift_div(G, 2)))    # F'/x - N^2 G/x^2
    Um = scal(2, mul(num, inv(F, deg), deg))
    return Up[:deg-3], Um[:deg-3]                                  # reliable degrees

def in_beta(p, N):          # substitute x = beta/N
    return [c/Fr(N)**i for i, c in enumerate(p)]

def mlog(p, deg):
    return scal(-1, log_series(p, deg))

def fmt(p, k=6):
    return ', '.join(str(c) for c in p[:k])

if __name__ == "__main__":
    print("Upsilon_pm as exact power series in beta (coefficients of beta^0..beta^5), and -log Upsilon_pm")
    for N in (1, 2, 3, 4, 5, 6, 7):
        deg = DEG if N <= 6 else 10
        Up, Um = upsilons(N, deg)
        Upb, Umb = in_beta(Up, N), in_beta(Um, N)
        name = "U(1)" if N == 1 else f"SU({N})"
        d = len(Upb)-1
        if N == 2:      # iota = 1: the unit plaquette function is chi itself, Upsilon = Var(chi) = (1/2) Var(sqrt2 Re chi); Im chi = 0
            U = scal(Fr(1, 2), Upb); ml = mlog(U, d)
            print(f"{name} (iota=1): Upsilon = Var(chi) = [{fmt(U)}];  Var(Im chi) series = [{fmt(Umb)}]")
            print(f"   -log Upsilon = [{fmt(ml)}]   beta^3 coeff = {ml[3]}")
            continue
        mlp, mlm = mlog(Upb, d), mlog(Umb, d)
        print(f"{name}:")
        print(f"   Upsilon_+ = [{fmt(Upb)}]")
        print(f"   Upsilon_- = [{fmt(Umb)}]")
        print(f"   -log Upsilon_+ = [{fmt(mlp)}]   beta^3 coeff = {mlp[3]}")
        print(f"   -log Upsilon_- = [{fmt(mlm)}]   beta^3 coeff = {mlm[3]}")
        sys.stdout.flush()
    # Claims of Lemma 11.1 (ii): Upsilon_+ = 1 + 2 r3 t + (r4 - 3/4) t^2 + (r5 - 5 r3) t^3/3, Upsilon_- = 1 - 2 r3 t + (3/4 - r4) t^2 + (m3 + r3) t^3/3
    print("\nLemma 11.1(ii) formulas with the stated r3, r4, r5, m3 (t = beta/N), compared with the Bessel-determinant series:")
    for N in (3, 4, 5, 6, 7):
        r3 = Fr(1, 4)*(N == 3); r4 = Fr(3, 4) + Fr(1, 8)*(N == 4); r5 = Fr(15, 16)*(N == 3) + Fr(1, 16)*(N == 5)
        m3 = -Fr(3, 16)*(N == 3) - Fr(1, 16)*(N == 5)
        t = Fr(1, N)
        claim_p = [Fr(1), 2*r3*t, (r4-Fr(3, 4))*t**2, (r5-5*r3)*t**3/3]
        claim_m = [Fr(1), -2*r3*t, (Fr(3, 4)-r4)*t**2, (m3+r3)*t**3/3]
        deg = DEG if N <= 6 else 10
        Up, Um = upsilons(N, deg)
        Upb, Umb = in_beta(Up, N), in_beta(Um, N)
        ok_p = claim_p == Upb[:4]; ok_m = claim_m == Umb[:4]
        print(f"   SU({N}): Upsilon_+ agrees to beta^3: {ok_p}   Upsilon_- agrees to beta^3: {ok_m}")

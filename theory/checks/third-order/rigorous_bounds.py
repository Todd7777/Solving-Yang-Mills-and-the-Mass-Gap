# Part IV, Section 16 (Lemma mg:lem:upsbounds, Theorem mg:thm:thirdcap): the group data replaced by proven bounds, each beta^k times an increasing function,
# for the computer-assisted version of the third-order term at practical beta.
# Usage: python3 rigorous_bounds.py        (the bounds at beta_0(3), beta_0(4) for SU(3), SU(5), the proof of Upsilon >= 1 for SU(5), and the lower bounds for
#                                           U_+ - U_- of Corollary mg:cor:sectorscap(ii))
"""In "bound mode" the evaluation of Theorem mg:thm:lower5 uses, in place of the exact u, Upsilon, Ubar of the group:
    u      <= ubar   := ((1+iota) beta + beta^2 e^beta)/(2N^2)          (Lemma mg:lem:char(iii)), in every term where u enters increasingly,
    u      >= u_low  := (1+iota) beta e^{-beta}/(2N^2)                   (Lemma mg:lem:char(iii)), where 1/u enters (n_B, c_2, b_2),
    u_2    <= u2bar  := Lemma mg:lem:u2sun (SU(N), N >= 3) or e^beta - 1 - beta (Lemma mg:lem:capelement(iv)),
    sigma6 <= (2-iota) N^2 ubar^6 + u2bar^5 (e^beta - 1)                   (Lemma mg:lem:capelement(iv) with u <= ubar),
    e^{-2 beta} <= Upsilon <= e^{beta}                                     (Lemma mg:lem:upsbounds),
and, from Lemma mg:lem:dressed with these, Upsilon <= Ubar <= Upsilon + Delta, Delta := n_B^2/(e^{-2beta} - c_2) + 2 e^{2beta} sigma6/(e^{-2beta} - ubar^4(1+gamma)).
So Ubar_low := e^{-2 beta} replaces Ubar where it is a lower bound that suffices (the right sides of the conditions on delta', F, F_6; C_varpi), Ubar_up := e^beta +
Delta where Ubar enters increasingly (omega's old condition, eta, 'beyond', 'growth'), and the deficit is bounded by D <= 8(K+1)theta + e^{2beta} Delta + log(1 + delta')
(log(Ubar/Upsilon) <= Delta/Upsilon). Every one of these is beta^k times a function that visibly increases with beta (power series with non-negative coefficients,
their products, and quotients by decreasing positive functions), which is what the scaling of Remark mg:rem:nextorder(vii) needs: no monotonicity is checked on grids."""
import os, sys
from fractions import Fraction
here = os.path.dirname(os.path.abspath(__file__))
for sub in ('second-order', 'upper-bound'): sys.path.insert(0, os.path.join(here, '..', sub))
sys.path.insert(0, here)
from mpmath import mp, mpf, exp, expm1, sqrt, log

def bound_group(N, iota=0, sharp_u2=True):
    """the group data in bound mode (see the module docstring); N >= 3 with sharp_u2 uses Lemma mg:lem:u2sun"""
    from sun_group import u2_bound
    def g(beta):
        beta = mpf(beta)
        ubar = ((1 + iota)*beta + beta**2*exp(beta))/(2*N**2)
        ulow = (1 + iota)*beta*exp(-beta)/(2*N**2)
        u2 = u2_bound(N, beta) if (sharp_u2 and N >= 3) else expm1(beta) - beta
        eps = expm1(2*beta); gam = expm1(2*beta) + beta**2*exp(beta)
        s6 = (2 - iota)*N**2*ubar**6 + u2**5*expm1(beta)
        nB = (u2/ulow)**2*eps + sqrt(2 - iota)*(N/ulow)*sqrt(u2**5*expm1(beta))
        c2 = (u2/ulow)**4*exp(beta)
        ups_low, ups_up = exp(-2*beta), exp(beta)
        assert ups_low > c2 and ups_low > ubar**4*(1 + gam), 'the denominators of Lemma mg:lem:dressed'
        Delta = nB**2/(ups_low - c2) + 2*exp(2*beta)*s6/(ups_low - ubar**4*(1 + gam))
        return dict(name=f'SU({N}) bounds', N=N, iota=iota, bounds=True, u=ubar, u_low=ulow, u2=u2, s6=s6,
                    Ubar_low=ups_low, Ubar_up=ups_up + Delta, Delta=Delta, Ups_low=ups_low)
    return g

# --- Upsilon >= 1 for SU(N) (used by the upper bound, Remark mg:rem:upperexact(ii)): U_+ - 1 = P(t)/Z(t)^2 with P = 2 Z_2 Z - 2 Z_1^2 - Z^2, t = beta/N,
# Z = G(0,0), Z_1 = (G(1,0)+G(0,1))/2, Z_2 = (G(2,0)+2G(1,1)+G(0,2))/4 in the notation of sun_group (all with non-negative coefficients). P is computed
# exactly to degree K; the rest is bounded with z_k <= N^k/k! (for Z), N^(k+1)/k! (Z_1), N^(k+2)/k! (Z_2) (as m(a,b) <= N^(a+b) and sum_{a+b=k} 1/(2^k a! b!) = 1/k!).
def upsilon_plus_minus_one(N, K=30):
    """exact coefficients p_k (k <= K) of P(t) = U_+ Z^2 - Z^2, and the coefficient majorants needed for the tail"""
    from sun_group import coeffs
    C = coeffs(N, K)
    Z = C[(0, 0)]; Z1 = [(x + y)/2 for x, y in zip(C[(1, 0)], C[(0, 1)])]; Z2 = [(x + 2*y + z)/4 for x, y, z in zip(C[(2, 0)], C[(1, 1)], C[(0, 2)])]
    def mul(A, B): return [sum(A[i]*B[k - i] for i in range(k + 1)) for k in range(K + 1)]
    ZZ, Z1Z1, Z2Z = mul(Z, Z), mul(Z1, Z1), mul(Z2, Z)
    P = [2*a - 2*b - c for a, b, c in zip(Z2Z, Z1Z1, ZZ)]
    return P, Z, Z1, Z2

def prove_upsilon_ge_1(N, beta_max, K=30):
    """a proof that U_+ >= 1 on (0, beta_max] for SU(N) (N >= 4, where the first coefficients of P vanish): P(t) = sum_k p_k t^k with p_k exact for
    k <= K, and |sum_{k>K} p_k t^k| <= tail(t), a bound by majorants; with k0 the first k with p_k != 0 (it must be > 0), U_+ - 1 >= 0 on (0, t_max] follows from
    p_k0 - sum_{k0<k<=K} |p_k| t_max^(k-k0) - tail(t_max)/t_max^k0 > 0 (every subtracted term increases with t)."""
    P, Z, Z1, Z2 = upsilon_plus_minus_one(N, K)
    k0 = next(k for k, p in enumerate(P) if p != 0)
    assert P[k0] > 0, f'the first coefficient of U_+ - 1 is negative ({P[k0]})'
    tm = mpf(beta_max)/N*(1 + mpf(10)**-20)       # t_max, rounded up
    # tail of each series beyond K, with z_k <= c N^k/k!
    from math import factorial
    def tail(cN, tt):           # sum_{k>K} cN N^k tt^k/k!
        return cN*sum((N*tt)**k/factorial(k) for k in range(K + 1, K + 60))*(1 + mpf(10)**-20)
    def full(series, cN, tt):   # the series at tt, majorised: exact part + tail
        return sum(mpf(c.numerator)/c.denominator*tt**k for k, c in enumerate(series)) + tail(cN, tt)
    # the degree > K part of a product AB of series with non-negative coefficients, written without cancellation:
    # [A_K B_K]_{deg>K} + A_K tail(B) + tail(A) B_K + tail(A) tail(B), with A_K, B_K the truncations (exact) and tail() the majorants above
    def trunc(A, tt): return sum(mpf(c.numerator)/c.denominator*tt**k for k, c in enumerate(A))
    def prod_tail(A, cA, B, cB, tt):
        hi = sum(mpf(sum(A[i]*B[k - i] for i in range(k - K, K + 1)).numerator)/sum(A[i]*B[k - i] for i in range(k - K, K + 1)).denominator*tt**k
                 for k in range(K + 1, 2*K + 1))
        tA, tB = tail(cA, tt), tail(cB, tt)
        return hi + trunc(A, tt)*tB + tA*trunc(B, tt) + tA*tB
    T = 2*prod_tail(Z2, N**2, Z, 1, tm) + 2*prod_tail(Z1, N, Z1, N, tm) + prod_tail(Z, 1, Z, 1, tm)
    lead = mpf(P[k0].numerator)/P[k0].denominator
    rest = sum(abs(mpf(p.numerator)/p.denominator)*tm**(k - k0) for k, p in enumerate(P) if k > k0)
    assert T >= 0
    margin = lead - rest - T/tm**k0
    return dict(k0=k0, lead=P[k0], t_max=tm, rest=rest, tail=T/tm**k0, margin=margin, ok=margin > 0)

def split_lower(N, beta_max, K=30):
    """a lower bound for U_+ - U_- (Corollary mg:cor:sectorscap(ii)). With t = beta/N, Z, Z_1, Z_2 as above and Z_im := int (Im chi)^2 e^{t Re chi} dU
    (= -(G(2,0) - 2G(1,1) + G(0,2))/4 with G(i,j) := int chi^i chibar^j e^{t Re chi} dU): U_+ = 2(Z_2 Z - Z_1^2)/Z^2 and U_- = 2 Z_im Z/Z^2, so
    U_+ - U_- = 2Q/Z^2 with Q := Z_2 Z - Z_1^2 - Z_im Z. The coefficients of Q are exact to degree K; the coefficients of Z_im may have either sign, so every
    product's part of degree > K is bounded with absolute values: |coefficient of t^k| <= N^k/k!, N^{k+1}/k!, N^{k+2}/k!, N^{k+2}/k! for Z, Z_1, Z_2, Z_im
    (|Re chi|, |Im chi| <= N). With k0 the first k with q_k != 0 (it must be > 0): Q(t) >= t^k0 margin on (0, t_max], margin := q_k0 - sum_{k0<k<=K} |q_k|
    t_max^(k-k0) - tail(t_max)/t_max^k0. Then, as Z <= e^beta and U_+ <= e^beta (Lemma mg:lem:upsbounds(i)), log(U_+/U_-) >= (U_+ - U_-)/U_+ >= 2 e^{-3 beta} Q."""
    from sun_group import coeffs
    from math import factorial
    C = coeffs(N, K)
    Z = C[(0, 0)]
    Z1 = [(x + y)/2 for x, y in zip(C[(1, 0)], C[(0, 1)])]
    Z2 = [(x + 2*y + z)/4 for x, y, z in zip(C[(2, 0)], C[(1, 1)], C[(0, 2)])]
    Zi = [-(x - 2*y + z)/4 for x, y, z in zip(C[(2, 0)], C[(1, 1)], C[(0, 2)])]
    def mul(A, B): return [sum(A[i]*B[k - i] for i in range(k + 1)) for k in range(K + 1)]
    Q = [a - b - c for a, b, c in zip(mul(Z2, Z), mul(Z1, Z1), mul(Zi, Z))]
    k0 = next(k for k, q in enumerate(Q) if q != 0)
    assert Q[k0] > 0, f'the first coefficient of U_+ - U_- is negative ({Q[k0]})'
    tm = mpf(beta_max)/N*(1 + mpf(10)**-20)
    def fr(c): return mpf(c.numerator)/c.denominator
    def tail(cN, tt): return mpf(cN)*sum((N*tt)**k/factorial(k) for k in range(K + 1, K + 80))*(1 + mpf(10)**-20)
    def absval(A, tt): return sum(abs(fr(c))*tt**k for k, c in enumerate(A))
    def prod_tail(A, cA, B, cB, tt):     # |degree > K part of AB| <= sum_k |[A_K B_K]_k| tt^k + |A_K| tail(B) + tail(A) |B_K| + tail(A) tail(B)
        hi = sum(abs(fr(sum(A[i]*B[k - i] for i in range(k - K, K + 1))))*tt**k for k in range(K + 1, 2*K + 1))
        tA, tB = tail(cA, tt), tail(cB, tt)
        return hi + absval(A, tt)*tB + tA*absval(B, tt) + tA*tB
    T = prod_tail(Z2, N**2, Z, 1, tm) + prod_tail(Z1, N, Z1, N, tm) + prod_tail(Zi, N**2, Z, 1, tm)
    lead = fr(Q[k0])
    rest = sum(abs(fr(q))*tm**(k - k0) for k, q in enumerate(Q) if k > k0)
    margin = lead - rest - T/tm**k0
    return dict(k0=k0, lead=Q[k0], t_max=tm, rest=rest, tail=T/tm**k0, margin=margin, ok=margin > 0)

def sector_series(N, K=6):
    """the exact Taylor coefficients in beta, to order K, of -log U_+, -log U_- and log(U_+/U_-) (Fractions; U_+ = 2(Z_2 Z - Z_1^2)/Z^2, U_- = 2 Z_im Z/Z^2 as in
    split_lower, with t = beta/N): the series of Corollary mg:cor:sectorscap(i), (ii)."""
    from sun_group import coeffs
    C = coeffs(N, K + 2)
    Z = C[(0, 0)][:K + 1]
    Z1 = [(x + y)/2 for x, y in zip(C[(1, 0)], C[(0, 1)])][:K + 1]
    Z2 = [(x + 2*y + z)/4 for x, y, z in zip(C[(2, 0)], C[(1, 1)], C[(0, 2)])][:K + 1]
    Zi = [-(x - 2*y + z)/4 for x, y, z in zip(C[(2, 0)], C[(1, 1)], C[(0, 2)])][:K + 1]
    def mul(A, B): return [sum(A[i]*B[k - i] for i in range(k + 1)) for k in range(K + 1)]
    def inv(A):
        out = [Fraction(1)/A[0]]
        for k in range(1, K + 1): out.append(-sum(A[i]*out[k - i] for i in range(1, k + 1))/A[0])
        return out
    def neglog(A):                    # -log A for A(0) = 1: -sum_{j>=1} (-1)^{j+1} (A-1)^j / j
        assert A[0] == 1, A[0]
        x = [Fraction(0)] + A[1:]; out = [Fraction(0)]*(K + 1); pw = [Fraction(1)] + [Fraction(0)]*K
        for j in range(1, K + 1):
            pw = mul(pw, x)
            out = [o - Fraction((-1)**(j + 1), j)*q for o, q in zip(out, pw)]
        return out
    iZ2 = mul(inv(Z), inv(Z))
    Up = [2*c for c in mul([a - b for a, b in zip(mul(Z2, Z), mul(Z1, Z1))], iZ2)]
    Um = [2*c for c in mul(mul(Zi, Z), iZ2)]
    lp, lm = neglog(Up), neglog(Um)
    tob = lambda L: [c/Fraction(N)**k for k, c in enumerate(L)]          # coefficients of t^k -> of beta^k
    return dict(neglog_plus=tob(lp), neglog_minus=tob(lm), log_ratio=tob([b - a for a, b in zip(lp, lm)]))

def classic_series(K=8):
    """the exact Taylor coefficients in beta, to order K, of -log Upsilon for SU(2) (iota = 1: Upsilon = Var(chi) under e^{(beta/2) chi}, with the Haar moments
    int chi^{2k} = Catalan(k)), and of -log U_+, -log U_-, log(U_-/U_+) for U(1) (U_+ = Var(sqrt2 cos), U_- = Var(sqrt2 sin) under e^{beta cos}, with
    int cos^{2k} = binom(2k,k)/4^k): Corollary mg:cor:classic."""
    from math import comb, factorial
    F = Fraction
    def mul(A, B): return [sum(A[i]*B[k - i] for i in range(k + 1)) for k in range(K + 1)]
    def inv(A):
        out = [F(1)/A[0]]
        for k in range(1, K + 1): out.append(-sum(A[i]*out[k - i] for i in range(1, k + 1))/A[0])
        return out
    def neglog(A):
        x = [F(0)] + A[1:]; out = [F(0)]*(K + 1); pw = [F(1)] + [F(0)]*K
        for j in range(1, K + 1):
            pw = mul(pw, x); out = [o - F((-1)**(j + 1), j)*q for o, q in zip(out, pw)]
        return out
    cat = lambda k: F(comb(2*k, k), k + 1)
    # SU(2): moments m_j = int chi^j (0 for j odd); Z_i = sum_j (beta/2)^j m_{j+i}/j!
    m = lambda j: cat(j//2) if j % 2 == 0 else F(0)
    Zs = [[F(1, 2**j)*m(j + i)/factorial(j) for j in range(K + 1)] for i in range(3)]
    iZ = inv(Zs[0]); Ups = [a - b for a, b in zip(mul(Zs[2], iZ), mul(mul(Zs[1], iZ), mul(Zs[1], iZ)))]
    su2 = neglog(Ups)
    # U(1): moments c_j = int cos^j (0 for j odd), Z_i = sum_j beta^j c_{j+i}/j!;  E cos^2 = Z_2/Z, E cos = Z_1/Z
    c = lambda j: F(comb(j, j//2), 2**j) if j % 2 == 0 else F(0)
    W = [[c(j + i)/factorial(j) for j in range(K + 1)] for i in range(3)]
    iW = inv(W[0]); Ec = mul(W[1], iW); Ec2 = mul(W[2], iW)
    Up = [2*(a - b) for a, b in zip(Ec2, mul(Ec, Ec))]
    Um = [2*(F(1) if k == 0 else F(0)) - 2*a for k, a in enumerate(Ec2)]
    lp, lm = neglog(Up), neglog(Um)
    return dict(su2=su2, u1_plus=lp, u1_minus=lm, u1_log_ratio=[a - b for a, b in zip(lp, lm)])       # log(U_-/U_+) = (-log U_+) - (-log U_-)

def u1_split_lower(beta_max, K=30):
    """U(1): U_- - U_+ = 2 (I_1^2 - I_0 I_2)/I_0^2. The series P := I_1^2 - I_0 I_2 = sum_k p_k (beta/2)^{2k+2} is computed exactly to k = K; the rest is at most
    the part of degree > K of the products of the majorant series sum_k (beta/2)^{2k+nu}/(k!(k+nu)!) (non-negative), so P >= (beta/2)^2 (p_0 - sum_{0<k<=K} |p_k| x^k
    - tail/(beta/2)^2) with x = (beta_max/2)^2. Then log(U_-/U_+) >= (U_- - U_+)/U_- >= 2 e^{-3beta} P (I_0^2 <= e^{2beta}, U_- <= e^beta)."""
    from math import factorial
    F = Fraction
    def I(nu): return [F(1, factorial(k)*factorial(k + nu)) for k in range(K + 1)]          # coefficients of (beta/2)^{2k+nu}
    def mul(A, B): return [sum(A[i]*B[k - i] for i in range(k + 1)) for k in range(K + 1)]
    P = [a - b for a, b in zip(mul(I(1), I(1)), mul(I(0), I(2)))]                             # coefficients of (beta/2)^{2k+2}
    assert P[0] == F(1, 2), P[0]
    x = (mpf(beta_max)/2)**2*(1 + mpf(10)**-20)
    def full(nu, kmax=K + 60): return [mpf(1)/(factorial(k)*factorial(k + nu)) for k in range(kmax)]
    def tailprod(nu1, nu2):        # the part of degree > K (in x) of the product of the two non-negative series, at x
        A, B = full(nu1), full(nu2)
        return sum(A[i]*B[j]*x**(i + j) for i in range(len(A)) for j in range(len(B)) if i + j > K)
    T = tailprod(1, 1) + tailprod(0, 2)
    rest = sum(abs(mpf(p.numerator)/p.denominator)*x**k for k, p in enumerate(P) if k > 0)
    margin = mpf(1)/2 - rest - T
    return dict(lead=P[0], rest=rest, tail=T, margin=margin, ok=margin > 0)

if __name__ == "__main__":
    from constants_upper import beta0
    mp.dps = 60
    for N in (3, 5):
        for d in (3, 4):
            b0 = beta0(d); g = bound_group(N)(b0)
            print(f"SU({N}), d={d}, beta_0 = {float(b0):.6g}: ubar/beta = {float(g['u']/b0):.8f}, u_low/beta = {float(g['u_low']/b0):.8f}, u2bar/beta^2 = {float(g['u2']/b0**2):.6g},"
                  f" Delta/beta^4 = {float(g['Delta']/b0**4):.4g}, Ubar in [{float(g['Ubar_low']):.8f}, {float(g['Ubar_up']):.8f}]")
    mp.dps = 120
    r = prove_upsilon_ge_1(5, mpf('4.03e-4'))         # the range of Lemma mg:lem:upsbounds(ii), which contains (0, beta_0(3)] (the ninth referee)
    print(f"SU(5): U_+ - 1 = P(t)/Z^2 with P = {r['lead']} t^{r['k0']} + ...; on (0, 4.03e-4/5]: {float(r['lead']):.6g} - {float(r['rest']):.3g} - {float(r['tail']):.3g}"
          f" = {float(r['margin']):.6g} > 0: {r['ok']}  (so Upsilon >= U_+ >= 1 on (0, 4.03e-4], which contains (0, beta_0(3)] and (0, beta_0(4)])")
    for N in (3, 4, 5):               # the exact series of Corollary mg:cor:sectorscap(i), (ii)
        S = sector_series(N)
        for key, lab in (('neglog_plus', '-log U_+'), ('neglog_minus', '-log U_-'), ('log_ratio', 'log(U_+/U_-)')):
            print(f"SU({N}): {lab} = " + ' + '.join(f"({c}) beta^{k}" for k, c in enumerate(S[key]) if c != 0) + " + O(beta^7)")
    S = classic_series()              # Corollary mg:cor:classic: SU(2) and U(1)
    for key, lab in (('su2', 'SU(2): -log Upsilon'), ('u1_plus', 'U(1): -log U_+'), ('u1_minus', 'U(1): -log U_-'), ('u1_log_ratio', 'U(1): log(U_-/U_+)')):
        print(f"{lab} = " + ' + '.join(f"({c}) beta^{k}" for k, c in enumerate(S[key]) if c != 0) + " + O(beta^9)")
    r = u1_split_lower(mpf('4.03e-4'))
    print(f"U(1): U_- - U_+ = 2(I_1^2 - I_0 I_2)/I_0^2, I_1^2 - I_0 I_2 = (beta/2)^2 ({r['lead']} + ...); on (0, 4.03e-4]: >= (beta/2)^2 ({float(r['lead']):.6g} - {float(r['rest']):.3g}"
          f" - {float(r['tail']):.3g}) = (beta/2)^2 {float(r['margin']):.12g} > 0: {r['ok']};  so log(U_-/U_+) >= 2 e^(-3beta) (I_1^2 - I_0 I_2) >= e^(-3beta) beta^2 {float(r['margin']/2):.12g}")
    for N in (3, 4, 5):               # the splitting of the two plaquette glueballs (Corollary mg:cor:sectorscap(ii))
        r = split_lower(N, beta0(3)); k0 = r['k0']
        print(f"SU({N}): U_+ - U_- = 2Q/Z^2, Q = {r['lead']} t^{k0} + ...; on (0, beta_0(3)/{N}]: Q >= t^{k0} ({float(r['lead']):.6g} - {float(r['rest']):.3g}"
              f" - {float(r['tail']):.3g}) = t^{k0} {float(r['margin']):.9g} > 0: {r['ok']};  so log(U_+/U_-) >= 2 e^(-3 beta) Q >= e^(-3 beta) beta^{k0}"
              f" {float(2*r['margin']/mpf(N)**k0):.9g}  (1/{float(mpf(N)**k0/(2*r['margin'])):.6f})")

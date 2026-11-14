# Referee check (eleventh pass on Section 16, Corollary mg:cor:classic), item 4: the series, recomputed independently of rigorous_bounds.classic_series.
# The authors expand Haar moments (Catalan numbers for SU(2), binom(2k,k)/4^k for U(1)). Here I use closed forms in modified Bessel functions instead:
#   SU(2): w_beta = e^{(beta/2) chi}, Z(beta) = int e^{(beta/2) chi} dU = 2 I_1(beta)/beta (Weyl: dU = (2/pi) sin^2 t dt, chi = 2 cos t), and as psi = chi_rho is
#          the unit eigenvector (iota = 1, int chi^2 dU = 1), Upsilon = Var(chi) = 4 (log Z)''.
#   U(1):  w_beta = e^{beta cos}, Z = I_0;  Upsilon_+ = Var(sqrt2 cos) = 2 (log I_0)'' = 1 + I_2/I_0 - 2 I_1^2/I_0^2;
#          Upsilon_- = Var(sqrt2 sin) = 2 <sin^2> = 1 - I_2/I_0 = 2 I_1/(beta I_0)   (I_0 - I_2 = (2/beta) I_1).
# (a) exact Taylor coefficients (Fractions) from the Bessel series, by power-series log/derivatives (my own code);
# (b) high-precision numerics (mpmath besseli, 80 digits) at several beta, against the series truncated after beta^8 (the remainder must be O(beta^10));
# (c) the Turan determinant P := I_1^2 - I_0 I_2 = sum_k p_k (beta/2)^{2k+2}: exact p_k to k = 60 (are they all positive?), and a rigorous lower bound on
#     (0, 4.03e-4] in Arb with my own tail bound |p_k| <= 2 (2^k/k!)^2 (from 1/(i!(i+nu)!) <= 1/(i!)^2 and sum_{i+j=k} 1/(i!j!) = 2^k/k!);
# (d) the chain log(U_-/U_+) >= (U_- - U_+)/U_- >= 2 e^{-3beta} P >= e^{-3beta}(1 - 1e-7) beta^2/4, checked numerically at many beta (sanity).
# Usage: python3 series_check_classic.py
import math, os, sys
from fractions import Fraction as F
from mpmath import mp, mpf, besseli, log, exp

K = 14                                        # internal order in beta (derivatives and the shift I_1/beta lose the top orders; results are used to beta^10)
KV = 10                                       # the orders that are valid and printed

def mul(A, B): return [sum(A[i]*B[k - i] for i in range(k + 1)) for k in range(K + 1)]
def inv(A):
    out = [F(1)/A[0]]
    for k in range(1, K + 1): out.append(-sum(A[i]*out[k - i] for i in range(1, k + 1))/A[0])
    return out
def plog(A):                                  # log A for A(0) = 1, via log A = int A'/A
    assert A[0] == 1
    dA = [(k + 1)*A[k + 1] for k in range(K)] + [F(0)]
    q = mul(dA, inv(A))
    return [F(0)] + [q[k - 1]/k for k in range(1, K + 1)]
def deriv(A): return [(k + 1)*A[k + 1] for k in range(K)] + [F(0)]
def bessel_I(nu):                             # I_nu(beta) = sum_k (beta/2)^{2k+nu}/(k!(k+nu)!), as a series in beta (exact)
    out = [F(0)]*(K + 3)
    for k in range(K + 2):
        e = 2*k + nu
        if e <= K + 2: out[e] += F(1, 2**e*math.factorial(k)*math.factorial(k + nu))
    return out[:K + 1]

def series():
    # SU(2): Z = 2 I_1/beta = sum_k (beta/2)^{2k}/(k!(k+1)!)
    Z = [F(0)]*(K + 1)
    for k in range(K//2 + 1): Z[2*k] = F(1, 4**k*math.factorial(k)*math.factorial(k + 1))
    Ups = [4*c for c in deriv(deriv(plog(Z)))]
    su2 = [-c for c in plog(Ups)]
    I0, I1, I2 = bessel_I(0), bessel_I(1), bessel_I(2)
    iI0 = inv(I0)
    Up = [2*c for c in deriv(deriv(plog(I0)))]
    # U_- = 2 I_1/(beta I_0): I_1/beta as a series
    I1b = [I1[k + 1] for k in range(K)] + [F(0)]
    Um = [2*c for c in mul(I1b, iI0)]
    # cross-check U_- = 1 - I_2/I_0 and U_+ = 1 + I_2/I_0 - 2 I_1^2/I_0^2 exactly
    Um2 = [(1 if k == 0 else 0) - c for k, c in enumerate(mul(I2, iI0))]
    Up2 = [(1 if k == 0 else 0) + a - 2*b for k, (a, b) in enumerate(zip(mul(I2, iI0), mul(mul(I1, I1), mul(iI0, iI0))))]
    assert Um[:KV + 1] == Um2[:KV + 1] and Up[:KV + 1] == Up2[:KV + 1], 'the two closed forms of U_+-'
    lp, lm = [-c for c in plog(Up)], [-c for c in plog(Um)]
    cut = lambda L: L[:KV + 1]
    return dict(su2=cut(su2), u1_plus=cut(lp), u1_minus=cut(lm), u1_log_ratio=cut([a - b for a, b in zip(lp, lm)]), Ups_su2=cut(Ups), Up=cut(Up), Um=cut(Um))

def fmt(L): return ' + '.join(f"({c}) beta^{k}" for k, c in enumerate(L) if c != 0)

if __name__ == "__main__":
    S = series()
    print("(a) exact series from the Bessel closed forms (to beta^10):")
    for key, lab in (('su2', 'SU(2): -log Upsilon'), ('u1_plus', 'U(1): -log U_+'), ('u1_minus', 'U(1): -log U_-'), ('u1_log_ratio', 'U(1): log(U_-/U_+)')):
        print(f"  {lab} = {fmt(S[key])} + O(beta^12)")
    claims = {'su2': (F(1, 8), F(-1, 192)), 'u1_plus': (F(3, 8), F(-13, 384)), 'u1_minus': (F(1, 8), F(-5, 384)), 'u1_log_ratio': (F(1, 4), F(-1, 48))}
    for key, (c2, c4) in claims.items():
        L = S[key]
        print(f"  {key}: odd coefficients all zero: {all(L[k] == 0 for k in range(1, KV + 1, 2))};  beta^2 = {L[2]} (claimed {c2}): {L[2] == c2};"
              f"  beta^4 = {L[4]} (claimed {c4}): {L[4] == c4}")
    # against the authors' classic_series
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    import rigorous_bounds as RB
    A = RB.classic_series(8)
    print("  identical to rigorous_bounds.classic_series to beta^8:", all(A[k][:9] == S[k][:9] for k in ('su2', 'u1_plus', 'u1_minus', 'u1_log_ratio')))
    # (b) numerics
    mp.dps = 80
    print("\n(b) high-precision numerics (mpmath besseli, 80 digits): (exact - series to beta^8)/beta^10 should stay bounded:")
    def num(beta):
        beta = mpf(beta)
        I0, I1, I2, I3 = (besseli(n, beta) for n in range(4))
        # SU(2): Upsilon = Var(chi) under e^{(beta/2)chi}: <chi> = 2 I_2/I_1, <chi^2> = 4 (I_1'' ... ) -> use Z = 2I_1/beta and numerical second derivative of log Z
        Z = lambda b: 2*besseli(1, b)/b
        Ups = 4*mp.diff(lambda b: log(Z(b)), beta, 2)
        Up = 1 + I2/I0 - 2*I1**2/I0**2
        Um = 2*I1/(beta*I0)
        return dict(su2=-log(Ups), u1_plus=-log(Up), u1_minus=-log(Um), u1_log_ratio=log(Um/Up), u_su2=I2/I1, u_u1=I1/I0)
    for beta in ('1e-1', '3e-2', '1e-2', '4.0262143352282e-4', '1e-4'):
        n = num(beta); b = mpf(beta)
        row = []
        for key in ('su2', 'u1_plus', 'u1_minus', 'u1_log_ratio'):
            ser = sum(mpf(c.numerator)/c.denominator*b**k for k, c in enumerate(S[key][:9]))
            row.append(f"{key} {mp.nstr((n[key] - ser)/b**10, 6)}")
        print(f"  beta = {beta}: " + ', '.join(row))
    # (c) the Turan determinant
    print("\n(c) P = I_1^2 - I_0 I_2 = sum_k p_k (beta/2)^{2k+2}:")
    KP = 60
    def I_(nu): return [F(1, math.factorial(k)*math.factorial(k + nu)) for k in range(KP + 1)]
    def mulP(A, B): return [sum(A[i]*B[k - i] for i in range(k + 1)) for k in range(KP + 1)]
    P = [a - b for a, b in zip(mulP(I_(1), I_(1)), mulP(I_(0), I_(2)))]
    print(f"  p_0..p_5 = {[str(p) for p in P[:6]]};  all p_k > 0 for k <= {KP}: {all(p > 0 for p in P)}")
    # closed form from the product formula I_mu I_nu = sum_k (x/2)^{mu+nu+2k} (mu+nu+2k)!/(k! (mu+nu+k)! (mu+k)! (nu+k)!) (the I-analogue of DLMF 10.8.3):
    # p_k = (2k+2)!/(k!(k+2)!) [1/((k+1)!)^2 - 1/(k!(k+2)!)] = (2k+2)!/(k! (k+2)! ((k+1)!)^2 (k+2)) > 0, so P >= (beta/2)^2/2 for EVERY beta > 0
    fac = math.factorial
    closed = [F(fac(2*k + 2), fac(k)*fac(k + 2)*fac(k + 1)**2*(k + 2)) for k in range(KP + 1)]
    print(f"  p_k = (2k+2)!/(k!(k+2)!((k+1)!)^2(k+2)) for all k <= {KP}: {closed == P}  (so I_1^2 - I_0 I_2 >= (beta/2)^2/2 for all beta > 0, by the product formula)")
    from flint import arb, ctx
    ctx.prec = 300
    bmax = arb('4.03e-4'); x = (bmax/2)**2
    rest = sum((arb(p.numerator)/arb(p.denominator)*x**k for k, p in enumerate(P) if 1 <= k <= 40), arb(0))
    # sum_{k>40} |p_k| x^k <= sum_{k>40} 2 (2^k/k!)^2 x^k <= 2 (2^41/41!)^2 x^41 / (1 - 4x/42^2): the ratio of consecutive terms is 4x/(k+1)^2 <= 4x/42^2
    r = 4*x/arb(42)**2; assert bool(1 - r > 0)
    tail = 2*(arb(2)**41/arb(math.factorial(41)))**2*x**41/(1 - r)
    # with all p_k > 0 the lower bound needs no subtraction at all: P >= p_0 (beta/2)^2 = (beta/2)^2/2. The authors' route subtracts |rest| + tail:
    margin = arb(1)/2 - rest - tail
    print(f"  on (0, 4.03e-4], authors' route: P >= (beta/2)^2 (1/2 - {rest.str(4, radius=False)} - {tail.str(3, radius=False)}) = (beta/2)^2 * {margin.str(15, radius=False)}"
          f" = (beta/2)^2 (1/2)(1 - {(1 - 2*margin).str(4, radius=False)});  >= (1/2)(1 - 3e-8): {bool(margin - (1 - arb('3e-8'))/2 > 0)}")
    print("  (since p_0 > 0 and p_k > 0 for the computed k, P >= (beta/2)^2/2 on (0, inf) would follow from positivity of all p_k, which I only checked to k = 60)")
    # (d) the chain, numerically (sanity), on a grid of (0, beta_0(3)]
    print("\n(d) the chain log(U_-/U_+) >= (U_- - U_+)/U_- >= 2 e^{-3 beta} P >= e^{-3beta}(1 - 1e-7) beta^2/4, at 2000 points of (0, 4.03e-4] (Numerical):")
    ok = True; worst = None
    for j in range(1, 2001):
        b = mpf('4.03e-4')*j/2000
        I0, I1, I2 = (besseli(n, b) for n in range(3))
        Up = 1 + I2/I0 - 2*I1**2/I0**2; Um = 2*I1/(b*I0); Pn = I1**2 - I0*I2
        c1 = log(Um/Up); c2 = (Um - Up)/Um; c3 = 2*exp(-3*b)*Pn; c4 = exp(-3*b)*(1 - mpf(10)**-7)*b**2/4
        good = c1 >= c2 >= c3 >= c4 and abs(Um - Up - 2*Pn/I0**2) < mpf(10)**-70
        ok = ok and good
        r = c1/(b**2/4)
        worst = r if worst is None or r < worst else worst
    print(f"  all hold: {ok};  min over the grid of log(U_-/U_+)/(beta^2/4) = {mp.nstr(worst, 12)}")

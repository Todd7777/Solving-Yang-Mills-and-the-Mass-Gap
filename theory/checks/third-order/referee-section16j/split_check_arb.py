# Referee check (tenth pass on Section 16, Corollary mg:cor:sectorscap), item 4: the lower bound for the splitting, independently of rigorous_bounds.split_lower.
#  - Q = Z_2 Z - Z_1^2 - Z_im Z (t = beta/N) with EXACT coefficients q_k, k <= KQ, from the referee's Toeplitz-determinant moments (series_check.py);
#  - the rest by the referee's own bound: |coefficient of t^k in each of Z_2 Z, Z_1^2, Z_im Z| <= N^(k+2) 2^k/k! (as |Re chi|, |Im chi| <= N, so the t^i
#    coefficients of Z, Z_1, Z_2, Z_im are at most N^i/i!, N^(i+1)/i!, N^(i+2)/i!, N^(i+2)/i! in absolute value), so |q_k| <= 3 N^2 (2N)^k/k! for every k;
#  - margin(t) := q_k0 - sum_{k0<k<=KQ} |q_k| t^(k-k0) - 3N^2 sum_{k>KQ} (2N)^k t^(k-k0)/k!, every subtracted term increasing in t, so Q(t) >= t^k0 margin(t_max)
#    on (0, t_max], t_max = beta_0(3)/N (an upper bound in Arb);
#  - then log(U_+/U_-) >= (U_+ - U_-)/U_+ = 2Q/(Z^2 U_+) >= 2 e^{-3 beta} Q (Z <= e^beta, U_+ <= e^beta), and the claim of (ii):
#    m_- - m_+ >= e^{-3beta}(1 - 1e-8) 2(beta/2N)^(N-2)/(N-2)! - (c'_- + c'_+) beta^4 > 0 for 0 < beta < beta_0(d), checked at beta_0(d) (the bracket after
#    dividing by beta^(N-2) decreases), with the table's constants, in Arb;
#  - also the numbers quoted in the proof: 3.329e-4 vs 1.58e-4, 3.332e-4 vs 1.36e-4, the factors 2.1 and 2.4, 0.642, 3.984 vs 3.83, 10.04, and 1.91, 5.01 of (i).
# Usage: python3 split_check_arb.py
import os, sys
from fractions import Fraction as Fr
from math import factorial
from flint import arb, ctx
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import series_check as SC
ctx.prec = 300
ONE = arb(1); E = ONE.exp()
def beta0(d):
    K = 4*(2*d - 3); return (ONE + ONE/(E**2*K*(K + 2))).log()/2
def A(fr): return arb(fr.numerator)/arb(fr.denominator)
CONST = {(3, 3): ('0.600', '0.042'), (3, 4): ('3.937', '0.047'), (4, 3): ('0.134', '0.131'), (4, 4): ('1.083', '0.148'), (5, 3): ('0.073', '0.319'), (5, 4): ('0.520', '0.361')}

if __name__ == "__main__":
    DEG = 16                                       # moments to total degree 16, q_k exact for k <= 14
    b03, b04 = beta0(3), beta0(4)
    print(f"beta_0(3) = {b03.str(12)}, beta_0(4) = {b04.str(12)}")
    allok = True
    for N in (3, 4, 5):
        s = SC.sector(N, DEG); Q = s['Q']; KQ = s['K']
        k0 = next(k for k, q in enumerate(Q) if q); qN = Q[k0]
        tmax = arb(b03.upper())/N; tmax = arb(tmax.upper())
        rest = sum((abs(A(Q[k]))*tmax**(k - k0) for k in range(k0 + 1, KQ + 1)), arb(0))
        z = 2*N*tmax; term = 3*N**2*(2*N)**(KQ + 1)*tmax**(KQ + 1 - k0)/factorial(KQ + 1)
        tail = term/(1 - z/(KQ + 2))
        margin = A(qN) - rest - tail
        ok1 = bool(margin - A(qN)*(1 - arb('1e-8')) > 0)
        allok = allok and ok1
        print(f"SU({N}): Q = ({qN}) t^{k0} + ..., q_N = 2^-(N-2)/(N-2)!: {qN == Fr(1, 2**(N - 2)*factorial(N - 2))};  on (0, beta_0(3)/N]:"
              f" margin = {margin.str(12)} = q_N (1 - {((A(qN) - margin)/A(qN)).str(3)}) >= q_N (1 - 1e-8): {ok1}")
        for d, b0 in ((3, b03), (4, b04)):
            cm, cp = (arb(x) for x in CONST[(N, d)])
            lead = 2*(b0/(2*N))**(N - 2)/factorial(N - 2)
            lhs = (-3*b0).exp()*(1 - arb('1e-8'))*lead/b0**(N - 2)          # e^{-3b}(1-1e-8) 2/((2N)^(N-2)(N-2)!)  (decreasing in beta)
            rhs = (cm + cp)*b0**(6 - N)                                    # (c'_- + c'_+) beta^(6-N)          (increasing in beta)
            ok2 = bool(lhs - rhs > 0); allok = allok and ok2
            extra = ''
            if N == 5: extra = f";  quoted: e^(-3b)(1-1e-8)/3000 = {lhs.str(5)}, (c'_-+c'_+) beta = {rhs.str(5)}, factor {(lhs/rhs).str(4)}"
            if N == 3: extra = f";  quoted: c'_-+c'_+ = {(cm + cp).str(4)} against 1/(648 beta_0) = {(1/(648*b0)).str(5)}; beta^3/1296 >= {(1/(1296*b0)).str(5)} beta^4"
            print(f"    d={d}: m_- - m_+ >= beta^(N-2) [e^(-3b)(1-1e-8) 2/((2N)^(N-2)(N-2)!) - (c'_-+c'_+) beta^(6-N)] > 0 on (0, beta_0(d)): {ok2}{extra}")
    print("all certified:", allok)

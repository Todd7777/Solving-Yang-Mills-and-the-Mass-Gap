# Fourth order (Part IV, Section 12): the predicted beta-series of the gap,
#     m_C = 4 log(1/u) - log Upsilon_C - h_C u^4 + O(beta^5),   h_+ = 26d - 70,  h_- = 10d - 22,
# for U(1), SU(2), ..., SU(5) with the (anti)fundamental representation, in exact rational arithmetic (the Bessel--Toeplitz series of
# ../third-order/upsilon_series.py, written by the referee of Section 11).  The gap is the smaller of the two sectors: C = + for
# SU(N) (N >= 2; for SU(2) there is only one sector), C = - for U(1) (Upsilon_- > Upsilon_+ at order beta^2).
# Printed: the coefficients c_k of m_C - 4 log(1/(u_1 beta)) = sum_k c_k beta^k (u = u_1 beta (1 + O(beta))), and those of
# m_C + 4 log u = -log Upsilon_C - h_C u^4 as a series in u.
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'third-order'))
from upsilon_series import F_G, deriv, inv, log_series, mul, add, scal, upsilons, in_beta, mlog

DEG = 12

def u_series(N, deg=DEG):
    """u(beta) = F'(x)/(N F(x)), x = beta/N (U(1): N = 1, F = I_0), as a series in beta."""
    F, G = F_G(N, deg)
    u = scal(Fr(1, N), mul(deriv(F), inv(F, deg), deg))
    return in_beta(u, N)

def compose(p, q, deg):
    """p(q(t)) for q[0] = 0."""
    out = [Fr(0)]*(deg + 1); powq = [Fr(1)] + [Fr(0)]*deg
    for c in p[:deg + 1]:
        out = add(out, scal(c, powq)); powq = mul(powq, q, deg)
    return out

def revert(q, deg):
    """inverse series of q(t) = q1 t + ..., q1 != 0."""
    r = [Fr(0), 1/q[1]] + [Fr(0)]*(deg - 1)
    for n in range(2, deg + 1):
        c = compose(q, r, n)[n]
        r[n] = -c/q[1]
    return r

def gap_series(N, C, d, order=4, deg=DEG):
    u = u_series(N, deg)
    assert u[0] == 0 and u[1] != 0
    v = [x/u[1] for x in u[1:]] + [Fr(0)]                      # u = u1 beta (1 + ...): v = u/(u1 beta)
    m4 = scal(-4, log_series(v[:deg], deg - 1))                 # -4 log(u/(u1 beta))
    Up, Um = upsilons(N, deg)
    if N == 2:
        U = scal(Fr(1, 2), in_beta(Up, N))                      # iota = 1: Upsilon = Var(chi)
    else:
        U = in_beta(Up if C == '+' else Um, N)
    ml = mlog(U, len(U) - 1)
    h = (26*d - 70) if C == '+' else (10*d - 22)
    u4 = mul(mul(u, u, deg), mul(u, u, deg), deg)
    n = min(len(m4), len(ml), len(u4)) - 1
    tot = [m4[k] + ml[k] - h*u4[k] for k in range(n + 1)]
    # in terms of u: -log Upsilon_C(beta(u)) - h u^4
    binv = revert(u[:n + 1], n)
    mlu = compose(ml[:n + 1], binv, n)
    tot_u = [mlu[k] - (h if k == 4 else 0) for k in range(n + 1)]
    return u[1], tot[:order + 1], tot_u[:order + 1], h

if __name__ == "__main__":
    cases = [(1, '-', 'U(1)'), (1, '+', 'U(1)'), (2, '+', 'SU(2)'), (3, '+', 'SU(3)'), (3, '-', 'SU(3)'), (4, '+', 'SU(4)'),
             (5, '+', 'SU(5)')]
    for d in (3, 4):
        print(f"d = {d}:  h_+ = {26*d-70}, h_- = {10*d-22}")
        for N, C, name in cases:
            u1, cb, cu, h = gap_series(N, C, d)
            print(f"  {name:6s} C{C}: u = ({u1}) beta (1+...);  m - 4 log(1/(u1 beta)) = "
                  + ' + '.join(f'({c}) b^{k}' for k, c in enumerate(cb) if c and k)
                  + f"   [beta^4: {float(cb[4]):+.6f}]")
            print(f"  {'':6s}     m + 4 log u = " + ' + '.join(f'({c}) u^{k}' for k, c in enumerate(cu) if c and k))
        sys.stdout.flush()

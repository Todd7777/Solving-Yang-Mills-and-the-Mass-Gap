# Part IV, the upper bound at practical beta: the tails T_Q, T_p of Lemma mg:lem:twopoint with the sets Y of its proof counted exactly
# (rel_closed_count.py) and the five-face caps bounded by characters.  Usage: python3 upper_exact.py [--group N] [d [beta ...]]
"""T'_S := (1+eps)^|C| e^{theta |P_S|} [ N^cap sigma_5 e^{5(K+1)theta} + sum_{n=5}^{n0} (N_n - [n=5] N^cap) w^n + (1+1/K)^J' sum_{n>n0} q_*^n ],
w = eps e^{(K+1)theta}, q_* = eK w, with theta the constant of the ratio bound: x/(1-x) as in Lemma mg:lem:twopoint ('crude'), theta_5 =
x^5/(1-x) (Lemma mg:lem:closed: any periodic lattice with sides >= 3), or the theta of Remark mg:rem:nextorder (vi) ('closed': infinite volume, d = 3,
closed polymers of up to 14 plaquettes counted exactly, crude from 15 on); sigma_5 = sum_{r != 0} d_r^2 u_r^5 <= (2-iota) N^2 u^5 + u_2^4 (e^beta - 1), u_2 <= e^beta - 1 - beta; for
(S, C) = (S_Q, Q) (|C| = 6, |P_S| = 24d - 42, J' = 24(d-2)) and (S_A, {p0}) (|C| = 1, |P_S| = 8d - 11, J' = K).  Then D'_2 = log(V_2/L_2) of
Theorem mg:thm:upper2 with T'_Q, T'_p in place of T_Q, T_p, and its error D'_2 + log Ups against the third-order term beta^3/432 of SU(3)."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
for sub in ('second-order', 'upper-bound'): sys.path.insert(0, os.path.join(here, '..', sub))
sys.path.insert(0, here)
from mpmath import mp, mpf, e, exp, log, sqrt, expm1, fsum
import third_order_effective as T
from constants_upper import beta0
mp.dps = 60
su3 = dict(T.GROUPS)['SU(3)']
# N_n(S) for n = 5..n0 (rel_closed_count.py), and the number of caps among the five-plaquette sets
COUNTS = {3: {'Q': ({5: 6, 6: 12, 7: 0, 8: 0, 9: 30, 10: 93, 11: 144, 12: 162}, 6), 'A': ({5: 2, 6: 8, 7: 0, 8: 0, 9: 10, 10: 47, 11: 60, 12: 96}, 2)},
          4: {'Q': ({5: 18, 6: 60, 7: 0, 8: 24, 9: 262, 10: 1143, 11: 2616}, 18), 'A': ({5: 4, 6: 32, 7: 0, 8: 0, 9: 60, 10: 506, 11: 760}, 4)}}
# Z^3 to 10 and Z^4 to 8: rel_closed_count.py; Z^3 11-12 and Z^4 9-10: the fifth referee, by two methods (referee-section16e/), and rel_closed_count.py where it finished;
# Z^4 11: admissible_split.c (the fifth referee's Redelmeier code in independent pieces, run_admissible_pieces.sh; orders/admissible_d4_n11.out)
def counts_to(d, n0):
    return {S: ({n: c for n, c in Nn.items() if n <= n0}, ncap) for S, (Nn, ncap) in COUNTS[d].items()}

from third_order_effective import dressed
CLOSED = {3: {10: 10, 11: 11, 12: 24, 13: 0, 14: 70}}
def theta(beta, d, g, kind):
    K = 4*(2*d - 3); eps = expm1(2*beta); x = e**2*K*eps
    if kind == 'crude': return x/(1 - x)
    if kind == 'theta5': return x**5/(1 - x)
    s6 = dressed(beta, g)[1]                  # 'closed': Lemma mg:lem:theta7 with the closed polymers of 10 to 14 plaquettes counted exactly
    return 2*(d - 2)*e**6*s6 + x**15/(1 - x) + fsum(c*(e*eps)**n for n, c in CLOSED[d].items())

def tails(beta, d, g, counts=COUNTS, kind='closed', ubar=False):
    N, iota, u = g['N'], g['iota'], g['u']
    if ubar: u = ((1 + iota)*beta + beta**2*exp(beta))/(2*N**2)      # u <= ubar (Lemma mg:lem:char(iii)); ubar/beta increases visibly
    K = 4*(2*d - 3); eps = expm1(2*beta); th = theta(beta, d, g, kind)
    w = eps*exp((K + 1)*th); q = e*K*w
    u2 = expm1(beta) - beta
    s5 = (2 - iota)*N**2*u**5 + u2**4*expm1(beta)
    out = {}
    for S, nC, PS, Jp in (('Q', 6, 24*d - 42, 24*(d - 2)), ('A', 1, 8*d - 11, K)):
        Nn, ncap = (counts[d] if d in counts else counts)[S]; n0 = max(Nn)
        body = ncap*s5*exp(5*(K + 1)*th) + fsum((c - (ncap if n == 5 else 0))*w**n for n, c in Nn.items()) \
               + (1 + mpf(1)/K)**Jp*q**(n0 + 1)/(1 - q)
        out[S] = (1 + eps)**nC*exp(th*PS)*body
    return out['Q'], out['A'], eps

def upper_err(beta, d, g, TQ=None, Tp=None, kind='closed', counts=COUNTS):
    N, iota, u, Ups, psi1 = g['N'], g['iota'], g['u'], g['Ups'], g['psi1']
    tau0 = mpf(1 + iota)/2
    TQe, Tpe, eps = tails(beta, d, g, counts=counts, kind=kind)
    TQ = TQe if TQ is None else TQ; Tp = Tpe if Tp is None else Tp
    a2 = (1 + Tp)*psi1 + Tp
    V2 = (1 + Tp)*(Ups + psi1**2) + Tp - max((1 - Tp)*psi1 - Tp, 0)**2
    L2 = (1 - TQ)/(1 + eps**6)*max(Ups - Tp*(1 + psi1)*sqrt(2 - iota)*N*u, 0)**2 - (1 + a2**2)*exp(16*beta)*(mpf(N)**4/tau0**2)*TQ
    return log(V2/L2) + log(Ups), TQ, Tp

def E2(beta, d, g, kind='theta5', counts=COUNTS):
    """an upper bound for the error D'_2 + log Ups that is beta^5 times an increasing function of beta (Remark mg:rem:upperexact):
    log(V_2/Ups) <= T'_p (2 + 3 a + 3 a^2) and -log(L_2/Ups^2) <= B/(1 - B), B := T'_Q + eps^6 + 2 sqrt(2-iota) N u T'_p (1 + a)
    + (1 + a'^2) e^{16 beta} N^4/tau_0^2 T'_Q, with a := sqrt(2-iota) N u >= a_1, a' := (1 + T'_p) a + T'_p >= a_2 and Ups >= 1"""
    N, iota, Ups = g['N'], g['iota'], g['Ups']
    assert Ups >= 1
    tau0 = mpf(1 + iota)/2
    u = ((1 + iota)*beta + beta**2*exp(beta))/(2*N**2)                # ubar >= u: E is increasing in u, and every factor is visibly beta^k times increasing
    TQ, Tp, eps = tails(beta, d, g, counts=counts, kind=kind, ubar=True)
    a = sqrt(2 - iota)*N*u; a2 = (1 + Tp)*a + Tp
    B = TQ + eps**6 + 2*sqrt(2 - iota)*N*u*Tp*(1 + a) + (1 + a2**2)*exp(16*beta)*(mpf(N)**4/tau0**2)*TQ
    return Tp*(2 + 3*a + 3*a**2) + B/(1 - B)

if __name__ == "__main__":
    args = sys.argv[1:]; c3 = mpf(1)/432
    if '--group' in args:                 # SU(N) by sun_group.sun (exact Haar moments), against its third-order term (beta^3/6000 for SU(5))
        i = args.index('--group'); Ng = int(args[i + 1]); del args[i:i + 2]
        from sun_group import sun
        su3 = sun(Ng); c3 = {3: mpf(1)/432, 5: mpf(1)/6000}[Ng]
    q = int(round(float(1/c3)))
    d = int(args[0]) if args else 3
    b0 = beta0(d)
    for n0 in sorted({n for S in COUNTS[d] for n in COUNTS[d][S][0] if n >= 8}):
        cn = {d: counts_to(d, n0)}
        g0 = su3(b0); E = E2(b0, d, g0, counts=cn); err, TQ, Tp = upper_err(b0, d, g0, kind='theta5', counts=cn)
        pts = sorted([mpf(10)**(-30 + j*(30 + mp.log10(b0))/400) for j in range(401)]); pts[-1] = b0
        vals = [E2(p, d, su3(p), counts=cn)/p**5 for p in pts]
        mono = all(x <= y for x, y in zip(vals, vals[1:]))
        print(f"d={d}, n0={n0}: at beta_0: T'_Q = {float(TQ/b0**4):.4g} beta^4, T'_p = {float(Tp/b0**4):.4g} beta^4, error {float(err/b0**4):.4g} beta^4; "
              f"E(beta_0) = {float(E/b0**4):.4g} beta_0^4 (with ubar); E/beta^5 increasing on 401 points of [1e-30, beta_0]: {mono}; "
              f"so the error is <= {float(E/b0**4):.4g} beta^4 on the whole range, against beta^3/{q} >= {float(c3/b0):.4g} beta^4: {E/b0**4 < c3/b0}")
        if E/b0**4 >= c3/b0:
            lo, hi = mpf('1e-7'), b0
            for _ in range(60):
                mid = sqrt(lo*hi)
                if E2(mid, d, su3(mid), counts=cn) < c3*mid**3: lo = mid
                else: hi = mid
            print(f"     E(beta) < beta^3/{q} for beta <= {float(lo):.4g}")

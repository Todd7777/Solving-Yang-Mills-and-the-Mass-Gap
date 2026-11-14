# Referee (Section 16h): the error of the upper bound of Theorem mg:thm:upper2 with the tails of Lemma mg:lem:uppertails, for SU(5), written from the
# formulas of the paper (Lemma mg:lem:uppertails, Remark mg:rem:upperexact (ii), (iv)); nothing is imported from upper_exact.py.
# Usage: python3 upper_su5.py         (a few seconds)
"""T'_S = (1+eps)^|C| e^{theta_5 |P_S|} ( N^cap sigma5bar e^{5(K+1)theta_5} + sum_{n=5}^{n0} N_n w^n + (1+1/K)^{J'} q^{n0+1}/(1-q) ),
theta_5 = x^5/(1-x), x = e^2 K eps, eps = e^{2beta}-1, w = eps e^{(K+1)theta_5}, q = eK w; N_n = admissible sets of n plaquettes that are not caps;
sigma5bar = (2-iota)N^2 ubar^5 + (e^beta-1-beta)^4 (e^beta-1), ubar = ((1+iota)beta + beta^2 e^beta)/(2N^2);  (S,C) = (S_Q,Q): |C| = 6, |P_S| = 24d-42,
J' = 24(d-2); (S_A,{p0}): |C| = 1, |P_S| = 8d-11, J' = K.
E = T'_A (2 + 3a + 3a^2) + B/(1-B),  B = T'_Q + eps^6 + 2a(1+a) T'_A + (1 + a'^2) e^{16 beta} (N^4/tau_0^2) T'_Q,  a = sqrt(2-iota) N ubar,
a' = (1 + T'_A) a + T'_A  (valid when Upsilon >= 1).  Also the actual error D_2 + log Ups of Theorem mg:thm:upper2 with T'_Q, T'_A, the exact u and
Upsilon = U_+ (referee's series, moments_check.py)."""
import os, sys
from mpmath import mp, mpf, e, exp, expm1, log, sqrt, fsum
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
mp.dps = 50
N, iota = 5, 0
tau0 = mpf(1 + iota)/2
# the counts of Remark mg:rem:upperexact (i), (iv): (N^cap, {n: non-cap admissible sets})
COUNTS = {3: {'Q': (6, {5: 0, 6: 12, 7: 0, 8: 0, 9: 30, 10: 93, 11: 144, 12: 162}), 'A': (2, {5: 0, 6: 8, 7: 0, 8: 0, 9: 10, 10: 47, 11: 60, 12: 96})},
          4: {'Q': (18, {5: 0, 6: 60, 7: 0, 8: 24, 9: 262, 10: 1143}), 'A': (4, {5: 0, 6: 32, 7: 0, 8: 0, 9: 60, 10: 506})}}

def beta0(d):
    K = 4*(2*d - 3); return log(1 + 1/(e**2*K*(K + 2)))/2

def tails(beta, d, counts, parts=False):
    K = 4*(2*d - 3); eps = expm1(2*beta); x = e**2*K*eps; th = x**5/(1 - x)
    w = eps*exp((K + 1)*th); q = e*K*w
    ub = ((1 + iota)*beta + beta**2*exp(beta))/(2*N**2)
    s5 = (2 - iota)*N**2*ub**5 + (expm1(beta) - beta)**4*expm1(beta)
    out, pr = {}, {}
    for S, nC, PS, Jp in (('Q', 6, 24*d - 42, 24*(d - 2)), ('A', 1, 8*d - 11, K)):
        ncap, Nn = counts[S]; n0 = max(Nn)
        pre = (1 + eps)**nC*exp(th*PS)
        caps = ncap*s5*exp(5*(K + 1)*th); exact = fsum(c*w**n for n, c in Nn.items()); crude = (1 + mpf(1)/K)**Jp*q**(n0 + 1)/(1 - q)
        out[S] = pre*(caps + exact + crude); pr[S] = (pre*caps, pre*exact, pre*crude)
    return (out, pr) if parts else out

def E(beta, d, counts):
    T = tails(beta, d, counts); TQ, TA = T['Q'], T['A']
    eps = expm1(2*beta); ub = ((1 + iota)*beta + beta**2*exp(beta))/(2*N**2)
    a = sqrt(2 - iota)*N*ub; a2 = (1 + TA)*a + TA
    B = TQ + eps**6 + 2*a*(1 + a)*TA + (1 + a2**2)*exp(16*beta)*(mpf(N)**4/tau0**2)*TQ
    return TA*(2 + 3*a + 3*a**2) + B/(1 - B)

def with11(counts, nQ, nA):
    c = {S: (v[0], dict(v[1])) for S, v in counts.items()}
    c['Q'][1][11] = nQ; c['A'][1][11] = nA
    return c

if __name__ == "__main__":
    from moments_check import weyl_moments, series
    mom, _, _ = weyl_moments(N, 19); us, Ups_s, _ = series(N, mom, 17)
    ev = lambda cs, b: fsum(mpf(x.numerator)/x.denominator*b**k for k, x in enumerate(cs))
    for d in (3, 4):
        b0 = beta0(d); c3 = 1/(6000*b0)
        cnt = COUNTS[d]
        T, pr = tails(b0, d, cnt, parts=True)
        print(f"d={d}, beta_0 = {mp.nstr(b0, 10)}, beta^3/6000 = {mp.nstr(c3, 5)} beta^4; n0 = {max(cnt['Q'][1])}")
        for S in ('Q', 'A'):
            print(f"   T'_{S} = {mp.nstr(T[S]/b0**4, 5)} beta^4 = caps {mp.nstr(pr[S][0]/b0**4, 3)} + exact {mp.nstr(pr[S][1]/b0**4, 4)} + crude tail {mp.nstr(pr[S][2]/b0**4, 4)}  (beta^4)")
        Eb = E(b0, d, cnt)
        # the actual error of Theorem mg:thm:upper2 with the exact u, Upsilon
        u = ev(us, b0); Ups = ev(Ups_s, b0); a1 = sqrt(2)*N*u; TQ, Tp = T['Q'], T['A']
        a2 = (1 + Tp)*a1 + Tp
        V2 = (1 + Tp)*(Ups + a1**2) + Tp - max((1 - Tp)*a1 - Tp, 0)**2
        L2 = (1 - TQ)/(1 + expm1(2*b0)**6)*max(Ups - Tp*(1 + a1)*sqrt(2 - iota)*N*u, 0)**2 - (1 + a2**2)*exp(16*b0)*(mpf(N)**4/tau0**2)*TQ
        print(f"   E(beta_0) = {mp.nstr(Eb/b0**4, 6)} beta^4;  actual D_2 + log Ups = {mp.nstr((log(V2/L2) + log(Ups))/b0**4, 6)} beta^4;  Ups - 1 = {mp.nstr(Ups - 1, 6)}")
        pts = [mpf(10)**(-30 + j*(30 + mp.log10(b0))/400) for j in range(401)]; pts[-1] = b0
        vals = [E(p, d, cnt)/p**5 for p in sorted(pts)]
        print(f"   E/beta^5 increasing on 401 points of [1e-30, beta_0]: {all(x <= y for x, y in zip(vals, vals[1:]))}")
        if d == 4:
            for nQ, nA in ((0, 0), (10**6, 10**6), (10**9, 10**9), (10**12, 10**12), (10**15, 10**15)):
                c11 = with11(cnt, nQ, nA)
                Ev = E(b0, d, c11); T11, pr11 = tails(b0, d, c11, parts=True)
                vals = [E(p, d, c11)/p**5 for p in sorted(pts)]
                print(f"   with N_11 = ({nQ:.0e}, {nA:.0e}) for (S_Q, S_A): E(beta_0) = {mp.nstr(Ev/b0**4, 6)} beta^4 (crude tail of T'_Q now {mp.nstr(pr11['Q'][2]/b0**4, 3)} beta^4);"
                      f" E/beta^5 increasing: {all(x <= y for x, y in zip(vals, vals[1:]))};  < beta^3/6000: {Ev/b0**4 < c3}")
            # the largest N_11 (same value for S_Q and S_A, and S_Q alone) for which E(beta_0(4)) < beta^3/6000
            for mode in ('both', 'Q'):
                lo, hi = mpf(0), mpf(10)**30
                for _ in range(200):
                    mid = (lo + hi)/2
                    ok = E(b0, d, with11(cnt, mid, mid if mode == 'both' else 0)) < c3*b0**4
                    lo, hi = (mid, hi) if ok else (lo, mid)
                K = 20
                rootedQ = (1 + mpf(1)/K)**48*(e*K)**11; rootedA = (1 + mpf(1)/K)**K*(e*K)**11
                print(f"   threshold ({mode}): E(beta_0(4)) < beta^3/6000 as long as N_11 <= {mp.nstr(lo, 4)}; the rooted bounds (1+1/K)^J' (eK)^11 are"
                      f" {mp.nstr(rootedQ, 4)} (S_Q) and {mp.nstr(rootedA, 4)} (S_A)")
        if d == 4:
            # N_11 bounded by the rooted binomial C(J' + 11K, 11) of the proof of Lemma mg:lem:twopoint (Part III, Lemma al:lem:rooted), which the crude
            # tail of Lemma mg:lem:uppertails replaces by its exponential majorant (1+1/K)^J' (eK)^11: no enumeration needed
            from mpmath import binomial
            K = 20; bQ = binomial(48 + 11*K, 11); bA = binomial(K + 11*K, 11)
            c11 = with11(cnt, bQ, bA); Ev = E(b0, d, c11)
            vals = [E(p, d, c11)/p**5 for p in sorted(pts)]
            print(f"   with N_11 := C(J'+11K, 11) = {mp.nstr(bQ, 4)} (S_Q), {mp.nstr(bA, 4)} (S_A), the proved rooted bound: E(beta_0(4)) = {mp.nstr(Ev/b0**4, 6)} beta^4"
                  f" < beta^3/6000 = {mp.nstr(c3, 5)} beta^4: {Ev/b0**4 < c3};  E/beta^5 increasing: {all(x <= y for x, y in zip(vals, vals[1:]))}")

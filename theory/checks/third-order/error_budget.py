# Part IV, Section 16 (Theorem mg:thm:thirdcap): the error budget of the evaluation at beta_* in bound mode.
# Usage: python3 error_budget.py        (the four cases SU(3), SU(5) x d = 3, 4; about 10 minutes)
"""The evaluation of Theorem mg:thm:lower5 at beta_* (third_order_effective.lower3e, bound mode) truncates its infinite series: the fronts at k < KMAX = 60, the clusters and
polymers at m, n < NMAX = 40. Every series has non-negative terms, so truncation only lowers the computed sums. This script bounds each omitted tail by a geometric series,
with majorants that hold for every k, m, n beyond the truncation:
  fronts:     N^v_k <= C (y/s)^k (Lemma mg:lem:frontsdecision: 0.783*4.081^k in Z^2, 0.734*7.355^k in Z^3);
  clusters:   binom(K_S m, m-1) <= e (e K_S)^(m-1)                       (binom(n,k) <= (en/k)^k and (m/(m-1))^(m-1) <= e);
  polymers:   binom(c + T n, n) <= (e (T + c/n))^n,  sum_k binom(4n+4, k) k^2 z^k <= (4n+4)^2 (1+z)^(4n+4);
and reports, for each family, the tail relative to the computed sum (tau). Then it checks, with every sum raised by its tail, the inequalities that hold with the least room:
the decoration constants (xi_F >= P_v Y_0(xi_C) + A_0, xi_C >= 4 P_v Y_0 + 8 A_0, h_J >= C_varpi k P_v Y_0 a^{k/2} e^{xi_F k}, the maximum being at k = 4) and the barred tree
constants at (bbar, abar) (the same two inequalities); these were raised by a relative 1e-15 above their fixed points (INFLATE). The marked constants solve monotone linear
systems whose coefficients the tails raise by a relative <= tau; their solutions rise by at most tau times the condition factor 1/det, also reported. All other hypotheses
have room >= 1e-12 (delta' and F are raised by 1e-12, 1e-9). Rounding: the evaluation runs at >= 40 digits."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from mpmath import mp, mpf, e, exp, binomial, fsum, log
import practical_beta as PB
import third_order_effective as T
from third_order_constants import tree_constants
from junction_constants import counts
sys.path.insert(0, os.path.join(here, '..', 'upper-bound'))
from constants_upper import beta0

KMAX, NMAX = 60, 40
DEC = {2: (mpf('0.783'), mpf('4.081')), 3: (mpf('0.734'), mpf('7.355'))}       # Lemma mg:lem:frontsdecision, D = d - 1

def geom(first, ratio):
    assert ratio < 1, ratio
    return first/(1 - ratio)

def tail_clusters(KS, y, p, M=NMAX):
    """sum_{m >= M} m^p binom(KS m, m-1) y^m <= e y M^p (e KS y)^(M-1) / (1 - (1+1/M)^p e KS y)"""
    r = e*KS*y
    return geom(e*y*M**p*r**(M - 1), (1 + mpf(1)/M)**p*r)

def tail_fronts(D, w, j, K=KMAX):
    """sum_{k >= K} k^j N^v_k w^k with N^v_k <= C Y^k"""
    C, Y = DEC[D]; r = Y*w
    return geom(C*K**j*r**K, (1 + mpf(1)/K)**j*r)

def tail_polymers(c, Tv, b, q_first, q_ratio, M=NMAX):
    """sum_{n >= M} binom(c + Tv n, n) b^n q_n with q_{n+1}/q_n <= q_ratio: first term (e(Tv + c/M) b)^M q_M, ratio e (Tv + c/M) b q_ratio"""
    r = e*(Tv + mpf(c)/M)*b
    return geom(r**M*q_first, r*q_ratio)

def budget(group, d, terms, closed, cf, tag):
    PB.GROUP = group; PB.BOUNDS = True
    PB.load_terms(d, terms, closed, cf)
    b0 = beta0(d)
    D, best, g = PB.deficit(d, b0)
    a, b, xF, xC, hJ, cT, F, F6 = best['a'], best['b'], best['xF'], best['xC'], best['hJ'], best['cT'], best['F'], best['F6']
    Pv, KS, Tv, N = counts(d, 'infinite')
    if T.NOVR and d in T.NOVR:
        for k, v in T.NOVR[d].items(): N[k] = mpf(v)
    if T.NV_DECISION:
        for k in N:
            if k >= 6 and not (k in T.NOVR.get(d, {})): N[k] = min(N[k], T.decision_bound(d - 1, k))
    Dd = d - 1
    gam = mp.expm1(2*b0) + b0**2*exp(b0)
    C = 1 + max(1/best['varpi'], (1 + gam)/g['Ubar_low'])
    def Y(z, p, bc): return fsum(m**p*binomial(KS*m, m - 1)*(bc*exp(z))**m for m in range(1, NMAX))
    def A(xf, j, aa): return exp(cT)*fsum(k**j*N[k]*aa**(mpf(k)/2)*exp(xf*k) for k in range(4, KMAX))
    rel = {}
    # actual decoration constants
    Y0, tY0 = Y(xC, 0, b), tail_clusters(KS, b*exp(xC), 0)
    A0, tA0 = A(xF, 0, a), exp(cT)*tail_fronts(Dd, a**mpf(0.5)*exp(xF), 0)
    rel['Y (actual b)'] = tY0/Y0; rel['A (actual a)'] = tA0/A0
    s1 = xF/(Pv*(Y0 + tY0) + A0 + tA0) - 1
    s2 = xC/(4*Pv*(Y0 + tY0) + 8*(A0 + tA0)) - 1
    hreq = max(C*a**(mpf(k)/2)*k*Pv*(Y0 + tY0)*exp(xF*k) for k in range(4, 12))
    assert all(a**(mpf(k + 1)/2)*(k + 1)*exp(xF*(k + 1)) < a**(mpf(k)/2)*k*exp(xF*k) for k in range(4, 12)), 'the maximum over k is at k = 4'
    s3 = hJ/hreq - 1
    # barred tree constants
    bb, ab = best['bars'][0], best['bars'][1]
    tc = tree_constants(bb, ab, cT, Pv, KS, N, KMAX, NMAX)
    xFb, xCb = tc[0], tc[1]
    Yb, tYb = Y(xCb, 0, bb), tail_clusters(KS, bb*exp(xCb), 0)
    Ab, tAb = A(xFb, 0, ab), exp(cT)*tail_fronts(Dd, ab**mpf(0.5)*exp(xFb), 0)
    rel['Y (bbar)'] = tYb/Yb; rel['A (abar)'] = tAb/Ab
    for p in (1, 2):
        rel[f'Y_{p} (bbar)'] = tail_clusters(KS, bb*exp(xCb), p)/Y(xCb, p, bb)
        rel[f'A_{p} (abar)'] = exp(cT)*tail_fronts(Dd, ab**mpf(0.5)*exp(xFb), p)/A(xFb, p, ab)
    s4 = xFb/(Pv*(Yb + tYb) + Ab + tAb) - 1
    s5 = xCb/(4*Pv*(Yb + tYb) + 8*(Ab + tAb)) - 1
    # B-type sums (f(k) = F k^2 e^{xF k} a^{k-4} for k >= 7): sum_k k^j N_k f(k) e^{xF k}, actual and barred
    for lab, aa, xf in (('actual', a, xF), ('bbar', ab, max(xF, xFb))):     # barred: f uses the actual xi_F, the sum e^{xi_F,bar k}: bound both by the larger
        Bsum = fsum(k**2*N[k]*F*k**2*exp(2*xf*k)*aa**(k - 4) for k in range(7, KMAX))
        tB = F*aa**-4*tail_fronts(Dd, exp(2*xf)*aa, 4)
        rel[f'B_2 ({lab})'] = tB/Bsum
    # polymer sums beyond n = NMAX, relative to the parts they belong to: 'four' (q_n = (n+1)(2d-3); in d=3 Z_{n,4} = 0 for all n >= 3, Remark 16.7(iv), so no tail),
    # the exits and Sigma^(6) (q_n <= max(F,F6) a^{-4} (4n+c)^2 (1+z)^{4n+c}, z = e^{xF} a, c = 4, 6)
    if d == 4:
        t4 = 2*(2*d - 3)*tail_polymers(4*Pv, Tv, b, NMAX + 1, mpf(NMAX + 2)/(NMAX + 1))
        rel['four (n >= 40)'] = t4/best['parts']['four']
    # q_n = sum_{k>=5} X(n,k) f(k) <= Fm a^{-4} sum_{k>=5} binom(L,k) k^2 z^k <= Fm L^2 binom(L,5) e^{5 xF} a (1+z)^(L-5), L = 4n+4 (4n+6 for Sigma^(6)),
    # as binom(L, 5+j) <= binom(L,5) binom(L-5,j); its ratio from n to n+1 is at most ((L+4)/L)^2 ((L+4)/(L-1))^5 (1+z)^4 at n = NMAX
    z = exp(xF)*a; Fm = max(F, F6)
    for lab, c0, base in (('exits (n >= 40)', 4, best['parts']['exits']), ('Sigma^(6) (n >= 40)', 6, best['sig6'])):
        L = 4*NMAX + c0
        q0 = Fm*L**2*binomial(L, 5)*exp(5*xF)*a*(1 + z)**(L - 5)
        qr = (mpf(L + 4)/L)**2*(mpf(L + 4)/(L - 1))**5*(1 + z)**4
        rel[lab] = 2*tail_polymers((4 if c0 == 4 else 6)*Pv, Tv, b, q0, qr)/base
    # the marked constants: minimal solutions of monotone 2x2 systems (third_order_constants.marked), at the barred values; inputs raised by tau raise them by
    # at most a relative tau (1 + 8 + ...)/det, with det, dt below (close to 1), far below the 1e-15 by which they are raised
    Yb1, Yb2 = Y(xCb, 1, bb) + tail_clusters(KS, bb*exp(xCb), 1), Y(xCb, 2, bb) + tail_clusters(KS, bb*exp(xCb), 2)
    Ab1 = A(xFb, 1, ab) + exp(cT)*tail_fronts(Dd, ab**mpf(0.5)*exp(xFb), 1); Ab2 = A(xFb, 2, ab) + exp(cT)*tail_fronts(Dd, ab**mpf(0.5)*exp(xFb), 2)
    det = (1 - Ab1)*(1 - 4*Pv*Yb1) - Pv*Yb1*8*Ab1; dt = (1 - 4*Pv*Yb2)*(1 - Ab2) - 8*Ab2*Pv*Yb2
    rel_det = min(det, dt)
    tau = max(rel.values())
    slack = dict(deco_xF=s1, deco_xC=s2, deco_hJ=s3, bar_xF=s4, bar_xC=s5)
    print(f"{tag}: D/beta^4 = {float(D/b0**4):.5g};  largest relative tail tau = {float(tau):.3g} ({max(rel, key=rel.get)});  room with the tails included: "
          + ', '.join(f'{k} {float(v):.3g}' for k, v in slack.items()) + f";  all > 0: {all(v > 0 for v in slack.values())};  marked systems: min det {float(rel_det):.6f}")
    return tau, slack

if __name__ == "__main__":
    mp.dps = 50
    CL3 = ({10: 10, 11: 11, 12: 24, 13: 0, 14: 70}, 15)
    CL4 = ({10: 60, 11: 66, 12: 320, 13: 0, 14: 1316, 15: 2520, 16: 12552}, 17)
    out = []; PB.CERTIFIED = True     # only counts made by two programs (practical_beta.CERTIFIED)
    for group in (3, 5):
        out.append(budget(group, 3, os.path.join(here, 'orders', 'd3_s9.json'), *CL3, f"SU({group}), d=3"))
        out.append(budget(group, 4, os.path.join(here, 'orders', 'd4_s10x_92.json'), *CL4, f"SU({group}), d=4"))
    print("all cases: largest tau", float(max(t for t, _ in out)), "; smallest room", float(min(min(s.values()) for _, s in out)))

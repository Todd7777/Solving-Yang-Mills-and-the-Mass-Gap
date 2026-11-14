# Part IV: the group data of SU(N) (fundamental representation) from exact Haar moments, for the evaluators of Sections 10-16 (any N >= 2).
# Usage: python3 sun_group.py [N]        (checks against second_order_constants.su2/su3, and prints the series of -log Upsilon)
"""Weight exp((beta/N) Re chi(U)) on SU(N), chi the trace (the normalisation of second_order_constants.su3). The Haar moments
    m(a, b) = E[chi^a conj(chi)^b] = sum f^lam f^mu over partitions lam of a, mu of b with at most N rows that give the same SU(N) representation
(lam - mu constant: Schur-Weyl duality, chi^a = sum_lam f^lam chi_lam, f^lam the number of standard tableaux) are exact integers. With t = beta/N and
    G(p, q) = E[chi^p conj(chi)^q e^{t Re chi}] = sum_{a,b} (t/2)^{a+b} m(a+p, b+q)/(a! b!)
(truncated at a + b <= KMAX; the tail is below (t/2)^KMAX N^KMAX/((KMAX/2)!)^2, about 1e-190 at beta = 4e-4 with KMAX = 40, N = 5):
    <Re chi> = (G(1,0) + G(0,1))/(2 G(0,0)),   <(Re chi)^2> = (G(2,0) + 2G(1,1) + G(0,2))/(4 G(0,0)),   <(Im chi)^2> = -(G(2,0) - 2G(1,1) + G(0,2))/(4 G(0,0)),
u = <Re chi>/N, Upsilon = max(U_+, U_-) with U_+ = 2 Var(Re chi), U_- = 2 <(Im chi)^2> (<Im chi> = 0 by U -> conj(U)), psi_1 = sqrt2 N u if U_+ >= U_-
(the dictionary of second_order_constants.su3). The coefficients are exact rationals, so the data are exact to the working precision."""
import os, sys
from fractions import Fraction
from functools import lru_cache
from math import factorial
from mpmath import mp, mpf, sqrt, log, taylor, mpmathify

def partitions(n, k, maxpart=None):
    """partitions of n into at most k parts, as non-increasing tuples"""
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    if k == 0: return
    for first in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - first, k - 1, first): yield (first,) + rest

def f_lambda(lam):
    """the number of standard Young tableaux of shape lam (hook length formula)"""
    n = sum(lam); cols = [sum(1 for r in lam if r > j) for j in range(lam[0])] if lam else []
    hooks = 1
    for i, r in enumerate(lam):
        for j in range(r): hooks *= (r - j - 1) + (cols[j] - i - 1) + 1
    return factorial(n)//hooks

@lru_cache(maxsize=None)
def by_class(N, a):
    """{normal form nu (lam minus its N-th part in every row): f^lam} over the partitions lam of a with at most N rows"""
    out = {}
    for lam in partitions(a, N):
        c = lam[N - 1] if len(lam) == N else 0
        nu = tuple(x - c for x in lam if x - c > 0)
        out[nu] = out.get(nu, 0) + f_lambda(lam)
    return out

@lru_cache(maxsize=None)
def moment(N, a, b):
    """E[chi^a conj(chi)^b] under Haar measure on SU(N), exactly"""
    A, B = by_class(N, a), by_class(N, b)
    return sum(v*B[nu] for nu, v in A.items() if nu in B)

@lru_cache(maxsize=None)
def coeffs(N, KMAX):
    """the exact coefficients c[(p,q)][k] of t^k in G(p, q), k <= KMAX"""
    out = {}
    for p, q in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]:
        c = [Fraction(0)]*(KMAX + 1)
        for a in range(KMAX + 1):
            for b in range(KMAX + 1 - a):
                m = moment(N, a + p, b + q)
                if m: c[a + b] += Fraction(m, 2**(a + b)*factorial(a)*factorial(b))
        out[(p, q)] = c
    return out

def u2_bound(N, beta):
    """an upper bound for u_2 = sup_{r not in {0, rho, rho-bar}} u_r (SU(N), N >= 3): with t = beta/N and Z >= 1 (Jensen),
    d_r u_r <= sum_{a,b} (t/2)^{a+b} m_r(a,b)/(a! b!), m_r(a,b) the multiplicity of r in rho^a (x) rho-bar^b, at most N^{a+b}/d_r. The irreps with
    m_r(a,b) > 0 for some a + b = 2 are Sym^2, Lambda^2 (Lambda^2 = rho-bar for N = 3), the adjoint and their conjugates, with m = 1 at (2,0), (2,0), (1,1);
    the tail a + b >= 3 is at most R/d_r, R = e^beta - 1 - beta - beta^2/2; every other r has d_r >= N and u_r <= R/d_r^2."""
    beta = mpf(beta); t = beta/N; R = mp.expm1(beta) - beta - beta**2/2
    cands = [(t**2/8, mpf(N*(N + 1))/2), (t**2/4, mpf(N*N - 1))]          # Sym^2, adjoint
    if N >= 4: cands.append((t**2/8, mpf(N*(N - 1))/2))                     # Lambda^2
    return max([(low + R/dr)/dr for low, dr in cands] + [R/N**2])

def sun(N, KMAX=40, sharp_u2=False):
    """the group data of SU(N) at beta, as second_order_constants.su3 (u, Ups, psi1, N, iota, Upm)"""
    C = coeffs(N, KMAX)
    def g(beta):
        beta = mpf(beta)
        extra = max(0, int(-6*mp.log10(abs(beta)))) if beta != 0 else 0
        with mp.workdps(mp.dps + 10 + extra):        # Upsilon - 1 can be O(beta^3) (SU(5)) or smaller: resolve it to the outer precision
            t = beta/N
            G = {pq: sum(mpf(c.numerator)/c.denominator*t**k for k, c in enumerate(cs) if c) for pq, cs in C.items()}
            Z = G[(0, 0)]
            mR = (G[(1, 0)] + G[(0, 1)])/(2*Z)
            Up = 2*((G[(2, 0)] + 2*G[(1, 1)] + G[(0, 2)])/(4*Z) - mR**2)
            Um = 2*(-(G[(2, 0)] - 2*G[(1, 1)] + G[(0, 2)])/(4*Z))
            u = mR/N
        Up, Um, u = +Up, +Um, +u                                  # rounded to the working precision
        if N == 2:                             # chi real: the unit function is chi itself, Upsilon = Var(chi) = U_+/2 (second_order_constants.su2)
            return dict(name='SU(2)', u=u, Ups=Up/2, psi1=2*u, N=2, iota=1)
        out = dict(name=f'SU({N})', u=u, Ups=max(Up, Um), psi1=(sqrt(2)*N*u if Up >= Um else 0), N=N, iota=0, Upm=(Up, Um))
        if sharp_u2: out['u2'] = u2_bound(N, beta)     # in place of e^beta - 1 - beta (Lemma mg:lem:capelement(iv)) in Lemma mg:lem:dressed
        return out
    return g

def series_minus_log_ups(N, order=6, KMAX=24):
    """the Taylor coefficients of -log Upsilon in beta (C-even sector U_+ and C-odd U_-), from the exact moments, as floats at high precision"""
    mp.dps = 60
    g = sun(N, KMAX)
    out = {}
    if N == 2: return {'Ups': taylor(lambda x: -log(g(x)['Ups']), mpf(0), order)}
    for idx, name in ((0, 'U+'), (1, 'U-')):
        out[name] = taylor(lambda x: -log(g(x)['Upm'][idx]), mpf(0), order)
    return out

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..', 'second-order')); sys.path.insert(0, os.path.join(here, '..', 'upper-bound'))
    from second_order_constants import su2, su3
    Ns = [int(sys.argv[1])] if len(sys.argv) > 1 else [2, 3, 4, 5, 6]
    mp.dps = 40
    print("moments E[|chi|^2k], SU(5), k <= 6:", [moment(5, k, k) for k in range(7)], "; E[chi^5] =", moment(5, 5, 0), "; E[chi^3] (SU(3)) =", moment(3, 3, 0))
    for N, ref in ((2, su2), (3, su3)):
        if N not in Ns: continue
        for beta in ('1e-6', '1e-4', '4e-4', '0.01', '0.1'):
            gs, gr = sun(N)(mpf(beta)), ref(mpf(beta))
            du = abs(gs['u']/gr['u'] - 1); dU = abs(gs['Ups']/gr['Ups'] - 1)
            print(f"SU({N}) beta={beta}: |u/u_ref - 1| = {float(du):.1e}, |Ups/Ups_ref - 1| = {float(dU):.1e}, sector {'+' if gs['psi1'] else '-'} vs {'+' if gr['psi1'] else '-'}")
    for N in Ns:
        s = series_minus_log_ups(N)
        print(f"SU({N}): " + ';  '.join(f"-log {name} = " + ' + '.join(f'{mp.nstr(c, 12)} b^{k}' for k, c in enumerate(cs) if abs(c) > mpf(10)**-40)
                                         for name, cs in s.items()))

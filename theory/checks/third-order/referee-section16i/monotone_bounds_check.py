# Referee sanity check (ninth pass on Section 16) of Theorem mg:thm:thirdcap(b): the bound-mode quantities are beta^k times increasing functions.
# This is NOT a proof (the claim is "by inspection"; the README argues each one); it checks on 2000 log-spaced points of [1e-24, beta_*] that
#   ubar/beta, beta/ulow (= 2N^2 e^beta/(1+iota)), u2bar/beta^2, s6bar/beta^6, nB/beta^3, c2/beta^4, Delta/beta^6, b2bar/beta^2, theta/beta^6 (d = 3, 4)
# increase, and that e^{-2 beta} - abar decreases, with the referee's implementation (my_eval.group_bounds).
# Usage: python3 monotone_bounds_check.py
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..'))
from mpmath import mp, mpf, e, exp, expm1, fsum, nstr
import my_eval as ME

mp.dps = 80
def qty(N, beta, d):
    g = ME.group_bounds(N, beta)
    K = 4*(2*d - 3); eps = g['eps']; x = e**2*K*eps
    closed, CF = ME.SETTINGS[d][1], ME.SETTINGS[d][2]
    th = 2*(d - 2)*e**6*g['s6'] + fsum(c*(e*eps)**n for n, c in closed.items()) + x**CF/(1 - x)
    lam = exp(2*(K + 1)*th); a = lam*g['ubar']
    u2, ulow = g['u2'], g['ulow']
    nB = (u2/ulow)**2*eps + mp.sqrt(2)*(N/ulow)*mp.sqrt(u2**5*expm1(beta))          # Lemma mg:lem:dressed, iota = 0
    c2 = (u2/ulow)**4*exp(beta)
    b2 = lam**2*(2*g['s6'] + g['s6']**2)/ulow**4
    return {'ubar/b': g['ubar']/beta, 'b/ulow': beta/ulow, 'u2/b^2': u2/beta**2, 's6/b^6': g['s6']/beta**6, 'nB/b^3': nB/beta**3,
            'c2/b^4': c2/beta**4, 'Delta/b^6': g['Delta']/beta**6, 'b2/b^2': b2/beta**2, 'theta/b^6': th/beta**6, '-(e^-2b - a)': -(exp(-2*beta) - a)}

if __name__ == "__main__":
    allok = True
    for N in (3, 5):
        for d in (3, 4):
            bs = mpf(ME.SETTINGS[d][4])
            pts = [mpf(10)**(-24 + i*(24 + mp.log10(bs))/1999) for i in range(2000)]; pts[-1] = bs
            vals = [qty(N, p, d) for p in pts]
            bad = [k for k in vals[0] if not all(v[k] <= w[k] for v, w in zip(vals, vals[1:]))]
            allok = allok and not bad
            print(f"SU({N}), d={d}: all ten quantities increasing on 2000 points of [1e-24, beta_*]: {not bad}" + (f"  (not: {bad})" if bad else ''))
    print("all increasing:", allok)

# Section 11 of Part IV: exact transfer-matrix spectra of Z_n gauge theory in d = 3 (spatial L x L torus, winding-zero sector,
# dual plaquette variables; the transfer matrix of ../second-order/zn_dual_spectra.py, written by the referee of Section 10).
# For Z_3 and Z_5 the plaquette variance -log Upsilon_C has a beta^3 term, so the claim m = 4 log(1/u) - log Upsilon_C + O(beta^4) of
# Corollary mg:cor:third is tested at third order: the residual after subtracting the exact -log Upsilon_C is about -8 u^4, whereas
# after subtracting only its expansion to second order a term c_3 beta^3 remains, with c_3 the third-order coefficient of -log Upsilon_C.
# The tori 3x3 have sides < 6 (outside the hypotheses of Theorem mg:thm:lower3); the referee's 4x4 and 5x5 runs (Section 10) agree.
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'second-order'))
from zn_dual_spectra import run
from mpmath import mp, mpf, taylor, log, cos, sin, pi, exp, fsum
mp.dps = 40

def mlogups_series(n, C, order=4):
    """Taylor coefficients of -log Upsilon_C(beta) for Z_n (weight exp(beta cos(2 pi m/n)))."""
    def F(beta):
        th = [2*pi*m/n for m in range(n)]; w = [exp(beta*cos(t)) for t in th]; Z = fsum(w)
        phi = [mp.sqrt(2)*(cos(t) if C == '+' else sin(t)) for t in th]
        m1 = fsum(wi*p for wi, p in zip(w, phi))/Z; m2 = fsum(wi*p*p for wi, p in zip(w, phi))/Z
        return -log(m2 - m1**2)
    return [float(c) for c in taylor(F, 0, order)]

if __name__ == "__main__":
    for n, L in ((3, 3), (5, 3)):
        ser = {C: mlogups_series(n, C) for C in ('+', '-')}
        for C in ('+', '-'):
            print(f"Z_{n}, C{C}: -log Upsilon_C = " + ' + '.join(f'({c:+.6f}) b^{k}' for k, c in enumerate(ser[C]) if k))
        for beta in (0.04, 0.02, 0.01):          # below 0.01 the residual (~1e-10) reaches the precision of the eigenvalues
            u, out = run(n, L, beta, [('+', '0'), ('-', '0')], iters=60, tol=1e-15, verbose=False)
            for (C, mom, m, pred, diff, it) in out:
                s = ser[C]; trunc2 = s[1]*beta + s[2]*beta**2
                r2 = (m - 4*np.log(1/u)) - trunc2
                print(f"  Z_{n} {L}x{L} beta={beta:<6} C{C}: m-4log(1/u)={m-4*np.log(1/u):+.10e}  residual (exact Upsilon)/u^4={diff/u**4:+.3f}"
                      f"   residual (Upsilon to 2nd order)/beta^3={r2/beta**3:+.5f}  (c_3={s[3]:+.5f})")
            sys.stdout.flush()

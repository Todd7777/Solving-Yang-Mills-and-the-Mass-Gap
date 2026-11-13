# Referee (Section 16, fourth pass): the scaling of Remark mg:rem:nextorder (vii), checked with the referee's own evaluator (interp_frozen.evaluate).
# Usage: python3 scaling_frozen.py [beta_* ('b0' or a number)]
"""At beta_* the parameters are those of the authors' search (lower3e via interp_frozen.capture: the smaller delta' of varpi = b^2, b^(5/2)).
For beta = t beta_* and kappa = beta_*: varpi = t^p varpi_*, F = (1+kappa) t F_*, F_6 = (1+kappa) t F_6*, delta' = (1+kappa)^2 t^4 delta'_*, the same
xi_F, xi_C, h_J, (bbar, abar) and kind of majorant, c_J = e^{h_J}/(1-omega) recomputed, the marked constants the minimal solution at beta.
Every condition is evaluated by interp_frozen.evaluate (the decoration conditions for k = 4..59, the six inequalities, the barred constants at the
fixed (bbar, abar), the hypotheses of Theorem mg:thm:lower5), and D(beta)/beta^4 is compared with (1+kappa)^3 D(beta_*)/beta_*^4.
Also printed: the relative margins of the F- and F_6-conditions (about kappa = 4.0e-4 at t = 1, growing to 3e-2 as t -> 0). That margin comes from
Ubar - a = 1 + beta/9 + ..., which increases with beta, and Ubar_* - 1 = kappa/6 (monotone_check.py); the bound 1 <= Ubar <= 1 + kappa quoted in the
remark would not suffice by itself (it leaves a deficit of order kappa a_* = 9e-9, larger than the 1e-9 margin built into F_*, F_6*)."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import interp_frozen as I
from mpmath import mp, mpf, log, exp

if __name__ == "__main__":
    mp.dps = 60
    bs = I.beta0(I.d) if len(sys.argv) < 2 or sys.argv[1] == 'b0' else mpf(sys.argv[1])
    g = dict(I.T.GROUPS)['SU(3)'](bs)
    best = None
    for vp in (2, mpf(5)/2):
        try: P, lo = I.capture(bs, g, vp)
        except AssertionError as ex: continue
        if best is None or lo['delta'] < best[1]['delta']: best = (P, lo)
    P, lo = best
    Ds = (lo['loss'] - log(g['Ups']))/bs**4
    k = bs; p = P['varpi_pow']
    q0 = I.basic(bs, g)
    print(f"beta_* = {float(bs):.10g}: D_*/beta_*^4 = {float(Ds):.6g}, p = {float(p):g}, bars {P['bars']}, Ubar_* - 1 = {float(q0['Ubar']-1):.4e} (kappa = {float(k):.4e}),"
          f" a_* = {float(q0['a']):.4e}, delta'_* = {float(P['dl']):.4e}")
    mp.dps = 130
    tol = mpf(10)**-40
    for ts in ('1', '0.999', '0.99', '0.95', '0.9', '0.7', '0.5', '0.3', '0.1', '1e-2', '1e-3', '1e-4', '1e-6', '1e-8', '1e-12', '1e-20'):
        t = mpf(ts); beta = t*bs; gb = dict(I.T.GROUPS)['SU(3)'](beta)
        Pt = dict(P); Pt.update(varpi=t**p*P['varpi'], F=(1+k)*t*P['F'], F6=(1+k)*t*P['F6'], dl=(1+k)**2*t**4*P['dl'], cT=None, C=None)
        out, ex = I.evaluate(beta, gb, Pt)
        bad = [n for n, (l, r) in out.items() if n != 'C_frozen' and not (l <= r or (n.startswith('six_') and l - r <= tol*abs(r)))]
        marg = {n: (r - l)/abs(r) for n, (l, r) in out.items() if n != 'C_frozen' and not n.startswith('six_')}
        worst = min(marg.items(), key=lambda kv: kv[1])
        bound = (1 + k)**3*Ds; Dq = ex['D']/beta**4
        # the F- and F6-conditions with the right side recomputed as if Ubar_* = 1 + kappa were attained at beta_* (the remark's stated bound):
        print(f"t = {ts:>6} (beta = {float(beta):.4e}): violated: {bad if bad else 'none'};  smallest margin {float(worst[1]):.3e} ({worst[0]});"
              f"  margins exits {float(marg['exits']):.3e}, six {float(marg['six']):.3e}, sigma {float(marg['sigma']):.3e}, dec_hJ {float(marg['dec_hJ']):.3e},"
              f" dec_xiF {float(marg['dec_xiF']):.3e};  D/beta^4 = {float(Dq):.6g} <= {float(bound):.6g}: {Dq <= bound}")
        sys.stdout.flush()

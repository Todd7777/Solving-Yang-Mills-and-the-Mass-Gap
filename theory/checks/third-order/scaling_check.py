# Section 16 of Part IV (Remark mg:rem:nextorder (vii)): the scaling of the admissible parameters below a computed point, d=3, S''_9.
# Usage: python3 scaling_check.py [--d4 TERMS.json] [beta_* [JOBS]]        (default beta_* = beta_0(3); about 10 minutes with 3 jobs)
"""At beta_* the parameters of Theorem mg:thm:lower5 are found as usual (practical_beta.deficit). For beta = t beta_*, 0 < t <= 1, and k := beta_* take
    varpi = t^p varpi_*  (varpi_* = b(beta_*)^p),   F = (1+k) t F_*,   F_6 = (1+k) t F_6*,   delta' = (1+k)^2 t^4 delta'_*,
and keep xi_F, xi_C, h_J and (bbar, abar) of beta_* (these constants come out of lower3e raised by a relative 1e-15, so that their conditions hold
strictly; c_J = e^{h_J}/(1 - omega) is at most its value at beta_*, as omega = O(t^(4-p))). Every term of the sums in the hypotheses then has at
least the order in t of its side: t^4 in Sigma_4^(4) (two: varpi b_2 ~ t^(p+2); exits: a f(6) ~ t^4 and more; one plaquette ~ t^5; the majorant ~ t^9),
t in the conditions on F and F_6 (need, Sigma^(6) ~ t), and the decoration conditions improve (C_varpi a^{k/2} Y_0 ~ t^(3-p)). The factors 1 + k absorb
the decrease of Ubar = 1 + beta/6 + ... on the right sides. As 8(K+1)theta/beta^6 and log(Ubar/Ups)/beta^4 increase with beta, then
D(beta)/beta^4 <= (1+k)^3 D(beta_*)/beta_*^4 (the third factor covers delta' <= (1 + delta') log(1 + delta')). This script checks all hypotheses (lower3e
with fixed parameters) and that inequality at 17 values of t; it is a check of the scaling, which is argued term by term in the remark."""
import os, sys
from multiprocessing import Pool
from mpmath import mp, mpf, log
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import practical_beta as PB
import third_order_effective as T

# the setting: d=3 with S''_9 (default), or d=4 with the set given by --d4 FILE (closed polymers of Z^4 with 10, 11 plaquettes counted exactly)
if '--d4' in sys.argv:
    i = sys.argv.index('--d4'); D4FILE = sys.argv[i + 1]; del sys.argv[i:i + 2]
    d = 4
    PB.load_terms(d, D4FILE, {10: 60, 11: 66}, 12)
else:
    d = 3
    PB.load_terms(d, os.path.join(here, 'orders', 'd3_s9.json'), {10: 10, 11: 11, 12: 24, 13: 0, 14: 70}, 15)
sys.path.insert(0, os.path.join(here, '..', 'upper-bound'))
from constants_upper import beta0
su3 = dict(T.GROUPS)['SU(3)']

def params(beta_s):
    """the admissible parameters at beta_* (the smaller delta' of varpi = b^2, b^(5/2)) and its deficit"""
    g = su3(beta_s); best = None
    for vp in (2, mpf(5)/2):
        try: lo = T.lower3e(beta_s, d, g, varpi_pow=vp, chars=True)
        except AssertionError: continue
        if best is None or lo['delta'] < best['delta']: best = lo; best['p'] = vp
    best['D'] = best['loss'] - log(g['Ups'])
    return best

def scaled(args):
    beta_s, P, t = args
    k = beta_s; beta = t*beta_s; g = su3(beta)
    fixed = dict(varpi=P['varpi']*t**P['p'], xF=P['xF'], xC=P['xC'], hJ=P['hJ'], F=(1+k)*t*P['F'], F6=(1+k)*t*P['F6'], dl=(1+k)**2*t**4*P['delta'])
    try:
        lo = T.lower3e(beta, d, g, chars=True, fixed=fixed, bbar_grid=[P['bars'][0]], abar_grid=[P['bars'][1]])
    except AssertionError as err:
        return t, None, str(err)[:300]
    return t, (lo['loss'] - log(g['Ups']))/beta**4, None

if __name__ == "__main__":
    beta_s = mpf(sys.argv[1]) if len(sys.argv) > 1 else beta0(d)
    jobs = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    P = params(beta_s); Ds = P['D']/beta_s**4
    mp.dps = 130                     # D ~ beta^4 down to beta ~ 1e-24 needs about 100 digits (log(1 + delta'), log Ubar - log Ups)
    print(f"beta_* = {float(beta_s):.10g}: D/beta^4 = {float(Ds):.6g}, varpi = b^{float(P['p']):g}, (bbar, abar) = {P['bars']}")
    ts = [mpf(1), mpf('0.999'), mpf('0.99'), mpf('0.95'), mpf('0.9'), mpf('0.7'), mpf('0.5'), mpf('0.3'), mpf('0.1'), mpf('0.03'), mpf('1e-2'), mpf('1e-3'),
          mpf('1e-4'), mpf('1e-6'), mpf('1e-8'), mpf('1e-12'), mpf('1e-20')]
    with Pool(jobs) as pool:
        res = pool.map(scaled, [(beta_s, P, t) for t in ts], chunksize=1)
    k = beta_s; ok = True
    for t, r, err in res:
        bound = (1 + k)**3*Ds
        good = err is None and r <= bound
        ok = ok and good
        print(f"  t = {float(t):.3g} (beta = {float(t*beta_s):.4g}): " + (f"all hypotheses hold; D/beta^4 = {float(r):.6g} <= (1+k)^3 D_*/beta_*^4 = {float(bound):.6g}: {good}"
                                                                    if err is None else f"FAILS: {err}"))
    print(f"all scaled points admissible with D/beta^4 <= (1+k)^3 D(beta_*)/beta_*^4: {ok};  then D < beta^3/432 for all beta <= beta_* as long as "
          f"(1+k)^3 D_*/beta_*^4 = {float((1+k)**3*Ds):.5g} < 1/(432 beta_*) = {float(1/(432*beta_s)):.5g}: {(1+k)**3*Ds < 1/(432*beta_s)}")

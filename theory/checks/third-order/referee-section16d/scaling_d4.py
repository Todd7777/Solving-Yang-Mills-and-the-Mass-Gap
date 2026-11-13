# Referee (Section 16, fourth pass): the scaling of Remark mg:rem:nextorder (vii) in the other settings where the text now uses it, d = 4:
# Remark mg:rem:practical (iii) (M_5, theta_7) at beta_* = 2.6e-6, and Remark mg:rem:nextorder (iii) (S_6, theta_10) at beta_* = 1.2e-5.
# Usage: python3 scaling_d4.py        (uses the authors' lower3e(fixed=...) of third_order_effective.py; the d = 3 check is scaling_frozen.py)
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..'))
from mpmath import mp, mpf, log
import practical_beta as PB
import third_order_effective as T
su3 = dict(T.GROUPS)['SU(3)']

def params(d, beta_s):
    g = su3(beta_s); best = None
    for vp in (2, mpf(5)/2):
        try: lo = T.lower3e(beta_s, d, g, varpi_pow=vp, chars=True)
        except AssertionError: continue
        if best is None or lo['delta'] < best['delta']: best = lo; best['p'] = vp
    best['D'] = best['loss'] - log(g['Ups'])
    return best

def check(label, d, beta_s):
    mp.dps = 60
    P = params(d, beta_s); Ds = P['D']/beta_s**4; k = beta_s
    print(f"{label}: beta_* = {float(beta_s):.3g}, D_*/beta_*^4 = {float(Ds):.6g}, p = {float(P['p']):g}, majorant {P['bars']}"); sys.stdout.flush()
    mp.dps = 130
    r = 1 + mpf(10)**-15
    bb, ab = P['bars'][0], P['bars'][1]
    for ts in ('1', '0.9', '0.5', '0.1', '1e-3', '1e-6', '1e-12', '1e-20'):
        t = mpf(ts); beta = t*beta_s; g = su3(beta)
        fixed = dict(varpi=P['varpi']*t**P['p'], xF=r*P['xF'], xC=r*P['xC'], hJ=r*P['hJ'], F=(1+k)*t*P['F'], F6=(1+k)*t*P['F6'], dl=(1+k)**2*t**4*P['delta'])
        try:
            lo = T.lower3e(beta, d, g, chars=True, fixed=fixed, bbar_grid=[bb], abar_grid=[ab])
            Dq = (lo['loss'] - log(g['Ups']))/beta**4
            print(f"   t = {ts:>6}: all hypotheses hold; D/beta^4 = {float(Dq):.6g} <= (1+k)^3 D_*/beta_*^4 = {float((1+k)**3*Ds):.6g}: {Dq <= (1+k)**3*Ds}")
        except AssertionError as err:
            print(f"   t = {ts:>6}: FAILS {str(err)[:200]}")
        sys.stdout.flush()

if __name__ == "__main__":
    # Remark mg:rem:practical (iii): theta_7, M_5 with R a maximum (practical_beta's defaults)
    T.CLOSED_FROM = 8; T.EXACT = None; T.CORNERS = None
    check("d=4, Remark practical (iii) [M_5, theta_7]", 4, mpf('2.6e-6'))
    # Remark mg:rem:nextorder (iii): S_6 with theta_10
    T.CLOSED_FROM = 10; PB.order6(4)
    check("d=4, Remark nextorder (iii) [S_6, theta_10]", 4, mpf('1.2e-5'))

# Referee check (ninth pass on Section 16) of Theorem mg:thm:thirdcap(b): the scaling of Remark mg:rem:nextorder(vii) in bound mode, evaluated with the
# referee's evaluator (my_eval.evaluate: extended sums, correctly raised marked constants) at beta = t beta_*, with
#   varpi = t^p varpi_*, F = (1+k) t F_*, F_6 = (1+k) t F_6*, delta' = (1+k)^2 t^4 delta'_*, k = beta_*, and the same xi_F, xi_C, h_J, (bbar, abar).
# Prints the smallest room of the hypotheses and D/beta^4 against (1+k)^3 D(beta_*)/beta_*^4.
# Usage: python3 scaled_check.py GROUP d
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..'))
from mpmath import mp, mpf, nstr
import my_eval as ME

if __name__ == "__main__":
    G, d = int(sys.argv[1]), int(sys.argv[2])
    mp.dps = 60
    P = ME.authors_params(G, d)
    bs = P['beta']; k = bs
    mp.dps = 140
    r0 = ME.evaluate(G, d, P)
    Ds = r0['D']/bs**4; bound = (1 + k)**3*Ds
    print(f"SU({G}), d={d}: D(beta_*)/beta_*^4 = {nstr(Ds, 10)}, (1+k)^3 D_*/beta_*^4 = {nstr(bound, 10)}")
    ok = True
    for t in ('1', '0.9', '0.5', '0.1', '1e-3', '1e-6', '1e-12', '1e-20'):
        t = mpf(t)
        Q = dict(P, beta=t*bs, varpi=P['varpi']*t**(mpf(5)/2), F=(1 + k)*t*P['F'], F6=(1 + k)*t*P['F6'], dl=(1 + k)**2*t**4*P['dl'])
        r = ME.evaluate(G, d, Q)
        mins = min(r['rooms'].values()); mk = min(min(r['rooms_my_a']), min(r['rooms_my_b']))
        Db = r['D']/Q['beta']**4
        good = mins > 0 and mk > 0 and Db <= bound
        ok = ok and good
        print(f"  t = {nstr(t, 3):>6s}: smallest room of the hypotheses {nstr(mins, 4)} ({min(r['rooms'], key=r['rooms'].get)}), of my six inequalities {nstr(mk, 3)};"
              f" D/beta^4 = {nstr(Db, 10)} <= {nstr(bound, 10)}: {good}")
    print(f"  all scaled points pass with the referee's evaluator: {ok}")

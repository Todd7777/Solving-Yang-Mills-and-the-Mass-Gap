# Referee check (eleventh pass on Section 16, Corollary mg:cor:classic), item 1: the inputs of bound mode (Lemma mg:lem:boundmode) and of the upper bound,
# evaluated with the EXACT group data of SU(2) and U(1) (modified Bessel functions), as a sanity check of the lemmas the corollary invokes (Numerical; the
# lemmas themselves are proofs). For w_beta = exp((beta/N) Re chi_rho):
#   SU(2): e^{(beta/2) chi} = (2/beta) sum_j (2j+1) I_{2j+1}(beta) chi_j, so u_j = I_{2j+1}/I_1, u = I_2/I_1, u_2 = sup_{j>=1} u_j = I_3/I_1,
#          sigma_k = sum_{n>=2} n^2 (I_n/I_1)^k, kappa = d_1 u_1 = 3 I_3/I_1 (chi_{1/2}^2 = chi_0 + chi_1), Upsilon = 1 + kappa - N^2 u^2 (Lemma mg:lem:Upsilon(i));
#   U(1):  u_n = I_n/I_0, u = I_1/I_0, u_2 = I_2/I_0, sigma_k = 2 sum_{n>=1} (I_n/I_0)^k, alpha = phi_{rho rho} = int f = 0, gamma = phi_{rho rhobar} = I_2/I_0,
#          kappa = alpha + gamma, Upsilon_+ = 1 + kappa - 2u^2, Upsilon_- = 1 + alpha - gamma.
# Checked at beta in {beta_0(3), beta_0(4), 1e-4, 1e-6}:
#   (1) u_* = u (Corollary mg:cor:window): u_r <= u for every r != 0;  (2) u_low <= u <= ubar (Lemma mg:lem:char(iii));  (3) u_2 <= e^b - 1 - b (Lemma
#   mg:lem:capelement(iv));  (4) sigma_6 <= sigma6bar, sigma_5 <= sigma5bar;  (5) e^{-2b} <= Upsilon(_C) <= e^b (Lemma mg:lem:upsbounds(i)), and Upsilon < 1
#   (so Lemma mg:lem:upperE, which needs Upsilon >= 1, is indeed not available);  (6) Upsilon from Lemma mg:lem:Upsilon(i) = the variance computed directly;
#   (7) the dressed tube: Ubar_(C) from Lemma mg:lem:dressed (resp. mg:lem:dressedC) with the exact data, gamma_g <= gamma_2, is <= Upsilon_(C) + Delta (bound mode),
#   and the bound-mode inputs dominate the exact ones (n_B <= nbar_B, (u_2/u)^4 e^b <= cbar_2);  (8) a_1 = |<f_beta psi>| = 2u (SU(2)), sqrt2 u, 0 (U(1));
#   (9) Upsilon_- > Upsilon_+ for U(1) (so Upsilon = Upsilon_-, Lemma mg:lem:Upsilon(iii)).
# Usage: python3 applicability_check.py
import os, sys
from mpmath import mp, mpf, besseli, exp, expm1, sqrt, log, quad, cos, sin, pi
mp.dps = 50
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..'))
for sub in ('upper-bound',): sys.path.insert(0, os.path.join(here, '..', '..', sub))
from rigorous_bounds import bound_group
from constants_upper import beta0

def su2(beta):
    I = [besseli(n, beta) for n in range(0, 40)]
    u = I[2]/I[1]; ur = [I[n]/I[1] for n in range(2, 39)]          # u_j for d_j = n = 2j+1 >= 2
    u2 = max(ur[1:])
    s6 = sum(n**2*(I[n]/I[1])**6 for n in range(2, 39)); s5 = sum(n**2*(I[n]/I[1])**5 for n in range(2, 39))
    kappa = 3*I[3]/I[1]
    Ups = 1 + kappa - 4*u**2
    # direct variance of chi = 2 cos t under (2/pi) sin^2 t e^{beta cos t} dt
    Z = quad(lambda t: (2/pi)*sin(t)**2*exp(beta*cos(t)), [0, pi])
    m1 = quad(lambda t: (2/pi)*sin(t)**2*2*cos(t)*exp(beta*cos(t)), [0, pi])/Z
    m2 = quad(lambda t: (2/pi)*sin(t)**2*4*cos(t)**2*exp(beta*cos(t)), [0, pi])/Z
    a1 = quad(lambda t: (2/pi)*sin(t)**2*2*cos(t)*(exp(beta*cos(t))/Z - 1), [0, pi])
    return dict(u=u, ustar=max(ur), u2=u2, s6=s6, s5=s5, Ups={'': Ups}, Ups_direct={'': m2 - m1**2}, a1={'': a1}, a1_claim={'': 2*u}, N=2, iota=1)

def u1(beta):
    I = [besseli(n, beta) for n in range(0, 40)]
    u = I[1]/I[0]; ur = [I[n]/I[0] for n in range(1, 39)]
    u2 = max(ur[1:])
    s6 = 2*sum((I[n]/I[0])**6 for n in range(1, 39)); s5 = 2*sum((I[n]/I[0])**5 for n in range(1, 39))
    alpha, gamma = mpf(0), I[2]/I[0]
    Up, Um = 1 + alpha + gamma - 2*u**2, 1 + alpha - gamma
    Z = quad(lambda t: exp(beta*cos(t))/(2*pi), [-pi, pi])
    def var(h):
        m1 = quad(lambda t: h(t)*exp(beta*cos(t))/(2*pi), [-pi, pi])/Z
        m2 = quad(lambda t: h(t)**2*exp(beta*cos(t))/(2*pi), [-pi, pi])/Z
        return m2 - m1**2
    ap = quad(lambda t: sqrt(2)*cos(t)*(exp(beta*cos(t))/Z - 1)/(2*pi), [-pi, pi]); am = quad(lambda t: sqrt(2)*sin(t)*(exp(beta*cos(t))/Z - 1)/(2*pi), [-pi, pi])
    return dict(u=u, ustar=max(ur), u2=u2, s6=s6, s5=s5, Ups={'+': Up, '-': Um}, Ups_direct={'+': var(lambda t: sqrt(2)*cos(t)), '-': var(lambda t: sqrt(2)*sin(t))},
                a1={'+': abs(ap), '-': abs(am)}, a1_claim={'+': sqrt(2)*u, '-': mpf(0)}, N=1, iota=0)

def dressed_exact(beta, g, Ups):
    """Ubar of Lemma mg:lem:dressed (dressedC) with the exact u, u_2, sigma_6, and gamma_g <= gamma_2"""
    N, iota = g['N'], g['iota']; u, u2, s6 = g['u'], g['u2'], g['s6']
    eps = expm1(2*beta); gam = eps + beta**2*exp(beta)
    nB = (u2/u)**2*eps + sqrt(2 - iota)*(N/u)*sqrt(u2**5*expm1(beta))
    c2 = (u2/u)**4*exp(beta)
    assert Ups > c2
    v1 = Ups + nB**2/(Ups - c2)
    assert v1 > u**4*(1 + gam)
    return v1 + 2*exp(2*beta)*s6/(v1 - u**4*(1 + gam)), nB, c2

if __name__ == "__main__":
    allok = True
    for name, fn, (N, iota) in (('SU(2)', su2, (2, 1)), ('U(1)', u1, (1, 0))):
        for beta in (beta0(3), beta0(4), mpf('1e-4'), mpf('1e-6')):
            g = fn(beta); g['N'], g['iota'] = N, iota
            bg = bound_group(N, iota, sharp_u2=False)(beta)
            eps = expm1(beta)
            s5bar = (2 - iota)*N**2*bg['u']**5 + (expm1(beta) - beta)**4*expm1(beta)
            chk = {}
            chk['(1) u_* = u'] = g['ustar'] <= g['u']
            chk['(2) u_low <= u <= ubar'] = bg['u_low'] <= g['u'] <= bg['u']
            chk['(3) u_2 <= e^b-1-b'] = g['u2'] <= bg['u2'] and bg['u2'] == expm1(beta) - beta
            chk['(4) sigma_6 <= s6bar, sigma_5 <= s5bar'] = g['s6'] <= bg['s6'] and g['s5'] <= s5bar
            for C, U in g['Ups'].items():
                chk[f'(5) e^-2b <= Ups{C} <= e^b, Ups{C} < 1'] = exp(-2*beta) <= U <= exp(beta) and U < 1
                chk[f'(6) Ups{C} (Lemma 10.1(i)) = direct variance'] = abs(U - g['Ups_direct'][C]) < mpf(10)**-30
                Ub, nB, c2 = dressed_exact(beta, g, U)
                chk[f'(7) Ups{C} <= Ubar{C} <= Ups{C} + Delta'] = U <= Ub <= U + bg['Delta']
                chk[f'(8) a_1{C} = claimed'] = abs(g['a1'][C] - g['a1_claim'][C]) < mpf(10)**-30
            nBbar = (bg['u2']/bg['u_low'])**2*expm1(2*beta) + sqrt(2 - iota)*(N/bg['u_low'])*sqrt(bg['u2']**5*expm1(beta))
            chk['(7) n_B <= nbar_B, c_2 <= cbar_2'] = nB <= nBbar and c2 <= (bg['u2']/bg['u_low'])**4*exp(beta)
            if name == 'U(1)': chk['(9) Ups_- > Ups_+'] = g['Ups']['-'] > g['Ups']['+']
            ok = all(chk.values()); allok = allok and ok
            Ubar_info = ', '.join(f"Ups{C} = {mp.nstr(U, 12)}" for C, U in g['Ups'].items())
            print(f"{name}, beta = {mp.nstr(beta, 6)}: all hold: {ok};  u/beta = {mp.nstr(g['u']/beta, 10)} in [{mp.nstr(bg['u_low']/beta, 10)}, {mp.nstr(bg['u']/beta, 10)}],"
                  f" u_2/beta^2 = {mp.nstr(g['u2']/beta**2, 6)} <= {mp.nstr(bg['u2']/beta**2, 6)}, sigma_6/s6bar = {mp.nstr(g['s6']/bg['s6'], 6)}, {Ubar_info},"
                  f" Delta/beta^4 = {mp.nstr(bg['Delta']/beta**4, 4)}")
            if not ok: print('   failed:', [k for k, v in chk.items() if not v])
    print("all hold:", allok)

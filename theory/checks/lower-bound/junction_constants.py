# Theorem mg:thm:junctions and Remark mg:rem:junctions of Part IV (Section 8): the constants of Definition mg:def:jconst, iterated to a
# solution, and the resulting delta_J and loss for SU(2), d = 3, 4, in infinite and finite volume. The bounds implemented are those proved:
#   S_4 <= b (cap) + 4(P_v-1) b^2 (two plaquettes) + 2*24(d-2)(2d-3)*C(18,2) b^3 (three plaquettes, structured)
#          + sum_{n>=4} 2 C(4P_v+T_v n, n) (n+1)(2d-3) b^n                (two fronts, 4-exit, rooted counting)
#          + sum_{n>=1} 2 C(4P_v+T_v n, n) b^n sum_{k>=k0(n)} C(4n+4,k) f(k)   (two fronts, big exit)
#          + anchored doubly-marked tree sum                             (three or more fronts)
#   S_k (k>=5, big entry) <= k P_v [gC Y' + gF k Y] e^{eta k}            (anchored singly-marked tree sum)
#   f(4) = 1, f(k) = F k^2 e^{eta k} a^{k-4} (k>=5);  conditions S_4 < delta', S_k <= F k^2 e^{eta k} delta'.
#   decorations: h_T = sup_k a^{k/2} k P_v Y e^{eta k},  c_T = e^{h_T}/(1-a^4).
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'upper-bound'))
from mpmath import mp, mpf, e, exp, log, besseli, expm1, binomial
from constants_upper import consts, beta0
mp.dps = 30

def u_su2(b): return besseli(2, b) / besseli(1, b)

def counts(d, volume):
    n = d - 1; Pv = 2*n*(n-1)
    N = {4: Pv + (n if volume == 'finite' else 0), 5: (n if volume == 'finite' else 0)}
    N[6] = 16*n**5                                   # 6-edge fronts are 6-cycles: closed walks of length 6 / 2
    for k in range(7, 200): N[k] = n * binomial((4*d-6)*k, k-1)
    return Pv, 4*(2*d-5), 4*(Pv-1), N

def junction_bound(beta, d, u, volume='infinite', KMAX=60, NMAX=40, verbose=False, gamma='gamma2'):
    K = 4*(2*d-3)
    eps = expm1(2*beta); x = e**2*K*eps; th5 = x**5/(1-x); lam = exp(2*(K+1)*th5)
    a = lam*u; b = lam*((expm1(2*beta) + beta**2*exp(beta)) if gamma == 'gamma2' else expm1(4*beta))
    Pv, KS, Tv, N = counts(d, volume)
    def Yf(z, p=0): return sum(m**p*binomial(KS*m, m-1)*(b*exp(z))**m for m in range(1, NMAX))
    # junction tree-sum constants eta, zeta (front weight a^{k/2} e^{cT}), iterate with c_T
    cT = mpf(1); eta = mpf(0); zeta = mpf(0)
    for _ in range(200):
        A = exp(cT)*sum(N[k]*a**(mpf(k)/2)*exp(eta*k) for k in range(4, KMAX))
        eta_n = Pv*Yf(zeta) + A
        zeta_n = 4*Pv*Yf(zeta) + 8*A
        hT = max(a**(mpf(k)/2)*k*Pv*Yf(zeta)*exp(eta*k) for k in range(4, 12))
        cT_n = exp(hT)/(1-a**4)
        if abs(eta_n-eta) + abs(zeta_n-zeta) + abs(cT_n-cT) < mpf(10)**-25: eta, zeta, cT = eta_n, zeta_n, cT_n; break
        eta, zeta, cT = eta_n, zeta_n, cT_n
    A = exp(cT)*sum(N[k]*a**(mpf(k)/2)*exp(eta*k) for k in range(4, KMAX))
    A1 = exp(cT)*sum(k*N[k]*a**(mpf(k)/2)*exp(eta*k) for k in range(4, KMAX))
    Y, Y1, Y2 = Yf(zeta), Yf(zeta, 1), Yf(zeta, 2)
    def solve(F):
        f = lambda k: mpf(1) if k == 4 else F*k**2*exp(eta*k)*a**(k-4)
        B1 = sum(N[k]*f(k)*exp(eta*k) for k in range(4, KMAX))          # exits through a vertex, one layer
        # singly marked: gF(1-A1) = Pv gC Y1 + B1 ; gC(1-4Pv Y1) = 8 B1 + 8 gF A1
        M11, M12, r1 = 1-A1, -Pv*Y1, B1
        M21, M22, r2 = -8*A1, 1-4*Pv*Y1, 8*B1
        det = M11*M22-M12*M21; gF = (r1*M22-M12*r2)/det; gC = (M11*r2-M21*r1)/det
        # extra-front marks (weight a^{k/2}e^cT): same system with B1 -> A
        det2 = det; g2F = (A*M22-M12*8*A)/det2; g2C = (M11*8*A-M21*A)/det2
        A2 = exp(cT)*sum(k**2*N[k]*a**(mpf(k)/2)*exp(eta*k) for k in range(4, KMAX))
        Bf1 = sum(k*N[k]*f(k)*exp(eta*k) for k in range(4, KMAX))
        # kappa_C = gC g2C + 8(g2F Bf1 + gF A1 + kF A2) + 4 Pv kC Y2 ;  kappa_F = gF g2F + g2F Bf1 + gF A1 + kF A2 + Pv kC Y2
        a11, a12, r1 = 1-4*Pv*Y2, -8*A2, gC*g2C + 8*(g2F*Bf1 + gF*A1)
        a21, a22, r2 = -Pv*Y2, 1-A2, gF*g2F + g2F*Bf1 + gF*A1
        dt = a11*a22-a12*a21; kC = (r1*a22-a12*r2)/dt; kF = (a11*r2-a21*r1)/dt
        # anchored bounds at a 4-entry
        S3plus = 4*Pv*exp(4*eta)*(kC*Y2 + 4*(gC*g2F + g2C*gF)*Y1 + 16*kF*Y)
        # two-front parts
        S_cap = b; S_2 = (4*(Pv-1)-KS)*b**2; S_3 = 2*24*(d-2)*(2*d-3)*binomial(18, 2)*b**3
        S_4plus = sum(2*binomial(4*Pv+Tv*n, n)*(n+1)*(2*d-3)*b**n for n in range(4, NMAX))
        def nfr(n, k):
            if k == 5: return mpf(4*n+4)/5 if (N[5] > 0 and n >= 4) else 0
            return min(binomial(4*n+4, k), (4*n+4)*N[k])
        S_big = sum(2*binomial(4*Pv+Tv*n, n)*b**n*sum(nfr(n, k)*f(k) for k in range(5, min(4*n+5, 199))) for n in range(1, NMAX))
        S4 = S_cap + S_2 + S_3 + S_4plus + S_big + S3plus
        # big entries: S_k <= k Pv [gC Y1 + gF k Y] e^{eta k}; need <= F k^2 e^{eta k} delta'
        Fneed = lambda dl: max(Pv*(gC*Y1 + gF*k*Y)/(k*dl) for k in range(5, KMAX))
        return dict(S4=S4, parts=(S_cap, S_2, S_3, S_4plus, S_big, S3plus), Fneed=Fneed, gF=gF, gC=gC, B1=B1)
    F = mpf(1)
    for _ in range(100):
        r = solve(F); dl = r['S4']*(1+mpf(10)**-12)
        Fn = r['Fneed'](dl)
        if abs(Fn-F) < mpf(10)**-20*F: F = Fn; break
        F = Fn
    r = solve(F); delta = r['S4']
    return dict(delta=delta, b=b, a=a, F=F, eta=eta, zeta=zeta, cT=cT, parts=r['parts'], th5=th5, K=K, gF=r['gF'], gC=r['gC'])

def worst_case():
    """Hypotheses at beta_0(d) for every G and rho: u_* <= eps, ||g|| <= e^{4beta}-1 (or gamma_2 without trivial summand);
    counts for spatial sides >= 6 (Remark mg:rem:junctions(i))."""
    print("worst case over G, rho at beta_0(d) (u_* = eps), counts for spatial sides >= 6:")
    for gam, label in (('e4', 'any rho'), ('gamma2', 'rho without trivial summand')):
        for d in (3, 4, 5, 6):
            b0 = beta0(d)*(1-mpf(10)**-12)
            r = junction_bound(b0, d, expm1(2*b0), 'infinite', gamma=gam)
            loss = 8*(r['K']+1)*r['th5'] + log(1+r['delta'])
            print(f"  {label:28s} d={d}: delta_J={float(r['delta']):.3e} (delta_J/b={float(r['delta']/r['b']):.3f}), loss {float(loss):.3f}")

if __name__ == "__main__":
    worst_case()
    for d in (3, 4):
        b0 = beta0(d)
        for volume in ('infinite', 'finite'):
            print(f"d={d}, {volume} volume:")
            for beta in (mpf('1e-7'), mpf('1e-6'), mpf('1e-5'), b0/2, b0*(1-mpf(10)**-12)):
                r = junction_bound(beta, d, u_su2(beta), volume)
                loss = 8*(r['K']+1)*r['th5'] + log(1+r['delta'])
                D = consts(beta, d, 2, 1)['D']
                pc = ', '.join(f'{float(x/r["b"]):.2e}' for x in r['parts'])
                print(f"  beta={float(beta):.4e}: delta={float(r['delta']):.4e} (delta/b={float(r['delta']/r['b']):.5f}; parts/b: {pc}); "
                      f"F={float(r['F']):.3g}, eta={float(r['eta']):.2e}, gF={float(r['gF']):.3g}, gC={float(r['gC']):.3g}; loss {float(loss):.3e}; window {float(loss+D):.3e}")
                sys.stdout.flush()

# Written by the second independent referee of Section 15 of Part IV (2026-11-12) and added to the repository with credit; unchanged
# apart from this header. A compass, not part of the proofs.
"""Referee check of Lemma mg:lem:dressedC (the dressed tube in a sector of charge conjugation), and of the U(1) numbers
of Conjecture mg:conj:codd. A compass, not part of the proofs.

(A) For G = Z_n (n = 3, 4, 5, 7), rho the character k -> omega^k (iota = 0), charge conjugation c(U) = U^{-1}:
    the dressed tube P^q D_q P^q on the class functions of the plaquette holonomy in the two copies is
        T = Lambda (m (x) m) Lambda,   Lambda e_{aa'} = uhat_{aa'}^2 e_{aa'},   m_{ab} = <chi_a, (1+f) chi_b>
    (proof of Lemma mg:lem:dressed). We check that C (e_{aa'} -> e_{-a,-a'}) commutes with T and with the swap S,
    compute the top eigenvalue lambda_C of T on the range of P_- Q_C, and compare it with u^4 Upsilon_C and with
    u^4 bar-Upsilon_C, bar-Upsilon_C given by the formulas of Lemma mg:lem:dressed with Upsilon_C in place of Upsilon.
    We also check that the compression of T to (e^-_rho + C e^-_rhobar)/sqrt2 is u^4 Upsilon_C exactly.
(B) U(1): Upsilon_- = 1 - I_2/I_0, u = I_1/I_0; the series of -log Upsilon_- in u, and the coefficients claimed in
    Conjecture mg:conj:codd: -425/24 (d = 4), -185/24 (d = 3), and h_- u_1^4 = (10d-22)/16; and h_+ - h_- = 16(d-3).
Run: python3 sector_dressed_tube.py   (seconds)
"""
import numpy as np
import sympy as sp


def zn_data(n, beta):
    k = np.arange(n)
    w = np.exp(beta * np.cos(2 * np.pi * k / n))       # Wilson weight, N = d_rho = 1
    chi = lambda a: np.exp(2j * np.pi * a * k / n)
    c0 = w.mean()
    f = w / c0 - 1
    u = np.array([(f * np.conj(chi(a))).mean().real for a in range(n)])   # u_a = int f conj(chi_a), d_a = 1
    u[0] = 0.0
    m = np.array([[(np.conj(chi(a)) * (1 + f) * chi(b)).mean() for b in range(n)] for a in range(n)])
    return k, f, u, m, chi


def check_zn(n, beta):
    k, f, u, m, chi = zn_data(n, beta)
    uh = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            if a == 0 and b == 0:
                uh[a, b] = 0.0
            elif a == 0:
                uh[a, b] = u[b]
            elif b == 0:
                uh[a, b] = u[a]
            else:
                uh[a, b] = u[a] * u[b]
    L = np.diag((uh ** 2).ravel())
    T = L @ np.kron(m, m) @ L
    idx = lambda a, b: a * n + b
    S = np.zeros((n * n, n * n)); Cc = np.zeros((n * n, n * n))
    for a in range(n):
        for b in range(n):
            S[idx(b, a), idx(a, b)] = 1
            Cc[idx((-a) % n, (-b) % n), idx(a, b)] = 1
    comm = max(np.abs(T @ Cc - Cc @ T).max(), np.abs(T @ S - S @ T).max(), np.abs(S @ Cc - Cc @ S).max())
    I = np.eye(n * n)
    Pm = (I - S) / 2
    usup = u[1:].max()
    uu = u[1]                      # u = u_rho, rho = 1
    # V on span(chi_1, chi_{-1}) (Lemma mg:lem:Upsilon(i)); Upsilon_+ for sqrt2 Re chi, Upsilon_- for sqrt2 Im chi
    def var(psi):
        return ((psi ** 2) * (1 + f)).mean() - ((psi * (1 + f)).mean()) ** 2
    Ups = {+1: var(np.sqrt(2) * np.cos(2 * np.pi * k / n)), -1: var(np.sqrt(2) * np.sin(2 * np.pi * k / n))}
    # constants of Lemma mg:lem:dressed
    others = [u[a] for a in range(1, n) if a not in (1, n - 1)]
    u2 = max(others) if others else 0.0
    eps = np.exp(2 * beta) - 1
    sig6 = (u[1:] ** 6).sum()
    ff = 1 + f
    gam_g = np.abs(np.outer(ff, ff) - 1).max()
    nB = (u2 ** 2 / uu ** 2) * eps + np.sqrt(2) * (1 / uu) * np.sqrt(u2 ** 5 * (np.exp(beta) - 1))
    res = {}
    for C in (+1, -1):
        Q = (I + C * Cc) / 2
        Pr = Pm @ Q
        A = Pr @ T @ Pr
        lam = np.linalg.eigvalsh((A + A.T.conj()) / 2).max()
        v = np.zeros(n * n, complex)
        # e^-_a = (e_{a0} - e_{0a})/sqrt2 ;  (e^-_1 + C e^-_{-1})/sqrt2
        for a, coef in ((1, 1.0), (n - 1, C)):
            v[idx(a, 0)] += coef / 2
            v[idx(0, a)] -= coef / 2
        comp = (v.conj() @ T @ v).real
        Ub = Ups[C]
        ups1 = Ub + nB ** 2 / (Ub - (u2 / uu) ** 4 * np.exp(beta))
        Ubar = ups1 + 2 * np.exp(2 * beta) * sig6 / (ups1 - uu ** 4 * (1 + gam_g))
        res[C] = (lam, uu ** 4 * Ub, comp, uu ** 4 * Ubar)
    lam_all = np.linalg.eigvalsh(Pm @ T @ Pm).max()
    return comm, usup <= uu + 1e-15, res, lam_all, Ups


def u1_values():
    b, x = sp.symbols('beta x', positive=True)
    order = 10
    I = lambda nu: sp.series(sp.besseli(nu, b), b, 0, order + 2).removeO()
    I0, I1, I2 = I(0), I(1), I(2)
    u = sp.series(I1 / I0, b, 0, order).removeO()
    Um = sp.series(1 - I2 / I0, b, 0, order).removeO()
    Up = sp.series(1 + I2 / I0 - 2 * (I1 / I0) ** 2, b, 0, order).removeO()
    # invert u(beta) = x
    binv = sp.symbols('binv')
    c = sp.symbols('c1:8')
    ansatz = 2 * x + c[0] * x ** 3 + c[1] * x ** 5 + c[2] * x ** 7
    eq = sp.series(u.subs(b, ansatz) - x, x, 0, 8).removeO()
    sol = sp.solve([eq.coeff(x, j) for j in (3, 5, 7)], c[:3], dict=True)[0]
    bx = ansatz.subs(sol)
    Lm = sp.series(-sp.log(Um.subs(b, bx)), x, 0, 7).removeO()
    Lp = sp.series(-sp.log(Up.subs(b, bx)), x, 0, 7).removeO()
    return sp.expand(Lm), sp.expand(Lp), sp.series(Um - 2 * u / b, b, 0, 8)


if __name__ == '__main__':
    print('(A) Lemma mg:lem:dressedC for Z_n (two copies, class functions of the plaquette holonomy):')
    for n in (3, 4, 5, 7):
        for beta in (0.3, 0.1, 0.03):
            comm, ustar, res, lam_all, Ups = check_zn(n, beta)
            line = '  Z_%d beta=%.2f: [T,C],[T,S],[S,C] max %.1e; u_*=u: %s;' % (n, beta, comm, ustar)
            for C in (+1, -1):
                lam, lo, comp, hi = res[C]
                ok = lo * (1 - 1e-12) <= lam <= hi * (1 + 1e-12)
                line += ' C=%+d: lam/(u^4 Ups_C)-1 = %.2e, compression/(u^4 Ups_C)-1 = %.1e, (Ubar_C/Ups_C-1) = %.2e, in [u^4Ups_C, u^4Ubar_C]: %s;' % (
                    C, lam / lo - 1, comp / lo - 1, hi / lo - 1, ok)
                assert ok and abs(comp / lo - 1) < 1e-10
            line += ' odd-sector top (C-blind)/max_C lam_C - 1 = %.1e' % (lam_all / max(res[1][0], res[-1][0]) - 1)
            print(line)
    print()
    print('(B) U(1) and the hopping numbers:')
    Lm, Lp, chk = u1_values()
    print('  -log Upsilon_- = %s + O(u^7)    [x = u]' % Lm)
    print('  -log Upsilon_+ = %s + O(u^7)' % Lp)
    print('  Upsilon_- - 2u/beta =', chk)
    c4 = Lm.coeff(sp.symbols('x', positive=True), 4)
    for d in (3, 4, 5):
        hm, hp = 10 * d - 22, 26 * d - 70
        print('  d=%d: h_- = %d, h_+ = %d, h_+ - h_- = %d = 16(d-3) = %d; jumps 2(d-3) + coplanar bulges/hops 8 + perpendicular bulges 8(d-3) = %d;'
              ' u^4 coefficient of m+4log u (C=-) = %s; h_- u_1^4 = %s'
              % (d, hm, hp, hp - hm, 16 * (d - 3), 2 * (d - 3) + 8 + 8 * (d - 3), c4 - hm, sp.Rational(hm, 16)))

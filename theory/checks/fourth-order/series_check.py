# series_check.py -- written by the independent referee of Section 12 of Part IV (2026-11-12), added to the repository with credit;
# unchanged apart from this header. A compass, not part of the proofs.
# Referee check (independent of theory/checks): exact rational power series of u, Upsilon_+-, and the
# predicted fourth-order glueball masses, from moments <chi^j chibar^k> computed by Weyl integration
# (constant term of chi^j chibar^k |Delta|^2 on the maximal torus), for U(1), SU(2), SU(3).
#   w = exp((beta/N) Re chi),  c0 = <w>,  u = <chibar w>/(N c0),
#   Upsilon_+ = Var_w(sqrt2 Re chi) (iota = 0), Upsilon = Var_w(chi) (SU(2)), Upsilon_- = Var_w(sqrt2 Im chi).
# Also: LMP07's abbreviations a, b, c (eq. 2.9) are compared with Upsilon, Upsilon_-, Upsilon_+.
from fractions import Fraction as Fr
from math import comb, factorial
import sys

D = 14  # maximal total degree j+k of the moments (series in beta to order ~D-3)

# ---------- Laurent polynomials as dicts {exponent tuple: int} ----------
def pmul(a, b):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            out[e] = out.get(e, 0) + ca * cb
    return {e: c for e, c in out.items() if c}

def ppow(a, n, one):
    r = dict(one)
    for _ in range(n):
        r = pmul(r, a)
    return r

def moments(group, D):
    """M[j][k] = <chi^j chibar^k> (exact Fractions), j + k <= D."""
    if group == 'U1':
        chi = {(1,): 1}; chib = {(-1,): 1}; one = {(0,): 1}; meas = {(0,): 1}; norm = 1
    elif group == 'SU2':
        chi = {(1,): 1, (-1,): 1}; chib = dict(chi); one = {(0,): 1}
        meas = {(0,): 2, (2,): -1, (-2,): -1}; norm = 2          # |z - 1/z|^2, Weyl factor 1/2
    elif group == 'SU3':
        # z1 = (1,0), z2 = (0,1), z3 = (-1,-1)
        z = [(1, 0), (0, 1), (-1, -1)]
        chi = {e: 1 for e in z}; chib = {(-a, -b): 1 for (a, b) in z}; one = {(0, 0): 1}
        Delta = dict(one)
        for i in range(3):
            for j in range(i + 1, 3):
                Delta = pmul(Delta, {z[i]: 1, z[j]: -1})
        Deltab = {(-a, -b): c for (a, b), c in Delta.items()}
        meas = pmul(Delta, Deltab); norm = 6
    else:
        raise ValueError
    zero = tuple(0 for _ in next(iter(one)))
    chipow = [dict(one)]
    for _ in range(D): chipow.append(pmul(chipow[-1], chi))
    chibpow = [dict(one)]
    for _ in range(D): chibpow.append(pmul(chibpow[-1], chib))
    M = [[None] * (D + 1) for _ in range(D + 1)]
    for j in range(D + 1):
        A = pmul(chipow[j], meas)
        for k in range(D + 1 - j):
            B = chibpow[k]
            ct = 0
            for e, c in B.items():
                ne = tuple(-x for x in e)
                if ne in A: ct += A[ne] * c
            M[j][k] = Fr(ct, norm)
    return M

# ---------- truncated power series in beta ----------
def smul(p, q, n):
    r = [Fr(0)] * (n + 1)
    for i, a in enumerate(p[:n + 1]):
        if a:
            for j, b in enumerate(q[:n + 1 - i]):
                r[i + j] += a * b
    return r

def sinv(p, n):
    r = [Fr(0)] * (n + 1); r[0] = 1 / p[0]
    for k in range(1, n + 1):
        r[k] = -sum(p[i] * r[k - i] for i in range(1, k + 1) if i < len(p)) / p[0]
    return r

def slog(p, n):  # p[0] = 1
    assert p[0] == 1
    dp = [i * p[i] for i in range(1, len(p))]
    q = smul(dp, sinv(p, n), n)
    return [Fr(0)] + [q[k - 1] / k for k in range(1, n + 1)]

def scompose(p, q, n):  # p(q(x)), q[0] = 0
    out = [Fr(0)] * (n + 1); pw = [Fr(1)] + [Fr(0)] * n
    for c in p[:n + 1]:
        out = [o + c * x for o, x in zip(out, pw)]
        pw = smul(pw, q, n)
    return out

def srevert(q, n):  # inverse of q = q1 x + ...
    r = [Fr(0), 1 / q[1]] + [Fr(0)] * (n - 1)
    for k in range(2, n + 1):
        r[k] = -scompose(q, r, k)[k] / q[1]
    return r

def expect_series(M, N, fun, n):
    """series in beta of <F exp((beta/N) Re chi)> where F = sum_{(a,b)} coef chi^a chibar^b (dict fun)."""
    # Re chi^m = 2^-m sum_j C(m,j) chi^j chibar^(m-j)
    out = []
    for m in range(n + 1):
        s = Fr(0)
        for j in range(m + 1):
            for (a, b), c in fun.items():
                s += Fr(comb(m, j), 2 ** m) * c * M[j + a][m - j + b]
        out.append(s / factorial(m) / Fr(N) ** m)
    return out

def analyse(group, N, iota, n, hs):
    M = moments(group, D)
    c0 = expect_series(M, N, {(0, 0): 1}, n)
    cb = expect_series(M, N, {(0, 1): 1}, n)              # <chibar w>
    u = smul(cb, sinv(c0, n), n); u = [x / N for x in u]
    # Re chi = (chi + chib)/2 ; Im chi = (chi - chib)/(2i)
    ReW = expect_series(M, N, {(1, 0): Fr(1, 2), (0, 1): Fr(1, 2)}, n)
    Re2W = expect_series(M, N, {(2, 0): Fr(1, 4), (1, 1): Fr(1, 2), (0, 2): Fr(1, 4)}, n)
    Im2W = expect_series(M, N, {(2, 0): Fr(-1, 4), (1, 1): Fr(1, 2), (0, 2): Fr(-1, 4)}, n)
    ImW = expect_series(M, N, {(1, 0): Fr(1, 2), (0, 1): Fr(-1, 2)}, n)  # = 2i <Im chi w>, must vanish
    assert all(x == 0 for x in ImW), "Im chi has nonzero mean?"
    ic0 = sinv(c0, n)
    mR = smul(ReW, ic0, n)
    if iota == 1:   # SU(2): chi real, Upsilon = Var_w(chi) = Var_w(2 Re chi)/... chi = 2 Re chi? no: chi real => chi = Re chi
        Ups = {'+': [x - y for x, y in zip(smul(Re2W, ic0, n), smul(mR, mR, n))]}
    else:
        Up = [2 * (x - y) for x, y in zip(smul(Re2W, ic0, n), smul(mR, mR, n))]
        Um = [2 * x for x in smul(Im2W, ic0, n)]
        Ups = {'+': Up, '-': Um}
    res = {'u': u, 'Ups': Ups, 'M': M}
    # series in u
    k = n - 1
    binv = srevert(u[:k + 1], k)
    res['beta_of_u'] = binv
    res['mlogUps_u'] = {}
    for C, U in Ups.items():
        ml = [-x for x in slog(U[:k + 1], k)]
        res['mlogUps_u'][C] = scompose(ml, binv, k)
        res.setdefault('mlogUps_b', {})[C] = ml
    return res

def fmt(p, var, kmax):
    return ' + '.join(f'({c}){var}^{i}' for i, c in enumerate(p[:kmax + 1]) if c and i)

if __name__ == '__main__':
    n = 8
    out = {}
    for group, N, iota in (('U1', 1, 0), ('SU2', 2, 1), ('SU3', 3, 0)):
        r = analyse(group, N, iota, n, None)
        out[group] = r
        print(f'== {group}: u(beta) = {fmt(r["u"], "b", 7)}')
        for C, s in r['mlogUps_u'].items():
            print(f'   -log Upsilon_{C} in u : {fmt(s, "u", 6)}')
            print(f'   -log Upsilon_{C} in b : {fmt(r["mlogUps_b"][C], "b", 6)}')
    # moments sanity: SU(3) <chi^4 chibar> = 3, <chi^3> = 1, <|chi|^4> = 2
    M3 = out['SU3']['M']
    print('SU(3) moments: <chi^3> =', M3[3][0], ' <|chi|^4> =', M3[2][2], ' <chi^4 chibar> =', M3[4][1], ' <|chi|^6> =', M3[3][3])
    M2 = out['SU2']['M']
    print('SU(2) moments: <chi^2> =', M2[2][0], ' <chi^4> =', M2[4][0], ' <chi^6> =', M2[6][0])
    # Predictions: m + 4 log u = -log Upsilon_C - h u^4
    print('\nPredicted m + 4 log u (through u^4):')
    for d in (3, 4):
        hp, hm = 26 * d - 70, 10 * d - 22
        for group, C, h, label in (('SU2', '+', hp, 'SU(2) A1++'), ('SU2', '+', 10 if d == 4 else None, 'SU(2) E++'),
                                    ('SU3', '+', hp, 'SU(3) A1++'), ('SU3', '+', 10 if d == 4 else None, 'SU(3) E++'),
                                    ('SU3', '-', hm, 'SU(3) T1+-'), ('U1', '-', hm, 'U(1) C-'), ('U1', '+', hp, 'U(1) C+')):
            if h is None: continue
            s = list(out[group]['mlogUps_u'][C][:5]); s[4] -= h
            print(f'  d={d} {label:11s} (h={h:2d}): ' + fmt(s, 'u', 4))
    # LMP07 eq. (2.9): a = 1+3v-4u^2 (SU(2), v = u_{adj}); b = 1-3u-6v+8w, c = 1+3u+6v+8w-18u^2 (SU(3), v = u_6, w = u_8)
    print('\nLMP07 (2.9) abbreviations versus Upsilon (exact series in beta):')
    r = out['SU2']; M = r['M']; N = 2
    c0 = expect_series(M, N, {(0, 0): 1}, n)
    # SU(2): chi_1 (adjoint) = chi^2 - 1 ; u_1 = <chi_1 w>/(3 c0)
    v = [x / 3 for x in smul(expect_series(M, N, {(2, 0): 1, (0, 0): -1}, n), sinv(c0, n), n)]
    a = [Fr(1)] + [Fr(0)] * n
    a = [x + 3 * y - 4 * z for x, y, z in zip(a, v, smul(r['u'], r['u'], n))]
    print('  SU(2): a - Upsilon =', [x - y for x, y in zip(a, r['Ups']['+'])][:n + 1])
    r = out['SU3']; M = r['M']; N = 3
    c0 = expect_series(M, N, {(0, 0): 1}, n); ic0 = sinv(c0, n)
    # chi_6 = chi^2 - chibar, so chibar_6 = chibar^2 - chi ; chi_8 = chi chibar - 1
    v = [x / 6 for x in smul(expect_series(M, N, {(0, 2): 1, (1, 0): -1}, n), ic0, n)]
    w = [x / 8 for x in smul(expect_series(M, N, {(1, 1): 1, (0, 0): -1}, n), ic0, n)]
    u = r['u']; one = [Fr(1)] + [Fr(0)] * n
    b = [o - 3 * x - 6 * y + 8 * z for o, x, y, z in zip(one, u, v, w)]
    c = [o + 3 * x + 6 * y + 8 * z - 18 * q for o, x, y, z, q in zip(one, u, v, w, smul(u, u, n))]
    print('  SU(3): v = u_6 =', fmt(v, 'b', 4), ';  w = u_8 =', fmt(w, 'b', 4))
    print('  SU(3): b - Upsilon_- =', [x - y for x, y in zip(b, r['Ups']['-'])][:n + 1])
    print('  SU(3): c - Upsilon_+ =', [x - y for x, y in zip(c, r['Ups']['+'])][:n + 1])
    # v, w as series in u (for the discussion of the -98/4 vs -99/4 discrepancy)
    k = n - 1; binv = r['beta_of_u']
    print('  SU(3): v(u) =', fmt(scompose(v[:k + 1], binv, k), 'u', 5))
    print('  SU(3): w(u) =', fmt(scompose(w[:k + 1], binv, k), 'u', 5))

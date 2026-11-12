# Written by the second independent referee of Section 15 of Part IV (2026-11-12) and added to the repository with credit; unchanged
# apart from this header. A compass, not part of the proofs.
"""Referee check of Lemma mg:lem:split (Part IV, Section 15), independent of upsilon_split.py.

No representation theory is assumed.  The generating function of SU(N),
    Z_N(alpha, gamma) := < exp(alpha chi + gamma conj(chi)) >_{SU(N)},   chi = Tr U,
is computed from the Heine/Toeplitz identity for U(N) summed over the powers of det U (which turns
the U(N) integral into the SU(N) one):
    Z_N = sum_{q in Z} det[ fhat_{j-k+q} ]_{j,k=1..N},
    fhat_n = sum_m alpha^(m+n) gamma^m / (m!(m+n)!)  (n >= 0),   fhat_{-n} = same with alpha <-> gamma,
the Fourier coefficients of exp(alpha z + gamma/z).  The determinants are expanded over permutations with
pruning by the minimal degree |j-k+q| of each entry, in exact rational arithmetic, up to total degree DEG.
From Z_N we read the moments <chi^a conj(chi)^b> = a! b! [alpha^a gamma^b] Z_N, and then the exact power
series in beta of
    Upsilon_+ = Var_y(sqrt2 Re chi),  Upsilon_- = Var_y(sqrt2 Im chi),  y = (beta/N) Re chi,
where Var_y is the variance under the density e^y/<e^y>.
Checks:
  (1) the moments used in the proof: zero unless a = b (mod N); <|chi|^{2k}> = k! for k <= 3; <chi^N> = 1;
  (2) Upsilon_+ - Upsilon_- = 2(<conj(chi)^2>_y - <conj(chi)>_y^2) (identity of the proof);
  (3) Upsilon_+ - Upsilon_- = beta^4/248832 + O(beta^5) for N = 6, O(beta^5) for N >= 7,
      and the remark after the lemma: beta^5/32269440 for N = 7, O(beta^6) for N = 8,
      against 2 (beta/2N)^(N-2)/(N-2)!;
  (4) Lemma mg:lem:Upsilon(iii) and mg:lem:Upsilon3 for N = 3, 4, 5 (the values used in Corollary mg:cor:sectors(ii));
  (5) the Gaussian cancellation: with the moments delta_ab a! the expression vanishes to the order computed.
Run: python3 split_toeplitz_exact.py   (about a minute)
"""
from fractions import Fraction as Fr
from math import factorial, comb
import sys

DEG = int(sys.argv[1]) if len(sys.argv) > 1 else 10  # total degree in (alpha, gamma) kept


def padd(p, q, c=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + c * v
        if r[k] == 0:
            del r[k]
    return r


def pmul(p, q, deg=DEG):
    r = {}
    for (a1, b1), v1 in p.items():
        for (a2, b2), v2 in q.items():
            if a1 + a2 + b1 + b2 <= deg:
                k = (a1 + a2, b1 + b2)
                r[k] = r.get(k, 0) + v1 * v2
    return {k: v for k, v in r.items() if v != 0}


def fhat(n, deg=DEG):
    """Fourier coefficient z^n of exp(alpha z + gamma / z), truncated at total degree deg."""
    p = {}
    m = 0
    while True:
        a, b = (m + n, m) if n >= 0 else (m, m - n)
        if a + b > deg:
            break
        p[(a, b)] = Fr(1, factorial(a) * factorial(b))
        m += 1
    return p


def det_series(N, q, deg=DEG):
    """det[fhat_{j-k+q}] by expansion over permutations, pruned by the minimal degree of the entries."""
    F = {n: fhat(n, deg) for n in range(-(N + abs(q) + 2), N + abs(q) + 3)}
    total = {}

    def rec(j, used, sign, prod, mindeg):
        nonlocal total
        if j == N:
            total = padd(total, prod, sign)
            return
        # columns not used; sign of the permutation by counting inversions incrementally
        for k in range(N):
            if used >> k & 1:
                continue
            n = j - k + q
            dmin = abs(n)
            if mindeg + dmin > deg:
                continue
            inv = bin(used >> (k + 1)).count('1')  # number of used columns > k  -> inversions
            newprod = pmul(prod, F[n], deg)
            if not newprod:
                continue
            rec(j + 1, used | (1 << k), sign * (-1) ** inv, newprod, mindeg + dmin)

    rec(0, 0, 1, {(0, 0): Fr(1)}, 0)
    return total


def Z_series(N, deg=DEG):
    Z = {}
    qmax = deg // N
    for q in range(-qmax, qmax + 1):
        Z = padd(Z, det_series(N, q, deg))
    return Z


def moments(N, deg=DEG):
    Z = Z_series(N, deg)
    return {(a, b): v * factorial(a) * factorial(b) for (a, b), v in Z.items()}


def ups_series(mom, N, order):
    """Exact series (list of Fractions, index = power of beta) of Upsilon_+, Upsilon_-, and E := <cb^2>_y - <cb>_y^2.
    y = (beta/N) Re chi = t (chi + cb), t = beta/(2N)."""
    def M(a, b):
        return mom.get((a, b), Fr(0))

    def avg_cbk_ey(k, c=0):
        # < chi^c cb^k e^y > as a series in beta up to 'order'
        s = [Fr(0)] * (order + 1)
        for n in range(order + 1):
            # (beta/2N)^n / n! * sum_j C(n,j) <chi^{j+c} cb^{n-j+k}>
            tot = Fr(0)
            for j in range(n + 1):
                a, b = j + c, n - j + k
                if a + b > DEG:
                    raise ValueError('degree too high: need DEG >= %d' % (a + b))
                tot += comb(n, j) * M(a, b)
            s[n] = tot / factorial(n) / Fr((2 * N) ** n)
        return s

    def smul(p, q):
        r = [Fr(0)] * (order + 1)
        for i, x in enumerate(p):
            if x == 0:
                continue
            for j, y in enumerate(q):
                if i + j > order:
                    break
                r[i + j] += x * y
        return r

    def sinv(p):
        r = [Fr(0)] * (order + 1)
        r[0] = 1 / p[0]
        for n in range(1, order + 1):
            r[n] = -sum(p[k] * r[n - k] for k in range(1, n + 1)) / p[0]
        return r

    Zs = avg_cbk_ey(0)
    iZ = sinv(Zs)
    cb1 = smul(avg_cbk_ey(1), iZ)      # <cb>_y  (real)
    cb2 = smul(avg_cbk_ey(2), iZ)      # <cb^2>_y
    abs2 = smul(avg_cbk_ey(1, 1), iZ)  # <|chi|^2>_y
    Vdiag = [x - y for x, y in zip(abs2, smul(cb1, cb1))]
    Voff = [x - y for x, y in zip(cb2, smul(cb1, cb1))]
    Up = [x + y for x, y in zip(Vdiag, Voff)]
    Um = [x - y for x, y in zip(Vdiag, Voff)]
    return Up, Um, Voff


def ups_direct(mom, N, order):
    """Independent route to Upsilon_+ - Upsilon_-: Var(sqrt2 Re chi) - Var(sqrt2 Im chi) from the moments directly,
    with Re chi = (chi+cb)/2, Im chi = (chi-cb)/(2i), so 2(Re^2 - Im^2) = chi^2 + cb^2 and <Im chi>_y = 0."""
    Up, Um, Voff = ups_series(mom, N, order)
    return [x - y for x, y in zip(Up, Um)], [2 * v for v in Voff]


def fmt(series, upto=None):
    terms = []
    for n, c in enumerate(series):
        if upto is not None and n > upto:
            break
        if c != 0:
            terms.append('%s*b^%d' % (c, n))
    return ' + '.join(terms) if terms else '0'


if __name__ == '__main__':
    print('DEG =', DEG, '(total degree of the moments kept); Upsilon series exact through beta^%d' % (DEG - 2))
    print()
    print('(1) moments <chi^a cb^b> of SU(N) from the Toeplitz sum (nonzero ones, a+b <= DEG):')
    allmom = {}
    for N in range(2, 10):
        mom = moments(N, DEG)
        allmom[N] = mom
        nz = sorted((k, v) for k, v in mom.items() if v != 0)
        bad_centre = [k for k, v in nz if (k[0] - k[1]) % N]
        diag = {k[0]: v for k, v in nz if k[0] == k[1]}
        offd = [(k, v) for k, v in nz if k[0] != k[1]]
        print('  SU(%d): diagonal <|chi|^2k>, k=0..%d: %s' % (N, DEG // 2, [int(diag.get(k, 0)) for k in range(DEG // 2 + 1)]))
        print('         off-diagonal nonzero: %s;  violating the centre rule: %s' % (offd if offd else 'none', bad_centre if bad_centre else 'none'))
        if N >= 3:
            assert all(diag[k] == factorial(k) for k in range(4)), 'k! for k<=3 fails'
        if N <= DEG:
            assert mom.get((N, 0)) == 1 and mom.get((0, N)) == 1, '<chi^N> = 1 fails'
        assert not bad_centre
    print('  => checked: centre rule; <|chi|^{2k}> = k! for k <= min(N, DEG/2) (Catalan numbers for SU(2)); <chi^N> = 1.')
    print()
    order = DEG - 2
    print('(2)-(4) exact series in beta (through beta^%d):' % order)
    for N in range(3, 10):
        Up, Um, Voff = ups_series(allmom[N], N, order)
        diff, twoVoff = ups_direct(allmom[N], N, order)
        assert diff == twoVoff
        pred = Fr(2, (2 * N) ** (N - 2) * factorial(N - 2))
        print('  SU(%d): Upsilon_+ = 1 + %s' % (N, fmt([0] + Up[1:])))
        print('         Upsilon_- = 1 + %s' % fmt([0] + Um[1:]))
        print('         Upsilon_+ - Upsilon_- = %s ;  2(beta/2N)^(N-2)/(N-2)! = %s * b^%d' % (fmt(diff), pred, N - 2))
        lead = next((n for n, c in enumerate(diff) if c != 0), None)
        if N - 2 <= order:
            assert lead == N - 2 and diff[N - 2] == pred, 'leading term'
    print()
    # Specific claims
    d6 = ups_series(allmom[6], 6, order)
    d6 = [x - y for x, y in zip(d6[0], d6[1])]
    d7 = ups_series(allmom[7], 7, order)
    d7 = [x - y for x, y in zip(d7[0], d7[1])]
    d8 = ups_series(allmom[8], 8, order)
    d8 = [x - y for x, y in zip(d8[0], d8[1])]
    print('Lemma mg:lem:split: N=6 coefficient of beta^4 = %s (claimed 1/248832 = %s); lower orders zero: %s'
          % (d6[4], Fr(1, 248832), all(c == 0 for c in d6[:4])))
    print('                    N>=7: beta^0..beta^4 zero: N=7 %s, N=8 %s' % (all(c == 0 for c in d7[:5]), all(c == 0 for c in d8[:5])))
    print('Remark after the lemma: N=7 beta^5 coefficient = %s (claimed 1/32269440 = %s); N=8 beta^5 = %s, beta^6 = %s'
          % (d7[5], Fr(1, 32269440), d8[5], d8[6] if order >= 6 else 'n/a'))
    assert d6[4] == Fr(1, 248832) and d7[5] == Fr(1, 32269440) and d8[5] == 0
    print('N=6, next term: Upsilon_+ - Upsilon_- = %s' % fmt(d6, order))
    print()
    print('(5) Gaussian cancellation: moments delta_ab a! (all degrees <= DEG), N only enters through y = (beta/N) Re chi:')
    for N in (6, 7):
        gm = {(a, b): (Fr(factorial(a)) if a == b else Fr(0)) for a in range(DEG + 1) for b in range(DEG + 1) if a + b <= DEG}
        Up, Um, Voff = ups_series(gm, N, order)
        print('  N=%d: Upsilon_+ - 1 = %s ; Upsilon_- - 1 = %s ; Voff = %s' % (N, fmt([0] + Up[1:]), fmt([0] + Um[1:]), fmt(Voff)))
    print()
    # Values used in Corollary mg:cor:sectors(ii) (SU(3)) and Lemma mg:lem:Upsilon3
    Up3, Um3, _ = ups_series(allmom[3], 3, order)
    print('SU(3): Upsilon_+ through beta^3: %s (Lemma Upsilon3: 1 + b/6 - 5 b^3/1296)' % fmt(Up3, 3))
    print('       Upsilon_- through beta^3: %s (Lemma Upsilon3: 1 - b/6 + O(b^3))' % fmt(Um3, 3))
    assert Up3[1] == Fr(1, 6) and Up3[2] == 0 and Up3[3] == Fr(-5, 1296) and Um3[1] == Fr(-1, 6) and Um3[2] == 0
    import sympy as sp
    b = sp.symbols('b')
    Lp = sp.series(sp.log(sum(sp.Rational(c.numerator, c.denominator) * b**n for n, c in enumerate(Up3))), b, 0, 4).removeO()
    Lm = sp.series(sp.log(sum(sp.Rational(c.numerator, c.denominator) * b**n for n, c in enumerate(Um3))), b, 0, 4).removeO()
    print('       m_- - m_+ = log Ups_+ - log Ups_- + O(b^4) = %s + O(b^4)  (Corollary mg:cor:sectors(ii): b/3 + O(b^2))' % sp.expand(Lp - Lm))
    Up4, Um4, _ = ups_series(allmom[4], 4, order)
    Up5, Um5, _ = ups_series(allmom[5], 5, order)
    print('SU(4): Upsilon_+ - 1 = %s ; SU(5): Upsilon_+ - 1 = %s (Lemma Upsilon(iii), Upsilon3: b^2/128, b^3/6000 leading)'
          % (fmt([0] + Up4[1:], 4), fmt([0] + Up5[1:], 4)))
    print('All assertions passed.')

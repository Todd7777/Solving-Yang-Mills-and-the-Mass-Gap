# Referee check (ninth pass on Section 16) of Lemma mg:lem:upsbounds(i): e^{-2 beta} <= Upsilon(beta) <= e^beta, and 1 <= c_0 <= e^beta,
# by direct Weyl integration (independent of the authors' code), for U(1), SU(2), SU(3) (fundamental), at small and large beta.
# Upsilon = max over real unit psi in span(Re chi, Im chi) (unit in L^2(Haar)) of Var(psi) under w_beta/c_0, w_beta = exp((beta/N) Re chi)
# (Section 10, before Lemma mg:lem:Upsilon). With an orthonormal basis e_1, e_2 of the real span, Upsilon is the top eigenvalue of the 2x2 covariance matrix.
# Usage: python3 ups_bounds_check.py
from mpmath import mp, mpf, quad, cos, sin, exp, pi, sqrt, matrix, eigsy, cos as mcos

mp.dps = 30

def tilted_cov(N, beta, integrate, fns):
    """the covariance matrix of the functions fns (orthonormal under Haar) under w/c_0, and c_0"""
    w = lambda *th: exp(beta/N*fns[0](*th))      # fns[0] is Re chi
    c0 = integrate(lambda *th: w(*th))
    E = lambda h: integrate(lambda *th: h(*th)*w(*th))/c0
    basis = fns[1]
    m = [E(f) for f in basis]
    C = matrix(len(basis), len(basis))
    for i, f in enumerate(basis):
        for j, g in enumerate(basis):
            C[i, j] = E(lambda *th: f(*th)*g(*th)) - m[i]*m[j]
    return C, c0

def u1(beta):
    integ = lambda h: quad(lambda t: h(t), [-pi, 0, pi])/(2*pi)
    re = lambda t: cos(t)
    basis = [lambda t: sqrt(2)*cos(t), lambda t: sqrt(2)*sin(t)]
    return tilted_cov(1, beta, integ, (re, basis))

def su2(beta):
    integ = lambda h: quad(lambda t: h(t)*2/pi*sin(t)**2, [0, pi])          # class functions of SU(2), chi = 2 cos t
    re = lambda t: 2*cos(t)
    basis = [lambda t: 2*cos(t)]                                            # chi real: the span is one-dimensional
    return tilted_cov(2, beta, integ, (re, basis))

def su3(beta):
    # Weyl: eigenvalues e^{i a}, e^{i b}, e^{-i(a+b)}; density |Delta|^2/(3! (2 pi)^2) on [-pi, pi]^2
    def dens(a, b):
        z = [mp.expj(a), mp.expj(b), mp.expj(-a - b)]
        v = abs((z[0] - z[1])*(z[0] - z[2])*(z[1] - z[2]))**2
        return v/(6*(2*pi)**2)
    integ = lambda h: quad(lambda a, b: h(a, b)*dens(a, b), [-pi, 0, pi], [-pi, 0, pi])
    re = lambda a, b: cos(a) + cos(b) + cos(a + b)
    im = lambda a, b: sin(a) + sin(b) - sin(a + b)
    basis = [lambda a, b: sqrt(2)*re(a, b), lambda a, b: sqrt(2)*im(a, b)]
    return tilted_cov(3, beta, integ, (re, basis))

if __name__ == "__main__":
    ok_all = True
    for name, fn, betas in (('U(1)', u1, ('0.01', '0.3', '1', '3', '10')), ('SU(2)', su2, ('0.01', '0.3', '1', '3', '10')),
                            ('SU(3)', su3, ('0.01', '0.3', '1', '3'))):
        if name == 'SU(3)': mp.dps = 15
        for b in betas:
            beta = mpf(b)
            C, c0 = fn(beta)
            ev = eigsy(C)[0] if C.rows > 1 else [C[0, 0]]
            ups = max(ev)
            ok = exp(-2*beta) <= ups <= exp(beta) and 1 <= c0 <= exp(beta)
            ok_all = ok_all and ok
            print(f"{name:6s} beta = {b:>5s}: c_0 = {mp.nstr(c0, 10):>14s} in [1, e^beta]; Upsilon = {mp.nstr(ups, 10):>12s},"
                  f" e^(-2beta) = {mp.nstr(exp(-2*beta), 6):>10s}, e^beta = {mp.nstr(exp(beta), 6):>10s}: {ok}")
    print("all within the bounds of Lemma mg:lem:upsbounds(i):", ok_all)

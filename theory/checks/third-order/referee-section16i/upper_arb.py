# Referee check (ninth pass on Section 16) of Theorem mg:thm:thirdcap(e) in ball arithmetic (Arb, python-flint): E(beta)/beta^4 at the upper end of the
# range, for SU(3), SU(5), d = 3, 4, with the formula of Remark mg:rem:upperexact(ii) and Lemma mg:lem:uppertails (as re-implemented in upper_check.py),
# evaluated on the ball beta_0(d) (and at an exact binary number above it), against the constants c_+ of the theorem.
# Usage: python3 upper_arb.py
import math
from flint import arb, ctx
ctx.prec = 300
ONE = arb(1); E = ONE.exp()
COUNTS = {3: {'A': (2, {6: 8, 9: 10, 10: 47, 11: 60, 12: 96}), 'Q': (6, {6: 12, 9: 30, 10: 93, 11: 144, 12: 162})},
          4: {'A': (4, {6: 32, 9: 60, 10: 506}), 'Q': (18, {6: 60, 8: 24, 9: 262, 10: 1143})}}

def Eb(beta, d, N, iota=0):
    K = 4*(2*d - 3); eps = (2*beta).expm1(); x = E**2*K*eps
    th = x**5/(1 - x)
    w = eps*((K + 1)*th).exp(); q = E*K*w
    ubar = ((1 + iota)*beta + beta**2*beta.exp())/(2*N**2)
    s5 = (2 - iota)*N**2*ubar**5 + (beta.expm1() - beta)**4*beta.expm1()
    tau0 = arb(1 + iota)/2
    T = {}
    for S, nC, PS, Jp in (('Q', 6, 24*d - 42, 24*(d - 2)), ('A', 1, 8*d - 11, K)):
        ncap, Nn = COUNTS[d][S]; Nn = dict(Nn)
        if d == 4: Nn[11] = math.comb(Jp + K*11, 11)
        n0 = max(Nn)
        body = ncap*s5*(5*(K + 1)*th).exp() + sum((arb(c)*w**n for n, c in Nn.items()), arb(0)) + (1 + ONE/K)**Jp*q**(n0 + 1)/(1 - q)
        T[S] = (1 + eps)**nC*(th*PS).exp()*body
    TQ, TA = T['Q'], T['A']
    a = arb(2 - iota).sqrt()*N*ubar; a2 = (1 + TA)*a + TA
    B = TQ + eps**6 + 2*a*(1 + a)*TA + (1 + a2**2)*(16*beta).exp()*(arb(N)**4/tau0**2)*TQ
    assert bool(1 - B > 0)
    return TA*(2 + 3*a + 3*a**2) + B/(1 - B)

if __name__ == "__main__":
    cplus = {(3, 3): '0.042', (3, 4): '0.047', (5, 3): '0.318', (5, 4): '0.361'}
    ok = True
    for N in (3, 5):
        for d in (3, 4):
            K = 4*(2*d - 3)
            b0 = (ONE + ONE/(E**2*K*(K + 2))).log()/2            # beta_0(d) as a ball
            bup = arb(b0.upper())*(1 + arb('1e-30'))              # an exact number above beta_0(d)
            bup = arb(bup.upper())
            v = Eb(bup, d, N)/bup**4
            good = bool(arb(cplus[(N, d)]) - v > 0)
            ok = ok and good
            print(f"SU({N}), d={d}: E/beta^4 at beta = beta_0(d)(1+1e-30) (exact, above beta_0): {v.str(10)}  <  c_+ = {cplus[(N, d)]}: {good}")
    print("all certified (with E(beta)/beta^5 increasing, Remark mg:rem:upperexact(ii)):", ok)

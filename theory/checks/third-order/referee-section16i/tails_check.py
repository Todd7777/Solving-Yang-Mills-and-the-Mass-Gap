# Referee check (ninth pass on Section 16) of Theorem mg:thm:thirdcap(c): the relative truncation tail of EVERY infinite series that
# third_order_effective.lower3e (with tree_constants and marked) truncates, at beta_* in bound mode with the certified counts, computed as
# (sum to KX or NX) - (sum as truncated by the authors), relative to the truncated sum; the remainders beyond KX, NX are < 1e-100 relative.
# Families marked [EB] are among those bounded by error_budget.py; the others are not listed there.
# Usage: python3 tails_check.py GROUP d
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..'))
from mpmath import mp, mpf, e, exp, binomial, fsum, nstr
import my_eval as ME
import third_order_effective as T
import third_order_constants as TC

KT, NT = 60, 40     # the authors' truncations (fronts k < 60, clusters/polymers m, n < 40)

def main(G, d):
    mp.dps = 60
    P = ME.authors_params(G, d)
    r = ME.evaluate(G, d, P)
    beta = P['beta']; a, b, cJ = r['a'], r['b'], r['cJ']
    xF, xC, F, F6 = P['xF'], P['xC'], P['F'], P['F6']
    bb, ab = P['bb'], P['ab']
    n_ = d - 1; Pv = 2*n_*(n_ - 1); KS = 4*(2*d - 5); Tv = 4*(Pv - 1)
    Nv = ME.front_counts(d, 4*ME.NX + 8)
    tc = TC.tree_constants(bb, ab, cJ, Pv, KS, {k: Nv[k] for k in range(4, 200)}, 60, 40)
    xFb, xCb = tc[0], tc[1]
    def f(k, aa=a):
        if k == 4: return mpf(1)
        if k == 5: return mpf(0)
        if k == 6: return 36*F6*exp(6*xF)*aa**2
        return F*k**2*exp(xF*k)*aa**(k - 4)
    def X(n, k, c=4): return mpf(0) if k == 5 else min(binomial(4*n + c, k), (4*n + c)*Nv[k])
    fam = {}
    def rel(name, term, lo, hi_trunc, hi_ext):
        s_t = fsum(term(i) for i in range(lo, hi_trunc)); s_x = fsum(term(i) for i in range(hi_trunc, hi_ext))
        fam[name] = s_x/s_t if s_t else mpf(0)
    for j in range(3):
        rel(f'Y_{j} (actual b){" [EB]" if j == 0 else ""}', lambda m, j=j: m**j*binomial(KS*m, m - 1)*(b*exp(xC))**m, 1, NT, ME.NX)
        rel(f'A_{j} (actual a){" [EB]" if j == 0 else ""}', lambda k, j=j: k**j*Nv[k]*a**(mpf(k)/2)*exp(xF*k), 4, KT, ME.KX)
        rel(f'Y_{j} (bbar) [EB]', lambda m, j=j: m**j*binomial(KS*m, m - 1)*(bb*exp(xCb))**m, 1, NT, ME.NX)
        rel(f'A_{j} (abar) [EB]', lambda k, j=j: k**j*Nv[k]*ab**(mpf(k)/2)*exp(xFb*k), 4, KT, ME.KX)
    for j in range(2):
        rel(f'B_{j} (actual) [EB, via k^4]', lambda k, j=j: k**j*Nv[k]*f(k)*exp(xF*k), 4, KT, ME.KX)
        rel(f'B_{j} (barred) [EB, via k^4]', lambda k, j=j: k**j*Nv[k]*f(k, ab)*exp(xFb*k), 4, KT, ME.KX)
    for kmin in (4, 8):
        rel(f'one-plaquette A(k>={kmin})', lambda k: Nv[k]*a**(mpf(k)/2), kmin, KT, ME.KX)
        rel(f'one-plaquette F(k>={kmin})', lambda k: Nv[k]*f(k), kmin, KT, ME.KX)
    crude4 = lambda n: 2*binomial(4*Pv + Tv*n, n)*(n + 1)*(2*d - 3)*b**n
    crude_ex = lambda n: 2*binomial(4*Pv + Tv*n, n)*b**n*fsum(X(n, k)*f(k) for k in range(6, 4*n + 5))
    rel('full: crude Z_(n,4), n>=4', crude4, 4, NT, ME.NX)
    rel('full: crude exits, n>=1', crude_ex, 1, NT, ME.NX)
    if d == 4: rel('four: crude Z_(n,4), n>=13 [EB]', crude4, 13, NT, ME.NX)
    first_crude = 6 if d == 3 else 5
    rel(f'exits: crude, n>={first_crude} [EB]', crude_ex, first_crude, NT, ME.NX)
    rel('Sigma^(6): four-edge fronts, n>=2', lambda n: (2*d - 3)*(2*n + 3)*binomial(6*Pv + Tv*n, n)*b**n, 2, NT, ME.NX)
    rel('Sigma^(6): exits, n>=1 [EB]', lambda n: 2*binomial(6*Pv + Tv*n, n)*b**n*fsum(X(n, k, 6)*f(k) for k in range(6, 4*n + 7)), 1, NT, ME.NX)
    print(f"SU({G}), d={d}, beta_* = {nstr(beta, 10)}: relative truncation tails of every series of lower3e (bound mode, certified counts):")
    for k_, v in sorted(fam.items(), key=lambda kv: -kv[1]):
        print(f"    {k_:42s} {nstr(v, 3)}")
    print(f"  largest: {nstr(max(fam.values()), 3)} ({max(fam, key=fam.get)}); largest among those error_budget.py does not list: "
          f"{nstr(max(v for k_, v in fam.items() if '[EB' not in k_), 3)}")

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))

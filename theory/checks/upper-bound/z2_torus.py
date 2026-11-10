# Exact check of the reflection-positivity lower bound for Z_2 gauge theory, d=3,
# spatial torus 3x3, infinite time extent (transfer matrix on the physical space).
# Physical space = span of Walsh characters chi_S, S a Z_2-cycle of the spatial graph (dim 2^10).
import numpy as np, itertools, math, sys
from constants_upper import consts
L = 3
def lid(x, y, i):  # link from (x,y) in direction i (0: x-dir, 1: y-dir)
    return ((x % L)*L + (y % L))*2 + i
nl = 2*L*L
def plaq_links(x, y):
    return [lid(x,y,0), lid(x+1,y,1), lid(x,y+1,0), lid(x,y,1)]
plaqs = [plaq_links(x,y) for x in range(L) for y in range(L)]
pmask = [sum(1<<l for l in p) for p in plaqs]
# cycle space generators: 8 plaquettes + 2 winding cycles
wind = [sum(1<<lid(x,0,0) for x in range(L)), sum(1<<lid(0,y,1) for y in range(L))]
gens = pmask[:-1] + wind
cycles = []
for bits in range(1<<len(gens)):
    m = 0
    for k,g in enumerate(gens):
        if bits>>k & 1: m ^= g
    cycles.append(m)
assert len(set(cycles)) == 1024
idx = {m:i for i,m in enumerate(cycles)}
n = len(cycles)
popc = np.array([bin(m).count('1') for m in cycles])
def shift(mask):  # multiplication by chi_mask (a cycle) permutes the basis
    M = np.zeros((n,n))
    for i,m in enumerate(cycles): M[idx[m ^ mask], i] = 1.0
    return M
Mp = [shift(m) for m in pmask]
p0 = plaqs[0]; E0 = pmask[0]   # A = sigma_{p0}; the four sides sit over the links of p0
A = Mp[0]
def run(beta):
    ch, sh = math.cosh(beta), math.sinh(beta)
    V = np.eye(n)
    for M in Mp: V = V @ (math.cosh(beta/2)*np.eye(n) + math.sinh(beta/2)*M)
    Kd = ch**(nl-popc) * sh**popc            # K_1 P on cycles
    T = V @ np.diag(Kd) @ V
    w, U = np.linalg.eigh(T)
    lam0, lam1 = w[-1], w[-2]
    Om = U[:,-1]; Om = Om*np.sign(Om[idx[0]])
    a = Om @ A @ Om
    psi = A @ Om - a*Om                        # (A - a) Omega
    var = psi @ psi
    cov = psi @ T @ psi / lam0
    # decomposition by F = S cap E0
    Fsets = np.array([m & E0 for m in cycles])
    terms = {}
    Vpsi = V @ psi
    for F in sorted(set(Fsets.tolist())):
        sel = (Fsets == F)
        terms[F] = np.sum(Kd[sel]*Vpsi[sel]**2)/lam0
    G = terms[E0]
    u = math.tanh(beta)
    gap = math.log(lam0/lam1)
    return dict(lam0=lam0, gap=gap, a=a, var=var, cov=cov, G=G, terms=terms, u=u)
print("beta     u^4         Cov/u^4    <G>/u^4    min term     Var       log(Var/Cov)-4log(1/u)  gap/log(1/u)")
for beta in [0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 1e-3, 4e-4, 1e-4]:
    r = run(beta)
    u = r['u']
    mn = min(r['terms'].values())
    s = sum(r['terms'].values())
    assert abs(s - r['cov']) <= 1e-6*abs(s), (s, r['cov'])
    print(f"{beta:<8g} {u**4:.4e}  {r['cov']/u**4:.6f}  {r['G']/u**4:.6f}  {mn:+.3e}  {r['var']:.6f}  {math.log(r['var']/r['cov'])-4*math.log(1/u):+.3e}   {r['gap']/math.log(1/u):.4f}")
# theorem's constants for Z_2 viewed as N=1, iota=1, d=3 (valid for beta < beta0(3) = 4.03e-4)
for beta in [4e-4, 1e-4]:
    c = consts(beta, 3, 1, 1)
    r = run(beta)
    u = r['u']
    print(f"beta={beta}: <G>/(tau0 u^4)={r['G']/u**4:.8f} >= L={float(c['L']):.8f};  Var={r['var']:.8f} <= V={float(c['V']):.8f};  log(Var/Cov)-4log(1/u)={math.log(r['var']/r['cov'])-4*math.log(1/u):+.3e} <= D={float(c['D']):.3e}")

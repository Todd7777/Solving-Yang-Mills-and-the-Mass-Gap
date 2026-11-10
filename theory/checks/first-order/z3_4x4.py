# Independent referee check (Section 9 referee report, 10 November 2026), added to the repository with its author's code unchanged
# below this header. Z_3 on the 4x4 spatial torus, which satisfies the hypotheses (sides >= 4), unlike the 3x3 torus of
# zn_first_order.py. Memory: several arrays of 3^16 doubles (about 2 GB in all). Usage: python3 z3_4x4.py 0.05
"""Referee check: Z_3 gauge theory, d=3, on the 4x4 spatial torus (sides >= 4, within the hypotheses of the theorems),
sector of zero electric winding (functions of the 16 plaquette variables with sum B = 0, 3^15 states).
Transfer matrix T = V^{1/2} K V^{1/2} as in zn_exact.py, applied by FFT on Z_3^16 (the constraint is preserved).
Power iteration in the zero-momentum, rotation-invariant sectors: C-odd (start sum_p sin(2 pi B_p/3)),
C-even (start sum_p cos(2 pi B_p/3), vacuum deflated).  Reports (m - 4 log(1/u))/u.
Usage: python3 z3_4x4.py beta
"""
import sys
import time
import numpy as np
import scipy.fft as sfft

n, L = 3, 4
P = L * L
beta = float(sys.argv[1])
h = np.arange(n)
w = np.exp(beta * np.cos(2 * np.pi * h / n))
c0 = w.mean()
onef = w / c0
ur = np.array([np.mean(w * np.exp(-2j * np.pi * r * h / n)).real / c0 for r in range(n)])
u = ur[1]
pid = lambda x, y: (x % L) * L + (y % L)
pairs = []
for x in range(L):
    for y in range(L):
        pairs.append((pid(x, y), pid(x, y - 1)))
        pairs.append((pid(x, y), pid(x - 1, y)))
shape = (n,) * P


def axis_vec(p, vals):
    s = [1] * P
    s[p] = n
    return vals.reshape(s)


t0 = time.time()
# V^{1/2}(B) and mask sum B = 0 mod 3
Vh = np.ones(shape)
S = np.zeros(shape, dtype=np.int8)
for p in range(P):
    Vh = Vh * axis_vec(p, np.sqrt(onef))
    S = (S + axis_vec(p, np.arange(n, dtype=np.int8))) % n
mask = (S == 0)
del S
Vh *= mask
# Khat(m) = prod_pairs ur[(m_p - m_q) mod 3]
Kh = np.ones(shape)
for (p, q) in pairs:
    d = (np.arange(n)[:, None] - np.arange(n)[None, :]) % n
    fac = ur[d]  # (n,n) as function of (m_p, m_q)
    s = [1] * P
    s[p] = n
    s[q] = n
    if p < q:
        Kh = Kh * fac.reshape(s)
    else:
        Kh = Kh * fac.T.reshape(s)
print(f"setup {time.time()-t0:.1f}s, u={u:.8e}", flush=True)


def T(x):
    y = Vh * x
    Y = sfft.fftn(y, workers=4)
    Y *= Kh
    y = sfft.ifftn(Y, workers=4).real
    return Vh * y


def Cflip(x):
    # B -> -B : index i -> (-i) mod 3 along every axis
    idx = np.array([0, 2, 1])
    for ax in range(P):
        x = np.take(x, idx, axis=ax)
    return x


def start(kind):
    s = np.zeros(shape)
    for p in range(P):
        s = s + axis_vec(p, (np.sin if kind == 'odd' else np.cos)(2 * np.pi * np.arange(n) / n))
    return s * mask


def power(x, deflate=None, cpar=None, iters=12, tol=1e-13):
    lam_old = None
    for it in range(iters):
        if cpar is not None:
            x = 0.5 * (x + cpar * Cflip(x))
        if deflate is not None:
            x = x - np.vdot(deflate, x) * deflate
        x /= np.linalg.norm(x)
        y = T(x)
        lam = np.vdot(x, y)
        x = y
        if lam_old is not None and abs(lam - lam_old) < tol * abs(lam):
            break
        lam_old = lam
    x /= np.linalg.norm(x)
    return lam, x


t0 = time.time()
lam0, v0 = power(start('even') + mask * 1.0, cpar=+1, iters=8)
print(f"vacuum lambda0={lam0:.15e} ({time.time()-t0:.0f}s)", flush=True)
lamo, _ = power(start('odd'), cpar=-1, iters=10)
mo = np.log(lam0 / lamo)
print(f"C-odd : lambda/lambda0={lamo/lam0:.12e}  (m-4log(1/u))/u = {(mo-4*np.log(1/u))/u:+.7f} ({time.time()-t0:.0f}s)", flush=True)
lame, _ = power(start('even'), deflate=v0, cpar=+1, iters=10)
me = np.log(lam0 / lame)
print(f"C-even: lambda/lambda0={lame/lam0:.12e}  (m-4log(1/u))/u = {(me-4*np.log(1/u))/u:+.7f} ({time.time()-t0:.0f}s)", flush=True)

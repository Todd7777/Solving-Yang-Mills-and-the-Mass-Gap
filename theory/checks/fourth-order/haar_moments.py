"""Monte Carlo check (Haar measure on SU(6), SU(7)) of the character moments used in upsilon_split.py: m11=1, m22=2, m33=6, <chi^N>=1, <chi^3>=0."""
import numpy as np
rng = np.random.default_rng(1)
def haar_su(N, n):
    Z = (rng.standard_normal((n, N, N)) + 1j * rng.standard_normal((n, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R, axis1=1, axis2=2); Q = Q * (d / np.abs(d))[:, None, :]
    det = np.linalg.det(Q); Q = Q / (det ** (1.0 / N))[:, None, None]
    return Q
for N in (6, 7):
    ch = np.concatenate([np.trace(haar_su(N, 200000), axis1=1, axis2=2) for _ in range(5)])
    a2 = np.abs(ch) ** 2
    print(N, 'm11 %.3f m22 %.3f m33 %.2f  Re<chi^N> %.3f  |<chi^3>| %.3f' % (a2.mean(), (a2**2).mean(), (a2**3).mean(), (ch**N).real.mean(), np.abs((ch**3).mean())))

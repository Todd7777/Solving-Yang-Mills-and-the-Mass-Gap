# Referee check (tenth pass on Section 16): a numerical sanity check (Numerical, not part of any proof) of
#   (a) Lemma mg:lem:upsbounds(i), the new sentence: e^{-2 beta} <= U_+-(beta) <= e^beta for SU(3), SU(4), SU(5) (iota = 0), also far outside the range;
#   (b) the exact series of series_check.py against the Weyl integration formula at small beta (U_+ - U_-, -log U_+-).
# U_+ = Var(sqrt2 Re chi), U_- = Var(sqrt2 Im chi) under exp((beta/N) Re chi) dU / c_0, by the Weyl integration formula on the maximal torus
# (eigenvalues e^{i theta_j}, sum theta_j = 0 mod 2 pi, density prop. to |Vandermonde|^2), with the trapezoidal rule on a uniform grid of the torus
# (exponentially accurate for these smooth periodic integrands).
# Usage: python3 weyl_check.py
import os, sys
import numpy as np
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)

def torus(N, M):
    """eigenvalue grid of the maximal torus of SU(N) (N-1 free angles, M points each) and the Weyl density |Delta|^2"""
    th = np.linspace(0, 2*np.pi, M, endpoint=False)
    grids = np.meshgrid(*([th]*(N - 1)), indexing='ij')
    angles = [g.ravel() for g in grids]
    angles.append(-sum(angles))
    z = np.exp(1j*np.array(angles))                       # N x points
    dens = np.ones(z.shape[1])
    for i in range(N):
        for j in range(i + 1, N): dens *= np.abs(z[i] - z[j])**2
    return z.sum(axis=0), dens/dens.sum()

def upsilons(N, beta, M):
    chi, w = torus(N, M)
    wt = w*np.exp((beta/N)*chi.real); wt /= wt.sum()
    re, im = np.sqrt(2)*chi.real, np.sqrt(2)*chi.imag
    Up = (wt*re**2).sum() - (wt*re).sum()**2
    Um = (wt*im**2).sum() - (wt*im).sum()**2
    norm_re, norm_im = (w*re**2).sum(), (w*im**2).sum()      # Haar: both must be 1 (iota = 0)
    return Up, Um, norm_re, norm_im

if __name__ == "__main__":
    import series_check as SC
    grids = {3: 96, 4: 48, 5: 24}
    print("(a) e^{-2 beta} <= U_+-(beta) <= e^beta (Lemma mg:lem:upsbounds(i)):")
    for N in (3, 4, 5):
        M = grids[N]
        for beta in (0.01, 0.3, 1.0, 3.0, 10.0):
            Up, Um, n1, n2 = upsilons(N, beta, M)
            ok = all(np.exp(-2*beta) <= U <= np.exp(beta) for U in (Up, Um))
            print(f"  SU({N}), beta = {beta:5}: U_+ = {Up:.10f}, U_- = {Um:.10f}  (Haar norms {n1:.12f}, {n2:.12f}; grid {M}^{N - 1});  within [e^-2b, e^b]: {ok}")
    print("(b) the exact series against the Weyl integration at small beta:")
    for N in (3, 4, 5):
        s = SC.sector(N, 12); M = grids[N]
        for beta in (0.05, 0.2):
            Up, Um, _, _ = upsilons(N, beta, M)
            ser = lambda L: float(sum(Fr(c)*Fr(beta).limit_denominator(10**6)**k for k, c in enumerate(L)))
            print(f"  SU({N}), beta = {beta}: U_+ - U_- (Weyl) = {Up - Um:.12e}, (series to beta^10) = {ser(s['Up']) - ser(s['Um']):.12e};"
                  f"  -log U_+: {-np.log(Up):.12e} vs {ser(s['neglog_plus']):.12e};  -log U_-: {-np.log(Um):.12e} vs {ser(s['neglog_minus']):.12e}")

"""Exact transfer matrix of two-dimensional SU(2) lattice Yang--Mills on a spatial circle.

In d = 2 the spatial lattice is a circle of L sites with links U_1, ..., U_L, there are no spatial
plaquettes (V = 1 in the notation of theory/lattice-foundations.tex, Definition def:T), and the
gauge-invariant Hilbert space consists of class functions of the holonomy h = U_1 U_2 ... U_L.
The functions chi_j(h) (j = 0, 1/2, 1, ...) are orthonormal, and by the convolution identity of
Peter--Weyl theory each link kernel w_beta(U V^{-1}) acts on chi_j(h) as multiplication by
c_j/d_j.  Hence the transfer matrix is diagonal:

    T chi_j(h) = (c_j / d_j)^L chi_j(h),     lambda_j / lambda_0 = (I_{2j+1}(beta) / I_1(beta))^L ,

and the exact finite-lattice mass gap is  m = L * log(I_1(beta) / I_2(beta)) = L * sigma(beta):
the lightest gauge-invariant excitation is a flux loop winding around space, with the exact 2-d
string tension sigma = -log u(beta).

The operator A = (1/2) chi_{1/2}(h) (half the trace of the spatial Polyakov loop) acts by
Clebsch--Gordan: chi_{1/2} chi_j = chi_{j-1/2} + chi_{j+1/2} (j >= 1/2), chi_{1/2} chi_0 = chi_{1/2}.
The exact correlator on a periodic lattice with L0 time slices is

    < A(tau) A(0) > = Tr(D^{L0 - tau} A D^tau A) / Tr(D^{L0}),   D = diag(lambda_j / lambda_0),

which this module evaluates in a truncated character basis (the truncation error is
super-exponentially small).  Used to test the full Monte Carlo -> correlator -> effective-mass
pipeline against an exactly known gap.
"""

from __future__ import annotations

import numpy as np
from scipy.special import ive


def eigenvalue_ratios(beta: float, L: int, jmax2: int = 40) -> np.ndarray:
    """lambda_j / lambda_0 for 2j = 0, 1, ..., jmax2."""
    nus = np.arange(1, jmax2 + 2)
    a = ive(nus, beta) / ive(1, beta)          # c_j / (d_j c_0) = I_{2j+1} / I_1
    return a ** L


def polyakov_operator(jmax2: int = 40) -> np.ndarray:
    """Matrix of A = (1/2) chi_{1/2}(h) in the orthonormal basis chi_j, 2j = 0..jmax2."""
    n = jmax2 + 1
    A = np.zeros((n, n))
    for k in range(n - 1):
        A[k, k + 1] = A[k + 1, k] = 0.5
    return A


def exact_gap(beta: float, L: int) -> float:
    """Exact mass gap -log(lambda_{1/2} / lambda_0) = L log(I_1/I_2)."""
    return float(L * np.log(ive(1, beta) / ive(2, beta)))


def exact_correlator(beta: float, L: int, L0: int, jmax2: int = 40) -> np.ndarray:
    """Exact <A(tau) A(0)> for tau = 0..L0-1 on the periodic L0 x L lattice (<A> = 0 exactly)."""
    lam = eigenvalue_ratios(beta, L, jmax2)
    A = polyakov_operator(jmax2)
    Z = np.sum(lam ** L0)
    out = np.empty(L0)
    for tau in range(L0):
        M = (lam ** (L0 - tau))[:, None] * A * (lam ** tau)[None, :]   # D^{L0-tau} A D^tau
        out[tau] = np.sum(M * A.T) / Z                                 # Tr(M A)
    return out

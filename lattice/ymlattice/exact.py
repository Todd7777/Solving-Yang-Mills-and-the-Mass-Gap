"""Exact results used to validate the simulation and the expansions in the theory notes.

Character expansion of the SU(2) plaquette weight
-------------------------------------------------
With ``P = (1/2) Tr U`` and Haar density ``(2/pi) sin^2(theta) dtheta`` on the class angle,

    exp(beta P) = sum_j d_j chat_j(beta) chi_j(U),   d_j = 2j+1,   chat_j(beta) = 2 I_{2j+1}(beta)/beta,

because ``int chi_j(U) e^{beta cos theta} dU = I_{2j}(beta) - I_{2j+2}(beta) = (2(2j+1)/beta) I_{2j+1}(beta)``.
The *fundamental coefficient ratio* ``u(beta) = chat_{1/2}/chat_0 = I_2(beta)/I_1(beta)`` is the
expansion parameter of the strong-coupling expansion (``u = beta/4 + O(beta^3)``).

Two-dimensional torus
---------------------
On a periodic two-dimensional lattice with ``V`` plaquettes the character expansion is exact:
``Z = sum_j chat_j(beta)^V`` (the Euler characteristic of the torus is zero), hence

    <P> = (1/V) d/dbeta log Z = sum_j chat_j^{V-1} chat_j' / sum_j chat_j^V  --->  I_2(beta)/I_1(beta)  (V -> inf).

Weak coupling
-------------
Equipartition for the quadratic fluctuations (``(d-1)(N^2-1)`` physical modes per site after
removing gauge directions) gives ``<1 - P> = (N^2 - 1)/(d beta) + O(beta^-2)`` in ``d`` dimensions,
i.e. ``1 - (N^2-1)/(4 beta)`` in four dimensions.
"""

from __future__ import annotations

import numpy as np
from scipy.special import iv, ive


# ----------------------------------------------------------------------------- SU(2)
def su2_chat(j2: int, beta: float) -> float:
    """``chat_j(beta) = 2 I_{2j+1}(beta)/beta`` for ``j2 = 2j`` (integer)."""
    return 2.0 * iv(j2 + 1, beta) / beta


def su2_u(beta: float) -> float:
    """Fundamental character coefficient ratio ``u(beta) = I_2(beta)/I_1(beta)``."""
    return ive(2, beta) / ive(1, beta)


def su2_plaquette_2d(beta: float, V: int | None = None, jmax2: int = 200) -> float:
    """Exact plaquette expectation of SU(2) on a 2-d torus with ``V`` plaquettes (``V=None``: infinite)."""
    if V is None:
        return su2_u(beta)
    # work with scaled Bessel functions to avoid overflow: chat_j = 2 e^beta ive(2j+1)/beta
    nus = np.arange(1, jmax2 + 2)                       # nu = 2j+1
    c = 2.0 * ive(nus, beta) / beta                      # chat_j * e^{-beta}
    # d/dbeta [2 I_nu(beta)/beta] = 2 I_nu'(beta)/beta - 2 I_nu(beta)/beta^2, I_nu' = (I_{nu-1}+I_{nu+1})/2
    cp = (ive(nus - 1, beta) + ive(nus + 1, beta)) / beta - 2.0 * ive(nus, beta) / beta ** 2
    logc = np.log(np.maximum(c, 1e-300))
    w = np.exp((V - 1) * (logc - logc[0]))               # chat_j^{V-1} / chat_0^{V-1}
    num = np.sum(w * cp)
    den = np.sum(w * c)
    return float(num / den)


# ----------------------------------------------------------------------------- SU(3)
def _su3_weyl_grid(n: int = 400):
    """Angles and Weyl density for class functions on SU(3).

    Eigenvalues ``e^{i t1}, e^{i t2}, e^{-i(t1+t2)}``; the Haar measure restricted to class
    functions is proportional to ``prod_{i<j} |e^{i t_i} - e^{i t_j}|^2 dt1 dt2``."""
    t = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    t1, t2 = np.meshgrid(t, t, indexing="ij")
    t3 = -(t1 + t2)
    e1, e2, e3 = np.exp(1j * t1), np.exp(1j * t2), np.exp(1j * t3)
    dens = np.abs(e1 - e2) ** 2 * np.abs(e1 - e3) ** 2 * np.abs(e2 - e3) ** 2
    trace = e1 + e2 + e3
    return trace, dens / dens.sum()


def su3_one_link(beta: float, n: int = 400):
    """One-link SU(3) integrals with weight ``exp((beta/3) Re Tr U)``.

    Returns ``(chat_0, <P>)`` where ``chat_0 = int dU exp((beta/3) Re Tr U)`` and
    ``<P> = <(1/3) Re Tr U>`` under that weight (= the exact 2-d infinite-volume plaquette)."""
    trace, w = _su3_weyl_grid(n)
    P = trace.real / 3.0
    boltz = np.exp(beta * P)
    Z = np.sum(w * boltz)
    return float(Z), float(np.sum(w * boltz * P) / Z)


def su3_plaquette_2d(beta: float) -> float:
    """Exact infinite-volume 2-d plaquette for SU(3)."""
    return su3_one_link(beta)[1]


def su3_u(beta: float) -> float:
    """Fundamental character coefficient ratio for SU(3): ``<(1/3) Tr U>`` under the one-link weight
    (equals ``chat_fund / chat_0`` since ``<Tr U> = <Tr U^+>`` is real)."""
    return su3_one_link(beta)[1]


# ----------------------------------------------------------------------------- series
def plaquette_weak(beta: float, N: int, d: int = 4) -> float:
    """One-loop weak-coupling plaquette ``1 - (N^2-1)/(d beta)``."""
    return 1.0 - (N * N - 1.0) / (d * beta)


def plaquette_strong_leading(beta: float, N: int) -> float:
    """Leading strong-coupling plaquette ``u(beta)`` (fundamental character coefficient ratio)."""
    return su2_u(beta) if N == 2 else su3_u(beta)


def _u_derivative(beta: float, N: int, h: float = 1e-5) -> float:
    if N == 2:
        # u = I_2/I_1 ; u' = I_2'/I_1 - I_2 I_1'/I_1^2 with I_n' = (I_{n-1} + I_{n+1})/2
        i0, i1, i2, i3 = (ive(n, beta) for n in range(4))
        return (0.5 * (i1 + i3)) / i1 - i2 * (0.5 * (i0 + i2)) / i1 ** 2
    return (su3_u(beta + h) - su3_u(beta - h)) / (2 * h)


def plaquette_strong_nlo(beta: float, N: int, d: int = 4) -> float:
    """Strong-coupling plaquette through the first polymer (the surface of a unit cube).

    ``log Z = N_p log c_0 + N_cubes z(cube) + ...`` with ``z(cube) = sum_{r != 0} (chat_r/chat_0)^6 d_r^2``
    (a closed surface of Euler characteristic 2 with all faces in the representation ``r``).  Keeping the
    lowest representations, ``z = 4u^6`` for SU(2) and ``z = 2 * 9 u^6 = 18u^6`` for SU(3) (fundamental and
    anti-fundamental).  With ``N_cubes/N_p = (d-2)/3`` cubes per plaquette and ``<P> = (1/N_p) d log Z / d beta``:

        <P> = u + (d-2)/3 * 6 z_6 u^5 u'(beta) + O(u^9),   z_6 = 4 (SU(2)), 18 (SU(3)).

    In four dimensions: ``u + 16 u^5 u'`` (SU(2)), ``u + 72 u^5 u'`` (SU(3))."""
    u = plaquette_strong_leading(beta, N)
    z6 = 4.0 if N == 2 else 18.0
    return u + (d - 2) / 3.0 * 6.0 * z6 * u ** 5 * _u_derivative(beta, N)


def strong_coupling_glueball_mass(beta: float, N: int = 2) -> float:
    """Leading strong-coupling 0++ mass ``a m = -4 log u(beta)`` (a plaquette propagating in time costs
    the four side plaquettes of a cube per time step)."""
    return -4.0 * np.log(plaquette_strong_leading(beta, N))


# ----------------------------------------------------------------------------- RG / scaling
def b0(N: int) -> float:
    """One-loop beta-function coefficient ``b0 = 11 N / (48 pi^2)`` (pure gauge)."""
    return 11.0 * N / (48.0 * np.pi ** 2)


def b1(N: int) -> float:
    """Two-loop coefficient ``b1 = 34 N^2 / (3 (16 pi^2)^2)``."""
    return 34.0 * N * N / (3.0 * (16.0 * np.pi ** 2) ** 2)


def a_lambda_lattice(beta: float, N: int) -> float:
    """Two-loop asymptotic-scaling prediction for ``a Lambda_L`` as a function of ``beta = 2N/g^2``:

        a Lambda_L = (b0 g^2)^{-b1/(2 b0^2)} exp(-1/(2 b0 g^2)).

    This is the function the lattice mass gap must track from above: ``a m(beta) >= c * a Lambda_L``
    is the quantitative form of the continuum mass gap along the Wilson-action trajectory."""
    g2 = 2.0 * N / beta
    return (b0(N) * g2) ** (-b1(N) / (2.0 * b0(N) ** 2)) * np.exp(-1.0 / (2.0 * b0(N) * g2))

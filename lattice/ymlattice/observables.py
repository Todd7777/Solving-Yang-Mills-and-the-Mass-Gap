"""Gauge-invariant observables: Wilson loops, Creutz ratios, Polyakov loops, APE smearing,
zero-momentum 0++ glueball correlators and effective masses."""

from __future__ import annotations

import numpy as np

from . import su2, su3
from .lattice import Lattice


def _group(lat: Lattice):
    return su2 if lat.group == "su2" else su3


def line(lat: Lattice, U: np.ndarray, mu: int, R: int) -> np.ndarray:
    """``L_mu^R(x) = U_mu(x) U_mu(x+mu) ... U_mu(x+(R-1)mu)``."""
    G = _group(lat)
    L = U[mu]
    for k in range(1, R):
        L = G.mul(L, lat.shift(U[mu], mu, +k))
    return L


def wilson_loop_field(lat: Lattice, R: int, T: int, mu: int, nu: int, U: np.ndarray | None = None) -> np.ndarray:
    """``(1/N) Re Tr`` of the ``R x T`` rectangular loop in the ``(mu, nu)`` plane based at every site."""
    G = _group(lat)
    U = lat.U if U is None else U
    Lmu = line(lat, U, mu, R)
    Lnu = line(lat, U, nu, T)
    a = G.mul(Lmu, lat.shift(Lnu, mu, +R))        # along mu then nu
    b = G.mul(Lnu, lat.shift(Lmu, nu, +T))        # along nu then mu; W = a b^+
    return G.retrace(G.mul_dag(a, b))


def wilson_loop(lat: Lattice, R: int, T: int, planes: str = "all", U: np.ndarray | None = None) -> float:
    """Average ``R x T`` Wilson loop.  ``planes='all'``: all ``(mu<nu)`` planes with both orientations
    (R along mu, T along nu and vice versa); ``planes='temporal'``: ``mu = 0`` only."""
    vals = []
    for mu in range(lat.d):
        for nu in range(lat.d):
            if mu == nu:
                continue
            if planes == "temporal" and mu != 0:
                continue
            if planes == "all" and nu < mu and R == T:
                continue
            vals.append(float(np.mean(wilson_loop_field(lat, R, T, mu, nu, U))))
    return float(np.mean(vals))


def creutz_ratio(W: dict, R: int, T: int) -> float:
    """``chi(R,T) = -log[ W(R,T) W(R-1,T-1) / (W(R,T-1) W(R-1,T)) ]`` -> string tension for large loops."""
    return -np.log(W[(R, T)] * W[(R - 1, T - 1)] / (W[(R, T - 1)] * W[(R - 1, T)]))


def polyakov_loop(lat: Lattice) -> complex:
    """Spatial average of ``(1/N) Tr`` of the temporal Polyakov line (direction 0)."""
    G = _group(lat)
    L = line(lat, lat.U, 0, lat.shape[0])
    L0 = np.take(L, 0, axis=0)                       # the loop is the same at every t; take t = 0
    if lat.group == "su2":
        return complex(np.mean(su2.retrace(L0)))
    tr = np.trace(L0, axis1=-2, axis2=-1) / 3.0
    return complex(np.mean(tr))


# ----------------------------------------------------------------------------- smearing
def staple_sum_dirs(lat: Lattice, U: np.ndarray, mu: int, dirs) -> np.ndarray:
    """Staple sum of link ``mu`` using only the directions in ``dirs`` (same orientation as
    :meth:`Lattice.staple_sum`)."""
    G = _group(lat)
    A = None
    for nu in dirs:
        if nu == mu:
            continue
        Unu_xmu = lat.shift(U[nu], mu, +1)
        Umu_xnu = lat.shift(U[mu], nu, +1)
        up = G.mul(Unu_xmu, G.dagger(G.mul(U[nu], Umu_xnu)))
        Unu_xmnu = lat.shift(U[nu], nu, -1)
        Umu_xmnu = lat.shift(U[mu], nu, -1)
        Unu_xmu_mnu = lat.shift(Unu_xmu, nu, -1)
        down = G.mul(G.dagger(G.mul(Umu_xmnu, Unu_xmu_mnu)), Unu_xmnu)
        A = up + down if A is None else A + up + down
    return A


def ape_smear(lat: Lattice, alpha: float = 0.5, n: int = 1, U: np.ndarray | None = None) -> np.ndarray:
    """APE smearing of the *spatial* links (directions 1..d-1), ``n`` iterations:

        U_i(x) -> Proj_{SU(N)}[ (1 - alpha) U_i(x) + alpha/(2(d-2)) * sum_{j spatial != i} staples_j(x, i) ]

    The staple sum ``A`` is oriented as in :meth:`Lattice.staple_sum` (``Tr U A`` = plaquettes), so the
    link-like object to add is ``A^+``.  Temporal links are untouched.  Returns a new link field."""
    G = _group(lat)
    U = (lat.U if U is None else U).copy()
    spatial = list(range(1, lat.d))
    nst = 2 * (lat.d - 2)
    for _ in range(n):
        new = U.copy()
        for mu in spatial:
            A = staple_sum_dirs(lat, U, mu, spatial)
            cand = (1.0 - alpha) * U[mu] + (alpha / nst) * G.dagger(A)
            if lat.group == "su2":
                new[mu] = su2.normalize(cand)
            else:
                new[mu] = su3.project_su3(cand)
        U = new
    return U


# ----------------------------------------------------------------------------- glueball
def spatial_plaquette_timeslices(lat: Lattice, U: np.ndarray | None = None) -> np.ndarray:
    """``O(t) = sum_{x spatial} sum_{i<j spatial} (1/N) Re Tr U_ij(x, t)`` (zero-momentum 0++ operator)."""
    G = _group(lat)
    U = lat.U if U is None else U
    out = np.zeros(lat.shape[0])
    for mu in range(1, lat.d):
        for nu in range(mu + 1, lat.d):
            a = G.mul(U[mu], lat.shift(U[nu], mu, +1))
            b = G.mul(U[nu], lat.shift(U[mu], nu, +1))
            p = G.retrace(G.mul_dag(a, b))
            out += p.reshape(lat.shape[0], -1).sum(axis=1)
    return out


def timeslice_correlator(O: np.ndarray) -> np.ndarray:
    """Connected zero-momentum correlator from time-slice operators ``O[meas, t]``:

        C(tau) = <(1/L0) sum_t O(t) O(t+tau)> - <O>^2 ."""
    O = np.asarray(O, dtype=float)
    L0 = O.shape[1]
    mean = O.mean()
    C = np.array([np.mean(O * np.roll(O, -tau, axis=1)) for tau in range(L0)])
    return C - mean * mean


def effective_mass(C: np.ndarray, periodic: bool = True) -> np.ndarray:
    """``m_eff(tau)`` from ``C(tau)/C(tau+1)``; for periodic time the cosh form
    ``C(tau) ~ cosh(m (tau - L0/2))`` is solved numerically."""
    C = np.asarray(C, dtype=float)
    L0 = len(C)
    m = np.full(L0 - 1, np.nan)
    for tau in range(L0 - 1):
        if C[tau] <= 0 or C[tau + 1] <= 0 or C[tau + 1] >= C[tau]:
            continue
        m0 = np.log(C[tau] / C[tau + 1])
        if not periodic:
            m[tau] = m0
            continue
        # solve cosh(m (tau - L0/2)) / cosh(m (tau+1 - L0/2)) = C[tau]/C[tau+1] by bisection
        target = C[tau] / C[tau + 1]
        f = lambda mm: np.cosh(mm * (tau - L0 / 2)) / np.cosh(mm * (tau + 1 - L0 / 2)) - target
        lo, hi = 1e-8, max(2 * m0, 1.0)
        if f(lo) * f(hi) > 0:
            m[tau] = m0
            continue
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if f(lo) * f(mid) <= 0:
                hi = mid
            else:
                lo = mid
        m[tau] = 0.5 * (lo + hi)
    return m

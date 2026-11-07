"""Analysis helpers for glueball / Wilson-loop runs produced by ``scripts/glueball.py``.

* :func:`load_runs` pools independent runs (same beta, lattice, smearing levels).
* :func:`corr_matrix` builds the connected, symmetrised correlation matrix of the time-slice
  operators at all smearing levels.
* :func:`gevp_lambda0` solves C(tau) v = lambda C(tau0) v and returns the largest eigenvalue.
* :func:`effective_masses` gives single-operator and GEVP effective masses with block-jackknife
  errors.

Reflection positivity makes every single-operator correlator a positive combination of
exponentials (log-convex), so effective masses approach the gap from ABOVE: Monte Carlo
effective masses estimate upper bounds on the gap, never lower bounds.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

from . import observables as obs
from .stats import jackknife


def load_runs(paths):
    """Pool runs; returns (meta, O, plaq, poly, W, sizes) with meta = (beta, shape, levels, group)."""
    data = [np.load(p, allow_pickle=True) for p in paths]

    def key(d):
        return (float(d["beta"]), tuple(int(x) for x in d["shape"]), tuple(int(x) for x in d["levels"]),
                str(d["group"]))

    k0 = key(data[0])
    for p, d in zip(paths, data):
        if key(d) != k0:
            raise ValueError(f"{p}: incompatible run {key(d)} vs {k0}")
    O = np.concatenate([d["O"] for d in data])
    plaq = np.concatenate([d["plaq"] for d in data])
    poly = np.concatenate([d["poly"] for d in data])
    maxloop = min(d["W"].shape[1] for d in data)
    W = np.concatenate([d["W"][:, :maxloop, :maxloop] for d in data])
    return k0, O, plaq, poly, W, [len(d["plaq"]) for d in data]


def corr_matrix(O: np.ndarray) -> np.ndarray:
    """O[meas, level, t] -> C[tau, i, j], connected, symmetrised, averaged over t."""
    n, nl, L0 = O.shape
    mean = O.mean(axis=(0, 2))
    C = np.empty((L0, nl, nl))
    for tau in range(L0):
        C[tau] = np.einsum("mit,mjt->ij", O, np.roll(O, -tau, axis=2)) / (n * L0)
    C = 0.5 * (C + np.swapaxes(C, 1, 2))
    return C - np.outer(mean, mean)[None]


def gevp_lambda0(C: np.ndarray, tau0: int = 0) -> np.ndarray:
    """Largest generalised eigenvalue of C(tau) v = lambda C(tau0) v, for every tau (NaN on failure)."""
    out = np.full(C.shape[0], np.nan)
    try:
        for tau in range(C.shape[0]):
            out[tau] = eigh(C[tau], C[tau0], eigvals_only=True)[-1]
    except (np.linalg.LinAlgError, ValueError):
        pass
    return out


def effective_masses(O: np.ndarray, tau0: int = 0, nblocks: int = 20):
    """Dict: level index -> (m_eff, err) for single operators, and 'gevp' -> (m_eff, err).

    Effective masses use the periodic (cosh) form and are returned for tau = 0..L0/2-1."""
    L0 = O.shape[2]
    half = L0 // 2

    def single(s, li):
        return obs.effective_mass(corr_matrix(s[:, li:li + 1, :])[:, 0, 0], periodic=True)[:half]

    res = {li: jackknife(O, lambda s, li=li: single(s, li), nblocks) for li in range(O.shape[1])}
    res["gevp"] = jackknife(O, lambda s: obs.effective_mass(gevp_lambda0(corr_matrix(s), tau0),
                                                            periodic=True)[:half], nblocks)
    return res


def creutz_ratios(W: np.ndarray, nblocks: int = 20) -> dict:
    """Diagonal Creutz ratios chi(R,R) with jackknife errors, from rows where all loops were measured."""
    ok = ~np.isnan(W).any(axis=(1, 2))
    Wm = W[ok]
    out = {}
    if len(Wm) < 2 * nblocks:          # too few measured configurations for a jackknife
        return out, int(ok.sum())
    for R in range(2, Wm.shape[1] + 1):
        def f(s, R=R):
            Wd = {k: np.mean(s[:, k[0] - 1, k[1] - 1]) for k in [(R, R), (R - 1, R - 1), (R, R - 1), (R - 1, R)]}
            return obs.creutz_ratio(Wd, R, R)
        try:
            v, e = jackknife(Wm, f, nblocks)
        except (ValueError, FloatingPointError):
            continue
        if np.isfinite(v) and np.isfinite(e):
            out[R] = (float(v), float(e))
    return out, int(ok.sum())


def best_creutz(chis: dict, max_rel_err: float = 0.1):
    """Largest R whose Creutz ratio has relative error below ``max_rel_err`` (None if none)."""
    good = [R for R, (v, e) in chis.items() if v > 0 and e / v < max_rel_err]
    return max(good) if good else None

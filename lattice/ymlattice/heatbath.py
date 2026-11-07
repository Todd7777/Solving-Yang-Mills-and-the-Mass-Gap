"""Exact SU(2) heat-bath sampling, vectorised.

We need ``a`` in SU(2) distributed as ``exp(k a0) dHaar(a)``, ``k >= 0``.  Writing the Haar
measure as ``(1/pi^2) sqrt(1 - a0^2) da0 dOmega`` (``dOmega`` = solid angle of the
vector part), ``a0`` has density proportional to ``sqrt(1 - a0^2) exp(k a0)`` on ``[-1, 1]``
and the vector part is uniform on the sphere of radius ``sqrt(1 - a0^2)``.

Two samplers are combined:

* Creutz (1980): ``a0 = 1 + log(r + (1 - r) e^{-2k}) / k`` samples ``exp(k a0)`` exactly by
  inversion; accept with probability ``sqrt(1 - a0^2)``.  Efficient for small ``k``.
* Kennedy--Pendleton (1985): with ``lambda^2 = -(log r1 + cos^2(2 pi r2) log r3) / (2k)``,
  accept if ``r4^2 <= 1 - lambda^2`` and set ``a0 = 1 - 2 lambda^2``.  Efficient for large ``k``.
  (``lambda^2`` is Gamma(3/2, rate 2k) distributed, which is exactly what the substitution
  ``a0 = 1 - 2 lambda^2`` requires; see the derivation in the docstring of ``_kp_step``.)

Both are exact; a check of the exact mean ``<a0> = I_2(k)/I_1(k)`` is in the test-suite.
"""

from __future__ import annotations

import numpy as np

KP_THRESHOLD = 2.0  # use Kennedy--Pendleton for k above this, Creutz below


def _creutz_step(k: np.ndarray, rng: np.random.Generator):
    r = rng.random(k.shape)
    small = k < 1e-12
    ksafe = np.where(small, 1.0, k)
    a0 = 1.0 + np.log(r + (1.0 - r) * np.exp(-2.0 * ksafe)) / ksafe
    # k == 0: exp(k a0) is flat, draw a0 uniformly
    a0 = np.where(small, 2.0 * r - 1.0, a0)
    accept = rng.random(k.shape) < np.sqrt(np.clip(1.0 - a0 * a0, 0.0, 1.0))
    return a0, accept


def _kp_step(k: np.ndarray, rng: np.random.Generator):
    """One Kennedy--Pendleton trial.

    Derivation: with ``a0 = 1 - 2 s``, ``sqrt(1 - a0^2) e^{k a0} da0`` becomes
    ``const * sqrt(s) sqrt(1 - s) e^{-2 k s} ds``.  ``sqrt(s) e^{-2ks}`` is the Gamma(3/2, 2k)
    density; Gamma(3/2) = Gamma(1) + Gamma(1/2), and ``-log r1 / (2k)`` is Gamma(1, 2k),
    ``-cos^2(2 pi r2) log r3 / (2k)`` is Gamma(1/2, 2k) (Box--Muller).  The remaining factor
    ``sqrt(1 - s)`` is the acceptance probability.
    """
    r1, r2, r3, r4 = (rng.random(k.shape) for _ in range(4))
    r1 = np.maximum(r1, 1e-300)
    r3 = np.maximum(r3, 1e-300)
    s = -(np.log(r1) + np.cos(2.0 * np.pi * r2) ** 2 * np.log(r3)) / (2.0 * k)
    accept = (s <= 1.0) & (r4 * r4 <= 1.0 - s)
    a0 = 1.0 - 2.0 * s
    return a0, accept


def sample_a0(k: np.ndarray, rng: np.random.Generator, max_iter: int = 10000) -> np.ndarray:
    """Vectorised sampling of ``a0`` with density ``~ sqrt(1 - a0^2) exp(k a0)``."""
    k = np.asarray(k, dtype=float)
    shape = k.shape
    k = k.ravel()
    out = np.empty(k.shape)
    todo = np.arange(k.size)
    use_kp = k > KP_THRESHOLD
    for _ in range(max_iter):
        if todo.size == 0:
            break
        kk = k[todo]
        kp = use_kp[todo]
        a0 = np.empty(kk.shape)
        acc = np.zeros(kk.shape, dtype=bool)
        if np.any(kp):
            a0[kp], acc[kp] = _kp_step(kk[kp], rng)
        if np.any(~kp):
            a0[~kp], acc[~kp] = _creutz_step(kk[~kp], rng)
        out[todo[acc]] = a0[acc]
        todo = todo[~acc]
    else:  # pragma: no cover
        raise RuntimeError("heat bath failed to converge")
    return out.reshape(shape)


def sample_su2(k: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Draw SU(2) elements ``a`` with density ``exp(k a0) dHaar``, ``k`` array of shape ``(n,)``."""
    k = np.asarray(k, dtype=float)
    a0 = sample_a0(k, rng)
    rad = np.sqrt(np.clip(1.0 - a0 * a0, 0.0, 1.0))
    cos_t = 2.0 * rng.random(k.shape) - 1.0
    sin_t = np.sqrt(np.clip(1.0 - cos_t * cos_t, 0.0, 1.0))
    phi = 2.0 * np.pi * rng.random(k.shape)
    a = np.empty(k.shape + (4,))
    a[..., 0] = a0
    a[..., 1] = rad * sin_t * np.cos(phi)
    a[..., 2] = rad * sin_t * np.sin(phi)
    a[..., 3] = rad * cos_t
    return a

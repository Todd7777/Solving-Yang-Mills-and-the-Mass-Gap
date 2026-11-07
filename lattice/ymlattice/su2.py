"""SU(2) as unit quaternions, vectorised with NumPy.

An SU(2) matrix is written ``U = a0 * 1 + i (a1 s1 + a2 s2 + a3 s3)`` with Pauli
matrices ``s_k`` and ``a0^2 + a1^2 + a2^2 + a3^2 = 1``.  Explicitly

    U = [[a0 + i a3,  a2 + i a1],
         [-a2 + i a1, a0 - i a3]].

Arrays have shape ``(..., 4)``.  All functions broadcast over leading axes.
Sums of SU(2) matrices are again proportional to SU(2) matrices, which is what
makes the SU(2) heat bath exact; the staple sum ``A`` satisfies ``A = xi * Ahat``
with ``xi = sqrt(det A) = |A|`` (Euclidean 4-norm) and ``Ahat`` in SU(2).
"""

from __future__ import annotations

import numpy as np


def mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Quaternion (= SU(2) matrix) product ``a @ b``.

    For ``a = a0 + i a.s`` and ``b = b0 + i b.s``:
    ``ab = (a0 b0 - a.b) + i (a0 b + b0 a - a x b).s``.
    """
    a0, av = a[..., :1], a[..., 1:]
    b0, bv = b[..., :1], b[..., 1:]
    s = a0 * b0 - np.sum(av * bv, axis=-1, keepdims=True)
    v = a0 * bv + b0 * av - np.cross(av, bv)
    return np.concatenate([s, v], axis=-1)


def dagger(a: np.ndarray) -> np.ndarray:
    """Hermitian conjugate = quaternion conjugate ``(a0, -a)``."""
    out = -a
    out[..., 0] = a[..., 0]
    return out


def mul_dag(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """``a @ b^dagger``."""
    return mul(a, dagger(b))


def dag_mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """``a^dagger @ b``."""
    return mul(dagger(a), b)


def retrace(a: np.ndarray) -> np.ndarray:
    """``(1/2) Re Tr a = a0`` (the normalised real trace)."""
    return a[..., 0]


def norm(a: np.ndarray) -> np.ndarray:
    """Euclidean 4-norm ``|a| = sqrt(det a)`` for a general (non-unit) quaternion."""
    return np.sqrt(np.sum(a * a, axis=-1))


def normalize(a: np.ndarray) -> np.ndarray:
    """Project a non-zero quaternion onto SU(2)."""
    return a / norm(a)[..., None]


def identity(shape=()) -> np.ndarray:
    out = np.zeros(tuple(shape) + (4,))
    out[..., 0] = 1.0
    return out


def random(shape, rng: np.random.Generator) -> np.ndarray:
    """Haar-random SU(2) elements: normalised 4-d Gaussian vectors."""
    g = rng.standard_normal(tuple(shape) + (4,))
    return normalize(g)


def to_matrix(a: np.ndarray) -> np.ndarray:
    """Convert quaternions ``(..., 4)`` to complex matrices ``(..., 2, 2)``."""
    a0, a1, a2, a3 = (a[..., k] for k in range(4))
    m = np.empty(a.shape[:-1] + (2, 2), dtype=complex)
    m[..., 0, 0] = a0 + 1j * a3
    m[..., 0, 1] = a2 + 1j * a1
    m[..., 1, 0] = -a2 + 1j * a1
    m[..., 1, 1] = a0 - 1j * a3
    return m


def from_matrix(m: np.ndarray) -> np.ndarray:
    """Inverse of :func:`to_matrix` (assumes ``m`` is in SU(2))."""
    a = np.empty(m.shape[:-2] + (4,))
    a[..., 0] = m[..., 0, 0].real
    a[..., 3] = m[..., 0, 0].imag
    a[..., 2] = m[..., 0, 1].real
    a[..., 1] = m[..., 0, 1].imag
    return a


def project_from_complex_2x2(w: np.ndarray) -> np.ndarray:
    """The SU(2) *projection* of an arbitrary complex ``2x2`` matrix ``w``.

    Returns the real 4-vector ``r`` such that ``Re Tr(a w) = sum_mu a_mu r_mu`` for every
    quaternion ``a``.  Writing ``r = |r| Rhat`` with ``Rhat`` in SU(2) one has
    ``Re Tr(a w) = (|r|/2) Tr(a Rhat^dagger) = |r| (a Rhat^dagger)_0``.
    This is the identity behind the Cabibbo--Marinari update of SU(N), N >= 3.
    """
    w00, w01, w10, w11 = w[..., 0, 0], w[..., 0, 1], w[..., 1, 0], w[..., 1, 1]
    r = np.empty(w.shape[:-2] + (4,))
    r[..., 0] = (w00 + w11).real
    r[..., 1] = -(w01 + w10).imag
    r[..., 2] = (w10 - w01).real
    r[..., 3] = -(w00 - w11).imag
    return r

"""SU(3) helpers: complex 3x3 arrays of shape ``(..., 3, 3)``.

Includes the embedding of SU(2) subgroups used by the Cabibbo--Marinari update,
reunitarisation, Haar-random elements and a polar projection used for APE smearing.
"""

from __future__ import annotations

import numpy as np

from . import su2

# The three SU(2) subgroups (row/column index pairs) used by Cabibbo--Marinari.
SUBGROUPS = ((0, 1), (0, 2), (1, 2))


def dagger(m: np.ndarray) -> np.ndarray:
    return np.conj(np.swapaxes(m, -1, -2))


def mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.matmul(a, b)


def mul_dag(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.matmul(a, dagger(b))


def dag_mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.matmul(dagger(a), b)


def retrace(m: np.ndarray) -> np.ndarray:
    """``(1/3) Re Tr m``."""
    return np.real(m[..., 0, 0] + m[..., 1, 1] + m[..., 2, 2]) / 3.0


def identity(shape=()) -> np.ndarray:
    out = np.zeros(tuple(shape) + (3, 3), dtype=complex)
    out[..., 0, 0] = out[..., 1, 1] = out[..., 2, 2] = 1.0
    return out


def reunitarize(m: np.ndarray) -> np.ndarray:
    """Gram--Schmidt projection back to SU(3): rows 0,1 orthonormalised, row 2 = conj(row0 x row1)."""
    r0 = m[..., 0, :]
    r0 = r0 / np.linalg.norm(r0, axis=-1, keepdims=True)
    r1 = m[..., 1, :]
    r1 = r1 - np.sum(np.conj(r0) * r1, axis=-1, keepdims=True) * r0
    r1 = r1 / np.linalg.norm(r1, axis=-1, keepdims=True)
    r2 = np.conj(np.cross(r0, r1))
    return np.stack([r0, r1, r2], axis=-2)


def random(shape, rng: np.random.Generator) -> np.ndarray:
    """Haar-random SU(3) via QR of a complex Ginibre matrix, phases fixed, det normalised."""
    shape = tuple(shape)
    z = rng.standard_normal(shape + (3, 3)) + 1j * rng.standard_normal(shape + (3, 3))
    q, r = np.linalg.qr(z)
    d = np.diagonal(r, axis1=-2, axis2=-1)
    ph = d / np.abs(d)
    q = q * ph[..., None, :]
    det = np.linalg.det(q)
    # divide by a cube root of det to land in SU(3)
    q = q / (det ** (1.0 / 3.0))[..., None, None]
    return q


def embed_su2(a: np.ndarray, sub: tuple[int, int]) -> np.ndarray:
    """Embed quaternions ``a`` ``(..., 4)`` into SU(3) in the subgroup ``sub``."""
    i, j = sub
    m = identity(a.shape[:-1])
    blk = su2.to_matrix(a)
    m[..., i, i] = blk[..., 0, 0]
    m[..., i, j] = blk[..., 0, 1]
    m[..., j, i] = blk[..., 1, 0]
    m[..., j, j] = blk[..., 1, 1]
    return m


def extract_block(w: np.ndarray, sub: tuple[int, int]) -> np.ndarray:
    """The complex 2x2 block of ``w`` at rows/cols ``sub``."""
    i, j = sub
    blk = np.empty(w.shape[:-2] + (2, 2), dtype=complex)
    blk[..., 0, 0] = w[..., i, i]
    blk[..., 0, 1] = w[..., i, j]
    blk[..., 1, 0] = w[..., j, i]
    blk[..., 1, 1] = w[..., j, j]
    return blk


def left_mul_su2(a: np.ndarray, sub: tuple[int, int], m: np.ndarray) -> np.ndarray:
    """Efficiently compute ``embed(a, sub) @ m`` (only two rows of ``m`` change)."""
    i, j = sub
    blk = su2.to_matrix(a)
    out = m.copy()
    out[..., i, :] = blk[..., 0, 0, None] * m[..., i, :] + blk[..., 0, 1, None] * m[..., j, :]
    out[..., j, :] = blk[..., 1, 0, None] * m[..., i, :] + blk[..., 1, 1, None] * m[..., j, :]
    return out


def project_su3(m: np.ndarray) -> np.ndarray:
    """Polar projection of an arbitrary complex 3x3 matrix onto SU(3).

    ``m = V S W^dagger`` (SVD); the closest unitary is ``V W^dagger``; the phase is removed by
    dividing by a cube root of its determinant.  Used for APE smearing.
    """
    v, _, wh = np.linalg.svd(m)
    u = np.matmul(v, wh)
    det = np.linalg.det(u)
    return u / (det ** (1.0 / 3.0))[..., None, None]

"""Numerical checks of the Dobrushin estimates of Part II (theory/infinite-volume-strong-coupling.tex).

Lemma iv:lem:tvratio: two probability densities whose ratio lies in [c e^-delta, c e^delta] are at total-variation
distance at most tanh(delta/2), and the bound is attained on a two-point space.

Lemma iv:lem:coeff: for the Wilson weight w = exp((beta/N) Re Tr U), the single-link conditional distributions of two
configurations that differ at one adjacent link are at total-variation distance at most tanh(beta).  For SU(2) the
conditional density of a link is proportional to exp(beta <U, B>) on the unit sphere S^3 (quaternions), where B is the
sum of the 2(d-1) staples (each a unit quaternion) and <.,.> is the Euclidean inner product; changing one adjacent link
replaces one staple.  Both densities depend on U only through its projection onto the plane spanned by the two staple
sums, and the uniform measure on S^3 projects to the uniform measure on the unit disc (Archimedes), so the distance is
a two-dimensional integral, computed here by quadrature.
"""
import numpy as np

from ymlattice import su2


def _disc_quadrature(n_r=400, n_phi=400):
    r = (np.arange(n_r) + 0.5) / n_r
    phi = 2 * np.pi * (np.arange(n_phi) + 0.5) / n_phi
    R, PHI = np.meshgrid(r, phi, indexing="ij")
    pts = np.stack([R * np.cos(PHI), R * np.sin(PHI)], axis=-1).reshape(-1, 2)
    wts = (R / (n_r * n_phi)).reshape(-1) * 2.0  # uniform probability measure on the unit disc
    return pts, wts


PTS, WTS = _disc_quadrature()


def _tv_staples(beta, B, Bp):
    """Total variation between the laws prop. to exp(beta<U,B>) and exp(beta<U,Bp>) on S^3."""
    e1 = B / np.linalg.norm(B) if np.linalg.norm(B) > 0 else np.array([1.0, 0, 0, 0])
    f = Bp - (Bp @ e1) * e1
    e2 = f / np.linalg.norm(f) if np.linalg.norm(f) > 1e-14 else np.array([0, 1.0, 0, 0]) - e1[1] * e1
    e2 = e2 / np.linalg.norm(e2)
    b, bp = np.array([B @ e1, B @ e2]), np.array([Bp @ e1, Bp @ e2])
    p = np.exp(beta * PTS @ b)
    q = np.exp(beta * PTS @ bp)
    p /= np.sum(p * WTS)
    q /= np.sum(q * WTS)
    return 0.5 * np.sum(np.abs(p - q) * WTS)


def test_disc_projection_of_haar_measure():
    """Sanity check of the quadrature: E[exp(beta a0)] for Haar SU(2) equals 2 I_1(beta)/beta."""
    from scipy.special import iv
    for beta in (0.5, 2.0, 5.0):
        val = np.sum(np.exp(beta * PTS[:, 0]) * WTS)
        assert abs(val - 2 * iv(1, beta) / beta) < 1e-4 * val


def test_tv_ratio_lemma_two_point_example_is_sharp():
    for delta in (0.1, 0.7, 2.5):
        pi = np.array([1.0, np.exp(delta)]) / (1 + np.exp(delta))
        h = np.array([np.exp(delta), np.exp(-delta)])
        assert abs(np.sum(h * pi) - 1) < 1e-12
        tv = 0.5 * np.sum(np.abs(h * pi - pi))
        assert abs(tv - np.tanh(delta / 2)) < 1e-12


def test_tv_ratio_lemma_random_examples():
    rng = np.random.default_rng(11)
    for _ in range(2000):
        n = rng.integers(2, 8)
        pi = rng.dirichlet(np.ones(n))
        delta = rng.uniform(0.01, 3.0)
        h = np.exp(rng.uniform(-delta, delta, size=n))
        h /= np.sum(h * pi)                      # ratio now in [c e^-delta, c e^delta] for some c > 0
        tv = 0.5 * np.sum(np.abs(h * pi - pi))
        assert tv <= np.tanh(delta / 2) + 1e-12


def test_single_link_dobrushin_coefficient_su2():
    """Lemma iv:lem:coeff for SU(2) in d = 4: TV <= tanh(beta) for one changed staple, incl. adversarial choices."""
    rng = np.random.default_rng(5)
    d = 4
    for beta in (0.02, 0.05, 0.2, 1.0, 2.5):
        worst = 0.0
        for trial in range(40):
            staples = su2.random((2 * (d - 1),), rng)
            if trial % 2 == 0:
                staples[1:] = staples[0]          # all staples aligned: largest field
            new = su2.random((), rng)
            if trial % 3 == 0:
                new = -staples[0]                 # reverse the changed staple
            B = staples.sum(axis=0)
            Bp = B - staples[0] + new
            worst = max(worst, _tv_staples(beta, B, Bp))
        assert worst <= np.tanh(beta) + 1e-6, (beta, worst)
        # the cruder bound of the foundations monograph is much weaker
        assert np.tanh(beta) < np.exp(4 * beta) - 1

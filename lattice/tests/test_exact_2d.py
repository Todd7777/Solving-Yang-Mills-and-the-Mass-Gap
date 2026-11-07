"""Two-dimensional torus: the character expansion is exact, so the Monte Carlo must reproduce
<P> = sum_j chat_j^{V-1} chat_j' / sum_j chat_j^V (SU(2)) and the one-link integral (SU(3))."""
import numpy as np
import pytest

from ymlattice import exact
from ymlattice.lattice import Lattice
from ymlattice.stats import binned_mean_err


def _run(group, beta, L=12, ntherm=200, nmeas=1500, n_or=2):
    lat = Lattice((L, L), beta=beta, group=group, seed=int(beta * 100) + L, start="hot")
    lat.thermalize(ntherm, n_overrelax=n_or)
    vals = []
    for _ in range(nmeas):
        lat.update(1, n_or)
        vals.append(lat.plaquette())
    return binned_mean_err(np.array(vals), 25)


@pytest.mark.parametrize("beta", [0.8, 2.5, 6.0])
def test_su2_two_dimensional_torus(beta):
    L = 12
    mean, err = _run("su2", beta, L)
    ex = exact.su2_plaquette_2d(beta, V=L * L)
    assert abs(mean - ex) < max(4 * err, 2.5e-3), (beta, mean, err, ex)


@pytest.mark.parametrize("beta", [1.5, 4.0, 9.0])
def test_su3_two_dimensional_torus(beta):
    L = 12
    mean, err = _run("su3", beta, L, nmeas=800)
    ex = exact.su3_plaquette_2d(beta)
    assert abs(mean - ex) < max(4 * err, 2.5e-3), (beta, mean, err, ex)


def test_exact_2d_limits():
    # infinite-volume limit of the torus formula is I_2/I_1; weak and strong coupling limits
    for beta in (0.5, 3.0, 10.0):
        assert np.isclose(exact.su2_plaquette_2d(beta, V=100000), exact.su2_u(beta), atol=1e-9)
    assert np.isclose(exact.su2_u(0.01), 0.01 / 4, rtol=1e-3)
    assert np.isclose(exact.su2_u(200.0), exact.plaquette_weak(200.0, 2, d=2), atol=2e-4)
    assert np.isclose(exact.su3_u(0.01), 0.01 / 18, rtol=2e-3)
    assert np.isclose(exact.su3_u(400.0), exact.plaquette_weak(400.0, 3, d=2), atol=3e-4)


def test_strong_coupling_nlo_reduces_to_leading():
    for N in (2, 3):
        assert np.isclose(exact.plaquette_strong_nlo(0.05, N), exact.plaquette_strong_leading(0.05, N), atol=1e-9)
        assert exact.plaquette_strong_nlo(1.0, N) > exact.plaquette_strong_leading(1.0, N)


@pytest.mark.parametrize("shape", [(5, 3), (7, 4)])
def test_su2_two_dimensional_odd_torus(shape):
    """Odd periodic lengths need the 3-colouring update (plain checkerboard is invalid there)."""
    beta = 2.5
    lat = Lattice(shape, beta=beta, group="su2", seed=sum(shape), start="hot")
    lat.thermalize(200, n_overrelax=1)
    vals = []
    for _ in range(4000):
        lat.update(1, 1)
        vals.append(lat.plaquette())
    mean, err = binned_mean_err(np.array(vals), 25)
    naive = np.std(vals, ddof=1) / np.sqrt(len(vals))
    ex = exact.su2_plaquette_2d(beta, V=shape[0] * shape[1])
    assert abs(mean - ex) < 4 * max(err, naive), (shape, mean, err, ex)

"""Four-dimensional checks against the strong-coupling expansion, the one-loop weak-coupling
expansion and widely quoted reference plaquette values.  Slow: run with ``pytest --runslow``."""
import numpy as np
import pytest

from ymlattice import exact
from ymlattice.lattice import Lattice
from ymlattice.stats import binned_mean_err


def _run(group, beta, L=6, ntherm=100, nmeas=300, n_or=2, seed=11):
    lat = Lattice((L, L, L, L), beta=beta, group=group, seed=seed, start="hot")
    lat.thermalize(ntherm, n_overrelax=n_or)
    vals = []
    for _ in range(nmeas):
        lat.update(1, n_or)
        vals.append(lat.plaquette())
    return binned_mean_err(np.array(vals), 20)


@pytest.mark.slow
def test_su2_strong_coupling_4d():
    beta = 0.5
    mean, err = _run("su2", beta, nmeas=200)
    u = exact.su2_u(beta)
    # <P> = u + 4(d-2) u^4 + ... ; at beta = 0.5, u ~ 0.124, the correction is ~ 2e-3
    assert abs(mean - u) < 4e-3, (mean, err, u)


@pytest.mark.slow
def test_su2_weak_coupling_4d():
    beta = 12.0
    mean, err = _run("su2", beta, nmeas=200)
    assert abs(mean - exact.plaquette_weak(beta, 2)) < 6e-3, (mean, err)


@pytest.mark.slow
def test_su3_weak_coupling_4d():
    beta = 30.0
    mean, err = _run("su3", beta, L=4, nmeas=200)
    assert abs(mean - exact.plaquette_weak(beta, 3)) < 6e-3, (mean, err)


@pytest.mark.slow
def test_su3_reference_plaquette_beta6():
    """<P>(beta = 6.0) = 0.5937 is one of the most frequently quoted numbers in lattice QCD."""
    mean, err = _run("su3", 6.0, L=6, ntherm=150, nmeas=300)
    assert abs(mean - 0.5937) < 0.006, (mean, err)

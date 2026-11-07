import numpy as np
from scipy.special import ive

from ymlattice import su2
from ymlattice.heatbath import sample_su2
from ymlattice.lattice import Lattice
from ymlattice import observables as obs


def test_quaternion_product_matches_matrix_product():
    rng = np.random.default_rng(1)
    a, b = su2.random((50,), rng), su2.random((50,), rng)
    lhs = su2.to_matrix(su2.mul(a, b))
    rhs = su2.to_matrix(a) @ su2.to_matrix(b)
    assert np.allclose(lhs, rhs)
    assert np.allclose(su2.to_matrix(su2.dagger(a)), np.conj(np.swapaxes(su2.to_matrix(a), -1, -2)))
    assert np.allclose(su2.from_matrix(su2.to_matrix(a)), a)
    m = su2.to_matrix(a)
    assert np.allclose(np.linalg.det(m), 1.0)
    assert np.allclose(su2.retrace(a), np.trace(m, axis1=-2, axis2=-1).real / 2)


def test_su2_projection_identity():
    rng = np.random.default_rng(2)
    w = rng.standard_normal((40, 2, 2)) + 1j * rng.standard_normal((40, 2, 2))
    a = su2.random((40,), rng)
    r = su2.project_from_complex_2x2(w)
    lhs = np.trace(su2.to_matrix(a) @ w, axis1=-2, axis2=-1).real
    assert np.allclose(lhs, np.sum(a * r, axis=-1))


def test_heatbath_mean_is_bessel_ratio():
    """<a0> under sqrt(1-a0^2) exp(k a0) equals I_2(k)/I_1(k) (exercises both Creutz and KP branches)."""
    rng = np.random.default_rng(3)
    n = 200_000
    for k in (0.0, 0.3, 1.9, 2.1, 6.0, 40.0):
        a = sample_su2(np.full(n, k), rng)
        assert np.allclose(np.sum(a * a, axis=-1), 1.0)
        mean = a[:, 0].mean()
        exact = ive(2, k) / ive(1, k) if k > 0 else 0.0
        err = a[:, 0].std() / np.sqrt(n)
        assert abs(mean - exact) < 4 * err + 1e-6, (k, mean, exact, err)
        # vector part isotropic
        assert abs(a[:, 1].mean()) < 5 * err and abs(a[:, 3].mean()) < 5 * err


def test_staple_sum_reproduces_plaquette_action():
    """sum_links Re Tr(U_l A_l) = 4 * sum_p Re Tr U_p (every plaquette contains four links)."""
    lat = Lattice((4, 4, 4, 4), beta=1.0, group="su2", seed=4, start="hot")
    tot = 0.0
    for mu in range(lat.d):
        tot += np.sum(su2.retrace(su2.mul(lat.U[mu], lat.staple_sum(mu))))
    assert np.isclose(tot, 4 * lat.nplaq * lat.plaquette())


def test_gauge_invariance_of_observables():
    lat = Lattice((4, 4, 4, 6), beta=2.0, group="su2", seed=5, start="hot")
    lat.thermalize(3)
    before = (lat.plaquette(), obs.wilson_loop(lat, 2, 3), obs.polyakov_loop(lat),
              obs.spatial_plaquette_timeslices(lat).copy())
    lat.gauge_transform()
    after = (lat.plaquette(), obs.wilson_loop(lat, 2, 3), obs.polyakov_loop(lat),
             obs.spatial_plaquette_timeslices(lat).copy())
    assert np.isclose(before[0], after[0])
    assert np.isclose(before[1], after[1])
    assert np.isclose(before[2], after[2])
    assert np.allclose(before[3], after[3])
    # APE smearing is gauge covariant: smeared observables are gauge invariant too
    Us = obs.ape_smear(lat, 0.5, 2)
    s_after = obs.spatial_plaquette_timeslices(lat, Us)
    lat.gauge_transform()
    Us2 = obs.ape_smear(lat, 0.5, 2)
    assert np.allclose(s_after, obs.spatial_plaquette_timeslices(lat, Us2))


def test_overrelaxation_preserves_action_and_group():
    lat = Lattice((4, 4, 4, 4), beta=2.3, group="su2", seed=6, start="hot")
    lat.thermalize(2)
    s0 = lat.action()
    lat.sweep("overrelax")
    assert np.isclose(lat.action(), s0, rtol=1e-9)
    assert np.allclose(np.sum(lat.U * lat.U, axis=-1), 1.0)


def test_wilson_loop_1x1_is_plaquette():
    lat = Lattice((4, 4, 4, 4), beta=2.3, group="su2", seed=7, start="hot")
    lat.thermalize(2)
    assert np.isclose(obs.wilson_loop(lat, 1, 1), lat.plaquette())

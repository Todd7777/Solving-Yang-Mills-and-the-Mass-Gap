import numpy as np

from ymlattice import su2, su3, exact
from ymlattice.lattice import Lattice, cabibbo_marinari
from ymlattice import observables as obs


def _is_su3(m, tol=1e-10):
    eye = np.eye(3)
    return np.allclose(m @ su3.dagger(m), eye, atol=tol) and np.allclose(np.linalg.det(m), 1.0, atol=tol)


def test_random_and_reunitarize_and_project():
    rng = np.random.default_rng(1)
    m = su3.random((30,), rng)
    assert _is_su3(m)
    noisy = m + 1e-3 * (rng.standard_normal(m.shape) + 1j * rng.standard_normal(m.shape))
    assert _is_su3(su3.reunitarize(noisy))
    assert _is_su3(su3.project_su3(noisy))
    assert np.allclose(su3.project_su3(m), m, atol=1e-8)


def test_embedding_and_left_mul():
    rng = np.random.default_rng(2)
    a = su2.random((10,), rng)
    m = su3.random((10,), rng)
    for sub in su3.SUBGROUPS:
        e = su3.embed_su2(a, sub)
        assert _is_su3(e)
        assert np.allclose(su3.left_mul_su2(a, sub, m), e @ m)


def test_single_link_heatbath_matches_weyl_quadrature():
    """Links with staple sum A = xi * 1 have weight exp((beta xi/3) Re Tr U): <P> is a one-link integral."""
    rng = np.random.default_rng(3)
    n = 4000
    for beta, xi in ((1.5, 1.0), (6.0, 1.0), (3.0, 2.0)):
        U = su3.random((n,), rng)
        A = xi * su3.identity((n,))
        vals = []
        for it in range(60):
            U = cabibbo_marinari(U, A, beta, rng, "heatbath", nhits=1)
            if it >= 20:
                vals.append(su3.retrace(U))
        vals = np.array(vals)
        mean = vals.mean()
        err = vals.mean(axis=1).std(ddof=1) / np.sqrt(len(vals))  # successive sweeps are nearly independent
        exact_p = exact.su3_one_link(beta * xi)[1]
        assert abs(mean - exact_p) < max(5 * err, 2e-3), (beta, xi, mean, exact_p, err)
        assert _is_su3(U, tol=1e-8)


def test_overrelaxation_preserves_action_su3():
    lat = Lattice((4, 4, 4, 4), beta=5.7, group="su3", seed=4, start="hot")
    lat.thermalize(2)
    s0 = lat.action()
    lat.sweep("overrelax")
    assert np.isclose(lat.action(), s0, rtol=1e-8)
    assert _is_su3(lat.U.reshape(-1, 3, 3), tol=1e-8)


def test_gauge_invariance_su3():
    lat = Lattice((4, 4, 4, 4), beta=5.7, group="su3", seed=5, start="hot")
    lat.thermalize(2)
    b = (lat.plaquette(), obs.wilson_loop(lat, 2, 2), obs.polyakov_loop(lat))
    lat.gauge_transform()
    a = (lat.plaquette(), obs.wilson_loop(lat, 2, 2), obs.polyakov_loop(lat))
    assert np.allclose(b, a)


def test_staple_sum_reproduces_plaquette_action_su3():
    lat = Lattice((4, 4, 4, 4), beta=1.0, group="su3", seed=6, start="hot")
    tot = 0.0
    for mu in range(lat.d):
        tot += np.sum(su3.retrace(lat.U[mu] @ lat.staple_sum(mu)))
    assert np.isclose(tot, 4 * lat.nplaq * lat.plaquette())

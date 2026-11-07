"""End-to-end test of the gap-measurement pipeline against the exactly solvable 2-d transfer matrix.

Monte Carlo of SU(2) on a periodic L0 x L lattice; the time-slice operator is half the trace of the
spatial Polyakov loop; its correlator must equal the exact transfer-matrix prediction of
ymlattice.transfer2d at every tau, and the effective mass must reproduce the exact gap L*log(I1/I2).
The odd spatial length L = 3 exercises the 3-colouring update (a plain checkerboard is invalid
there and was caught by this test).
"""
import numpy as np
import pytest

from ymlattice import Lattice, su2
from ymlattice import observables as obs
from ymlattice import transfer2d
from ymlattice.stats import jackknife


def test_exact_transfer_matrix_identities():
    beta, L, L0 = 2.0, 2, 12
    lam = transfer2d.eigenvalue_ratios(beta, L)
    assert lam[0] == 1.0 and np.all(np.diff(lam) < 0)
    assert np.isclose(-np.log(lam[1]), transfer2d.exact_gap(beta, L))
    C = transfer2d.exact_correlator(beta, L, L0)
    assert np.allclose(C, C[(-np.arange(L0)) % L0])           # periodic symmetry tau -> L0 - tau
    # large L0: pure exponential with the exact gap
    Cbig = transfer2d.exact_correlator(beta, L, 200)
    assert np.isclose(np.log(Cbig[3] / Cbig[4]), transfer2d.exact_gap(beta, L), rtol=1e-6)
    # <A^2> at tau=0 and large L0 is the vacuum expectation <(chi_{1/2}/2)^2> = 1/4
    assert np.isclose(Cbig[0], 0.25, rtol=1e-6)


@pytest.mark.parametrize("beta,L", [(2.0, 2), (4.0, 3)])
def test_monte_carlo_matches_exact_transfer_matrix(beta, L):
    L0 = 12
    lat = Lattice((L0, L), beta, "su2", seed=int(10 * beta) + L, start="hot")
    lat.thermalize(200, n_overrelax=1)
    nmeas = 12000
    O = np.empty((nmeas, L0))
    for n in range(nmeas):
        lat.update(1, 1)
        hol = obs.line(lat, lat.U, 1, L)          # spatial holonomy starting at every site
        O[n] = su2.retrace(hol[:, 0])            # (1/2) Tr at spatial site 0, all times
    corr = lambda s: np.array([np.mean(s * np.roll(s, -t, axis=1)) for t in range(L0)])
    C, Cerr = jackknife(O, corr, 20)
    # A 20-block jackknife error is itself uncertain by ~16%; guard against an unlucky underestimate
    # with the error computed as if the measurements were independent (no autocorrelation here).
    per_meas = np.array([np.mean(O * np.roll(O, -t, axis=1), axis=1) for t in range(L0)]).T
    naive = per_meas.std(axis=0, ddof=1) / np.sqrt(nmeas)
    err = np.maximum(Cerr, naive)
    exact = transfer2d.exact_correlator(beta, L, L0)
    for tau in range(4):
        assert abs(C[tau] - exact[tau]) < 4 * err[tau], (tau, C[tau], err[tau], exact[tau])
    m, merr = jackknife(O, lambda s: obs.effective_mass(
        np.array([np.mean(s * np.roll(s, -t, axis=1)) for t in range(L0)]), periodic=True), 20)
    m_exact = obs.effective_mass(exact, periodic=True)
    assert abs(m[0] - m_exact[0]) < 4 * merr[0] + 1e-3, (m[0], merr[0], m_exact[0])
    # at L0 = 12 the cosh effective mass already equals the exact gap to high accuracy
    assert np.isclose(m_exact[1], transfer2d.exact_gap(beta, L), rtol=1e-3)

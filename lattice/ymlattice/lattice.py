"""The lattice, the Wilson measure, and the Markov-chain updates.

The configuration is stored as ``self.U`` with shape

* SU(2): ``(d, L_0, ..., L_{d-1}, 4)``   (unit quaternions)
* SU(3): ``(d, L_0, ..., L_{d-1}, 3, 3)`` (complex matrices)

Periodic boundary conditions in every direction.  ``shift(f, mu, +1)[x] = f[x + mu]``.

Updates are done with a checkerboard decomposition: for a fixed direction ``mu`` all links of
even (then odd) site parity are updated simultaneously.  The staple of the link ``(x, mu)``
contains ``mu``-links only at the sites ``x +- nu``, which have the opposite parity, and
``nu``-links which are not being updated, so the simultaneous update is a valid
heat-bath / overrelaxation step for the conditional distributions.
"""

from __future__ import annotations

import numpy as np

from . import su2, su3
from .heatbath import sample_su2


def _parity_mask(shape) -> np.ndarray:
    grids = np.meshgrid(*[np.arange(L) for L in shape], indexing="ij")
    return (sum(grids) % 2).astype(bool)  # True = odd


def cabibbo_marinari(U: np.ndarray, A: np.ndarray, beta: float, rng: np.random.Generator,
                     kind: str = "heatbath", nhits: int = 1) -> np.ndarray:
    """Cabibbo--Marinari update of SU(3) links ``U`` ``(n,3,3)`` with staple sums ``A`` ``(n,3,3)``.

    The conditional weight of a link is ``exp((beta/3) Re Tr(U A))``.  For each SU(2) subgroup we
    write ``U -> a U`` and ``Re Tr(a U A) = Re Tr(a W)``, ``W = U A``; the 2x2 block ``w`` of ``W`` gives
    ``Re Tr(a w) = a . r = |r| (a Rhat^+)_0`` so that ``a' = a Rhat^+`` is an SU(2) heat-bath variable
    with ``k = beta |r| / 3``.  Overrelaxation uses ``a = Rhat^2`` (the measure-preserving involution
    ``a' -> a'^+``)."""
    W = np.matmul(U, A)
    for _ in range(nhits):
        for sub in su3.SUBGROUPS:
            w = su3.extract_block(W, sub)
            r = su2.project_from_complex_2x2(w)
            rn = su2.norm(r)
            Rhat = r / np.maximum(rn, 1e-300)[:, None]
            if kind == "heatbath":
                k = beta * rn / 3.0
                ap = sample_su2(k, rng)
                a = su2.mul(ap, Rhat)
            elif kind == "overrelax":
                a = su2.mul(Rhat, Rhat)
            else:
                raise ValueError(kind)
            U = su3.left_mul_su2(a, sub, U)
            W = su3.left_mul_su2(a, sub, W)
    return U


class Lattice:
    def __init__(self, shape, beta: float, group: str = "su2", seed: int | None = None,
                 start: str = "cold"):
        self.shape = tuple(int(L) for L in shape)
        self.d = len(self.shape)
        self.beta = float(beta)
        self.group = group.lower()
        if self.group not in ("su2", "su3"):
            raise ValueError("group must be 'su2' or 'su3'")
        self.N = 2 if self.group == "su2" else 3
        self.rng = np.random.default_rng(seed)
        self.odd = _parity_mask(self.shape)
        self.even = ~self.odd
        self.nsites = int(np.prod(self.shape))
        self.nplaq = self.nsites * self.d * (self.d - 1) // 2
        self.G = su2 if self.group == "su2" else su3
        if start == "cold":
            self.U = self.G.identity((self.d,) + self.shape)
        elif start == "hot":
            self.U = self.G.random((self.d,) + self.shape, self.rng)
        else:
            raise ValueError("start must be 'cold' or 'hot'")

    # ------------------------------------------------------------------ geometry
    def shift(self, f: np.ndarray, mu: int, sign: int) -> np.ndarray:
        """``shift(f, mu, +1)[x] = f[x + mu]``; ``shift(f, mu, -1)[x] = f[x - mu]`` (periodic).

        ``f`` has the lattice shape in its first ``d`` axes."""
        return np.roll(f, -sign, axis=mu)

    # ------------------------------------------------------------------ plaquettes / staples
    def plaquette_field(self, mu: int, nu: int) -> np.ndarray:
        """``U_p(x) = U_mu(x) U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+`` at every site."""
        G, U = self.G, self.U
        a = G.mul(U[mu], self.shift(U[nu], mu, +1))
        b = G.mul(U[nu], self.shift(U[mu], nu, +1))  # U_nu(x) U_mu(x+nu); U_p = a b^+ = U_mu(x) U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+
        return G.mul_dag(a, b)

    def staple_sum(self, mu: int) -> np.ndarray:
        """``A(x) = sum_{nu != mu} [ U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+ + U_nu(x+mu-nu)^+ U_mu(x-nu)^+ U_nu(x-nu) ]``

        so that ``sum_{p ∋ (x,mu)} Tr U_p = Tr( U_mu(x) A(x) )`` (plaquettes containing the link,
        all oriented so that the link enters as ``U_mu(x)``)."""
        G, U = self.G, self.U
        A = None
        for nu in range(self.d):
            if nu == mu:
                continue
            Unu_xmu = self.shift(U[nu], mu, +1)                      # U_nu(x+mu)
            Umu_xnu = self.shift(U[mu], nu, +1)                      # U_mu(x+nu)
            up = G.mul(Unu_xmu, G.dagger(G.mul(U[nu], Umu_xnu)))     # U_nu(x+mu) (U_nu(x) U_mu(x+nu))^+
            Unu_xmnu = self.shift(U[nu], nu, -1)                     # U_nu(x-nu)
            Umu_xmnu = self.shift(U[mu], nu, -1)                     # U_mu(x-nu)
            Unu_xmu_mnu = self.shift(Unu_xmu, nu, -1)                # U_nu(x+mu-nu)
            down = G.mul(G.dagger(G.mul(Umu_xmnu, Unu_xmu_mnu)), Unu_xmnu)  # (U_mu(x-nu) U_nu(x-nu+mu))^+ U_nu(x-nu)
            A = up + down if A is None else A + up + down
        return A

    # ------------------------------------------------------------------ observables
    def plaquette(self) -> float:
        """Average of ``P_p = (1/N) Re Tr U_p`` over all plaquettes."""
        tot = 0.0
        for mu in range(self.d):
            for nu in range(mu + 1, self.d):
                tot += float(np.mean(self.G.retrace(self.plaquette_field(mu, nu))))
        return tot / (self.d * (self.d - 1) / 2)

    def action(self) -> float:
        """``S = beta * sum_p (1 - P_p)``."""
        return self.beta * self.nplaq * (1.0 - self.plaquette())

    # ------------------------------------------------------------------ updates
    def _update_links_su2(self, mu: int, mask: np.ndarray, kind: str):
        A = self.staple_sum(mu)[mask]                 # (n, 4), A = xi * Ahat
        xi = su2.norm(A)
        Ahat = A / np.maximum(xi, 1e-300)[:, None]
        if kind == "heatbath":
            # weight exp((beta/2) Tr(U A)) = exp(beta xi (U Ahat)_0): draw a' ~ exp(k a'_0), U = a' Ahat^+
            k = self.beta * xi
            a = sample_su2(k, self.rng)
            Unew = su2.mul_dag(a, Ahat)
        elif kind == "overrelax":
            # U' = Ahat^+ U^+ Ahat^+ preserves (U Ahat)_0 and is a measure-preserving involution
            Uold = self.U[mu][mask]
            Unew = su2.mul(su2.dagger(Ahat), su2.mul(su2.dagger(Uold), su2.dagger(Ahat)))
        else:
            raise ValueError(kind)
        self.U[mu][mask] = Unew

    def _update_links_su3(self, mu: int, mask: np.ndarray, kind: str, nhits: int = 1):
        A = self.staple_sum(mu)[mask]                 # (n, 3, 3)
        self.U[mu][mask] = cabibbo_marinari(self.U[mu][mask], A, self.beta, self.rng, kind, nhits)

    def sweep(self, kind: str = "heatbath", nhits: int = 1):
        """One full sweep over all links (checkerboard order)."""
        for mu in range(self.d):
            for mask in (self.even, self.odd):
                if self.group == "su2":
                    self._update_links_su2(mu, mask, kind)
                else:
                    self._update_links_su3(mu, mask, kind, nhits=nhits)
        self._nsweeps = getattr(self, "_nsweeps", 0) + 1
        if self.group == "su3" and self._nsweeps % 20 == 0:
            self.U = su3.reunitarize(self.U)
        elif self.group == "su2" and self._nsweeps % 50 == 0:
            self.U = su2.normalize(self.U)

    def update(self, n_heatbath: int = 1, n_overrelax: int = 0):
        """A compound update: ``n_heatbath`` heat-bath sweeps then ``n_overrelax`` overrelaxation sweeps."""
        for _ in range(n_heatbath):
            self.sweep("heatbath")
        for _ in range(n_overrelax):
            self.sweep("overrelax")

    def thermalize(self, nsweeps: int, n_overrelax: int = 0):
        for _ in range(nsweeps):
            self.update(1, n_overrelax)

    # ------------------------------------------------------------------ gauge transformations
    def gauge_transform(self, g: np.ndarray | None = None):
        """Apply ``U_mu(x) -> g(x) U_mu(x) g(x+mu)^+`` (random ``g`` if none is given).

        Gauge-invariant observables must be unchanged; used by the tests."""
        if g is None:
            g = self.G.random(self.shape, self.rng)
        for mu in range(self.d):
            self.U[mu] = self.G.mul(g, self.G.mul_dag(self.U[mu], self.shift(g, mu, +1)))
        return g

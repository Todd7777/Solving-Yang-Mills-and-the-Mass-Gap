"""ymlattice: a small, tested lattice gauge theory toolkit for pure SU(2) and SU(3)
Yang--Mills theory with the Wilson action.

The package exists to support the Yang--Mills mass gap project: it gives us an
independent numerical probe of every conjecture we formulate and a sanity check
on every exact formula we derive (character expansions, strong- and
weak-coupling series, exact two-dimensional solutions).

Conventions (match theory/lattice-foundations.tex):

* Lattice: periodic hypercubic lattice with ``d`` directions and shape ``(L_0, ..., L_{d-1})``;
  direction 0 is Euclidean time.
* Link variables ``U_mu(x)`` in SU(N).  SU(2) elements are stored as unit quaternions
  ``(a0, a1, a2, a3)`` representing ``a0*1 + i*(a1*s1 + a2*s2 + a3*s3)``; SU(3) elements are
  stored as complex ``3x3`` arrays.
* Plaquette ``U_p = U_mu(x) U_nu(x+mu) U_mu(x+nu)^dagger U_nu(x)^dagger`` and
  ``P_p = (1/N) Re Tr U_p``.
* Wilson action ``S = beta * sum_p (1 - P_p)``, Boltzmann weight ``exp(-S)``.
  With this convention ``beta = 2N/g^2``.
"""

from .lattice import Lattice  # noqa: F401
from . import su2, su3, exact, observables, stats  # noqa: F401

__all__ = ["Lattice", "su2", "su3", "exact", "observables", "stats"]
__version__ = "0.1.0"

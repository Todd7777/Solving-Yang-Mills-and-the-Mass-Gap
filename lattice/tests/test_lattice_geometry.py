"""The simultaneous-update colour classes must be proper: no two mu-links updated together may
share a plaquette, i.e. a class must never contain both x and x +- e_nu for nu != mu."""
import numpy as np
import pytest

from ymlattice import Lattice


@pytest.mark.parametrize("shape", [(2, 2), (3, 3), (12, 3), (5, 3), (7, 4), (4, 4, 4, 4), (6, 5, 4, 3), (3, 3, 3, 3), (2, 3, 5)])
def test_colour_classes_are_proper_partitions(shape):
    lat = Lattice(shape, 1.0, "su2", seed=0)
    for mu in range(lat.d):
        classes = lat.colour_classes[mu]
        cover = np.zeros(shape, dtype=int)
        for mask in classes:
            cover += mask
            for nu in range(lat.d):
                if nu == mu:
                    continue
                for sign in (+1, -1):
                    assert not np.any(mask & np.roll(mask, sign, axis=nu)), (shape, mu, nu)
        assert np.all(cover == 1), "colour classes must partition the sites"


def test_even_lattices_keep_the_checkerboard():
    lat = Lattice((4, 6, 8, 2), 1.0, "su2", seed=0)
    for mu in range(4):
        assert len(lat.colour_classes[mu]) == 2
        assert np.array_equal(lat.colour_classes[mu][0], lat.even)


def test_rejects_length_one():
    with pytest.raises(ValueError):
        Lattice((4, 1), 1.0, "su2")

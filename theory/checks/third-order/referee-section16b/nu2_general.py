#!/usr/bin/env python3
"""nu_2 (Lemma mg:lem:threelow) for d = 3..7 by direct count, and a closed form fitted and checked.
nu_2 = number of unordered pairs of vertex-disjoint plaquettes != q sharing a vertex with q."""
from fractions import Fraction
from lattice16b import Lattice
from families16b import near, PV

vals = {}
for D in range(2, 7):
    lat = Lattice(D, 6)
    q = lat.plaq(lat.origin(), 0, 1)
    Nl = sorted(near(lat, q))
    nu2 = sum(1 for a in range(len(Nl)) for b in range(a + 1, len(Nl)) if not (PV(lat, Nl[a]) & PV(lat, Nl[b])))
    d = D + 1
    Pv = 2 * (d - 1) * (d - 2); KS = 4 * (2 * d - 5); Jv = 4 * (Pv - 1) - 2 * KS
    vals[d] = nu2
    print(f"d={d}: plaquettes sharing a vertex with q: {len(Nl)} (= 4(P_v-1)-K_S = {4*(Pv-1)-KS}), nu_2 = {nu2}, "
          f"N_3 = {3*nu2 + 12*(d-3)}, N_2 = {KS*KS//2 + KS*Jv + 8*(d-3)}, N'_2 = {KS*KS//2 + KS*Jv}")
# fit a polynomial of degree 4 in d through d=3..7 (Lagrange) and print its coefficients
ds = sorted(vals)
import itertools
def lagrange_coeffs(xs, ys):
    n = len(xs)
    coeffs = [Fraction(0)] * n
    for i in range(n):
        # basis polynomial
        num = [Fraction(1)]
        den = Fraction(1)
        for j in range(n):
            if j == i: continue
            num = [Fraction(0)] + num  # multiply by x
            for k in range(len(num) - 1):
                num[k] -= xs[j] * num[k + 1]
            den *= (xs[i] - xs[j])
        for k in range(n):
            coeffs[k] += ys[i] * num[k] / den
    return coeffs
c = lagrange_coeffs(ds, [vals[d] for d in ds])
print("interpolating polynomial nu_2(d) = " + " + ".join(f"({c[k]}) d^{k}" for k in range(len(c)) if c[k] != 0))

# closed form 8 (d-2)^2 (3d^2 - 14d + 17), checked at d = 3..7 and at d = 8 (not used in the fit)
for d in range(3, 8):
    assert vals[d] == 8 * (d - 2) ** 2 * (3 * d * d - 14 * d + 17)
lat = Lattice(7, 5)
q = lat.plaq(lat.origin(), 0, 1)
Nl = sorted(near(lat, q))
nu2_8 = sum(1 for a in range(len(Nl)) for b in range(a + 1, len(Nl)) if not (PV(lat, Nl[a]) & PV(lat, Nl[b])))
print(f"d=8 (L=5): nu_2 = {nu2_8}; closed form 8(d-2)^2(3d^2-14d+17) = {8 * 36 * (3 * 64 - 14 * 8 + 17)}")

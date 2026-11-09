# 0++ glueball, variational analysis: SU2 beta=2.3 lattice 8^3 x 16

- runs pooled: 2 (4000 + 4000 = 8000 measurements)
- plaquette = 0.60234 +- 0.00006; tau_int(plaquette) = 1.8 measurements
- Polyakov loop Re<L> = 0.0001 +- 0.0002 (zero in the confined phase)
- two-loop asymptotic scaling a*Lambda_L = 5.9144e-03

## Creutz ratios (4000 configurations with Wilson loops)

| chi(R,R) | value |
|---|---|
| chi(2,2) | 0.3178 +- 0.0003 |
| chi(3,3) | 0.2173 +- 0.0008 |
| chi(4,4) | 0.1775 +- 0.0067 |

## Effective masses

GEVP with tau0 = 0 over smearing levels [0, 2, 4, 8]; m_eff from the cosh ratio; block jackknife (20 blocks).

| tau | lvl 0 | lvl 2 | lvl 4 | lvl 8 | GEVP |
|---|---|---|---|---|---|
| 0 | 2.097 +- 0.029 | 1.512 +- 0.014 | 1.440 +- 0.014 | 1.475 +- 0.015 | 1.351 +- 0.013 |
| 1 | 1.538 +- 0.144 | 1.410 +- 0.068 | 1.342 +- 0.056 | 1.241 +- 0.046 | 1.278 +- 0.041 |
| 2 | 1.080 +- 0.290 | 1.241 +- 0.157 | 1.267 +- 0.127 | 1.159 +- 0.113 | 1.148 +- 0.112 |
| 3 | 1.585 +- 1.922 | 1.248 +- 0.708 | 1.572 +- 0.895 | 1.540 +- 0.640 | 1.261 +- 0.375 |
| 4 | - | - | - | 0.389 +- 0.870 | - |
| 5 | - | - | - | 0.718 +- 1.382 | 0.406 +- 0.542 |
| 6 | - | - | - | - | - |
| 7 | - | - | - | - | - |

GEVP effective mass at tau = 1: a*m = 1.278 +- 0.041  (an estimate of an upper bound on the gap); in units of the asymptotic-scaling scale: m/Lambda_L = 216 +- 7.
Ratio to the largest Creutz ratio: a*m / sqrt(chi(4,4)) = 3.03 (chi(4,4) over-estimates sigma at small R, so this under-estimates m/sqrt(sigma)).

Runs: two independent chains (seeds 101, 102), 4000 measurements each; commands in `results/glueball-su2.md`.

# 0++ glueball, variational analysis: SU2 beta=2.4 lattice 12^3 x 24

- runs pooled: 1 (2000 = 2000 measurements)
- plaquette = 0.62990 +- 0.00004; tau_int(plaquette) = 1.6 measurements
- Polyakov loop Re<L> = 0.0001 +- 0.0002 (zero in the confined phase)
- two-loop asymptotic scaling a*Lambda_L = 4.6004e-03

## Creutz ratios (500 configurations with Wilson loops)

| chi(R,R) | value |
|---|---|
| chi(2,2) | 0.2547 +- 0.0003 |
| chi(3,3) | 0.1501 +- 0.0006 |
| chi(4,4) | 0.1141 +- 0.0041 |
| chi(5,5) | 0.1101 +- 0.0215 |
| chi(6,6) | -0.0034 +- 0.0876 |

## Effective masses

GEVP with tau0 = 0 over smearing levels [0, 4, 8, 12]; m_eff from the cosh ratio; block jackknife (20 blocks).

| tau | lvl 0 | lvl 4 | lvl 8 | lvl 12 | GEVP |
|---|---|---|---|---|---|
| 0 | 2.293 +- 0.067 | 1.323 +- 0.025 | 1.192 +- 0.021 | 1.162 +- 0.020 | 1.114 +- 0.021 |
| 1 | 1.302 +- 0.173 | 1.036 +- 0.048 | 1.002 +- 0.040 | 0.986 +- 0.042 | 0.974 +- 0.039 |
| 2 | 1.303 +- 0.619 | 1.096 +- 0.199 | 0.973 +- 0.143 | 0.937 +- 0.123 | 0.877 +- 0.095 |
| 3 | - | 0.527 +- 0.187 | 0.892 +- 0.234 | 1.131 +- 0.304 | 0.753 +- 0.170 |
| 4 | - | 1.212 +- 0.770 | 1.010 +- 0.627 | 0.632 +- 0.485 | 0.397 +- 0.296 |
| 5 | 0.699 +- 1.171 | - | - | - | - |
| 6 | - | - | - | 0.263 +- 0.363 | 0.284 +- 0.385 |
| 7 | - | - | - | - | - |
| 8 | - | - | - | - | - |
| 9 | - | - | - | - | - |
| 10 | - | - | - | - | - |
| 11 | - | - | - | - | - |

GEVP effective mass at tau = 1: a*m = 0.974 +- 0.039  (an estimate of an upper bound on the gap); in units of the asymptotic-scaling scale: m/Lambda_L = 212 +- 8.

# 0++ glueball, variational analysis: SU2 beta=2.4 lattice 12^3 x 24

- runs pooled: 1 (350 = 350 measurements)
- plaquette = 0.62986 +- 0.00007; tau_int(plaquette) = 1.0 measurements
- Polyakov loop Re<L> = 0.0002 +- 0.0004 (zero in the confined phase)
- two-loop asymptotic scaling a*Lambda_L = 4.6004e-03

## Creutz ratios (88 configurations with Wilson loops)

| chi(R,R) | value |
|---|---|
| chi(2,2) | 0.2548 +- 0.0004 |
| chi(3,3) | 0.1492 +- 0.0015 |
| chi(4,4) | 0.1202 +- 0.0068 |
| chi(5,5) | 0.0996 +- 0.0496 |
| chi(6,6) | 0.2621 +- 0.2588 |

## Effective masses

GEVP with tau0 = 0 over smearing levels [0, 4, 8, 12]; m_eff from the cosh ratio; block jackknife (10 blocks).

| tau | lvl 0 | lvl 4 | lvl 8 | lvl 12 | GEVP |
|---|---|---|---|---|---|
| 0 | 2.235 +- 0.082 | 1.369 +- 0.046 | 1.221 +- 0.030 | 1.196 +- 0.039 | 1.153 +- 0.031 |
| 1 | 1.460 +- 0.471 | 1.155 +- 0.104 | 1.020 +- 0.104 | 0.947 +- 0.103 | 0.965 +- 0.084 |
| 2 | - | 1.385 +- 0.706 | 1.302 +- 0.229 | 1.217 +- 0.223 | 1.186 +- 0.232 |
| 3 | - | - | - | - | 0.559 +- 0.371 |
| 4 | - | - | - | - | - |
| 5 | - | - | - | - | 0.536 +- 0.657 |
| 6 | - | - | - | - | - |
| 7 | - | - | - | - | - |
| 8 | - | - | - | 0.640 +- 1.266 | 0.477 +- 0.514 |
| 9 | - | - | - | - | - |
| 10 | - | - | - | - | 1.861 +- 1.460 |
| 11 | - | - | - | - | - |

GEVP effective mass at tau = 1: a*m = 0.965 +- 0.084  (an estimate of an upper bound on the gap); in units of the asymptotic-scaling scale: m/Lambda_L = 210 +- 18.
Ratio to the largest Creutz ratio: a*m / sqrt(chi(6,6)) = 1.88 (chi(6,6) over-estimates sigma at small R, so this under-estimates m/sqrt(sigma)).

**Partial run:** a snapshot of the first 350 of 2000 planned measurements of one chain (seed 201), analysed with 10 jackknife
blocks; preliminary. Command in `results/glueball-su2.md`; the chain resumes from its checkpoint when rerun.

# 0++ glueball, variational analysis: SU2 beta=2.4 lattice 12^3 x 24

- runs pooled: 1 (700 = 700 measurements)
- plaquette = 0.62983 +- 0.00007; tau_int(plaquette) = 1.7 measurements
- Polyakov loop Re<L> = 0.0002 +- 0.0003 (zero in the confined phase)
- two-loop asymptotic scaling a*Lambda_L = 4.6004e-03

## Creutz ratios (175 configurations with Wilson loops)

| chi(R,R) | value |
|---|---|
| chi(2,2) | 0.2546 +- 0.0004 |
| chi(3,3) | 0.1497 +- 0.0010 |
| chi(4,4) | 0.1170 +- 0.0055 |
| chi(5,5) | 0.1110 +- 0.0349 |
| chi(6,6) | 0.0821 +- 0.1549 |

## Effective masses

GEVP with tau0 = 0 over smearing levels [0, 4, 8, 12]; m_eff from the cosh ratio; block jackknife (20 blocks).

| tau | lvl 0 | lvl 4 | lvl 8 | lvl 12 | GEVP |
|---|---|---|---|---|---|
| 0 | 2.167 +- 0.064 | 1.348 +- 0.034 | 1.215 +- 0.024 | 1.198 +- 0.027 | 1.146 +- 0.024 |
| 1 | 1.353 +- 0.302 | 1.089 +- 0.085 | 1.048 +- 0.071 | 1.017 +- 0.070 | 1.015 +- 0.050 |
| 2 | 1.164 +- 1.090 | 0.860 +- 0.212 | 0.984 +- 0.134 | 1.036 +- 0.145 | 0.941 +- 0.151 |
| 3 | 0.876 +- 2.759 | 0.652 +- 0.360 | 0.992 +- 0.502 | 1.297 +- 0.753 | 0.739 +- 0.360 |
| 4 | - | - | - | - | - |
| 5 | - | - | - | - | 0.362 +- 0.329 |
| 6 | - | - | - | 0.478 +- 0.833 | - |
| 7 | - | - | - | 0.571 +- 2.421 | - |
| 8 | - | - | - | - | - |
| 9 | - | - | - | - | - |
| 10 | - | - | - | - | 1.388 +- 1.321 |
| 11 | - | - | - | - | - |

GEVP effective mass at tau = 1: a*m = 1.015 +- 0.050  (an estimate of an upper bound on the gap); in units of the asymptotic-scaling scale: m/Lambda_L = 221 +- 11.
Ratio to the largest Creutz ratio: a*m / sqrt(chi(6,6)) = 3.54 (chi(6,6) over-estimates sigma at small R, so this under-estimates m/sqrt(sigma)).

**Partial run:** a snapshot of the first 700 of 2000 planned measurements of one chain (seed 201), analysed with 20 jackknife
blocks; preliminary. Command in `results/glueball-su2.md`; the chain resumes from its checkpoint when rerun.

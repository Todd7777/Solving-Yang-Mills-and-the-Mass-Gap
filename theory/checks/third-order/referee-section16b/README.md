# Referee checks for Part IV, Section 16: Lemmas 16.3 and 16.4 (`mg:lem:threelow`, `mg:lem:jsums5`)

Independent checks by the referee of the text of `theory/mass-gap-strong-coupling.tex` from "With $\theta_7$, the largest term of $\delta'$ at
practical $\beta$" through the end of the proof of Lemma `mg:lem:jsums5` (facts (F1), (F2); Lemma `mg:lem:threelow`; Lemma `mg:lem:jsums5`).
The code was written from the definitions of the paper (Definitions `mg:def:fronts`, `mg:def:junction`, Lemma `mg:lem:junction`) before the authors'
`three_front_count.py` was read; it uses none of the paper's lemmas. A compass, not part of the proofs.

Run from this directory with `python3` (standard library only, except `rmax_gain.py` and `exact_part_size.py`, which need `mpmath` and the authors'
modules through `toe_rmax.py`). Total CPU time of all runs: about 35 minutes. Outputs are in
`outputs/`. Everything is on the periodic spatial lattice $(\mathbb Z/L)^D$, $D=d-1$; $L=8$ is the paper's hypothesis "all spatial sides $\ge8$", and
larger $L$ is equivalent to $\mathbb Z^D$ for these configurations (they span at most $6$ in each coordinate).

## The conditions used (from the definitions only)

A candidate junction at slice $s$ is a pair $(P,\Phi)$, $P$ a set of plaquettes of the slice, $\Phi$ a set of front-pairs $(\varphi,\mu)$,
$\mu\in\{s-1,s\}$, with: every $\varphi$ a front (non-empty, connected, all degrees $\ge2$); distinct fronts of one layer vertex-disjoint; $n_J(e)\ge2$ on
$E(J)$; every front-pair sharing a vertex with $E(P)$ (no support); and $P\cup\Phi$ connected under "common vertex" (equivalent to the connectivity of
$H_s$). These conditions are necessary for a junction without support, so the counts are upper bounds for the counts of the paper; they agree.

## Scripts

| Script | What it does | Result | Runtime |
|---|---|---|---|
| `lattice16b.py`, `fronts16b.py` | Lattice; fronts through a vertex by a deficiency search from the definition. `python3 fronts16b.py D L KMAX` | $\mathbb Z^2$: $N^v_{4..9}=4,0,12,12,70,160$; $\mathbb Z^3$: $N^v_{4..8}=12,0,132,108,2034$ (as in Remark 16.5(ii)) | 1 s |
| `check_F12.py D L` | Every six-edge front through a vertex is a cycle and exactly one of planar domino / bent domino / skew hexagon; (F1); the three bullets of (F2); every seven-edge front is $\partial r_0\cup\partial c_0$ for exactly one pair sharing a link, and only the two faces perpendicular to $e_0$ of a bent pair share two links with it | all true for $D=2,3,4$, $L=8$ and $12$ ($\mathbb Z^3$: 36 planar, 72 bent, 24 skew; 36 coplanar and 72 bent seven-edge fronts) | seconds |
| `enum_junctions.py D L n sizes lam` | Exhaustive enumeration of the candidate junctions containing $p=(\partial q,\lambda)$ with $\lvert P\rvert=n$ and further front-pairs of the given sizes (deficiency search: branch at the smallest link covered once over all items containing it; at a closed state over all items touching it; every result re-checked by `is_junction`) | validation on known counts: $Z_{n,4}=1,8,0,0$ ($d=3$), $1,32,0,12$ ($d=4$); one-plaquette exits with $6,7,8$ edges $4,4,4$ ($d=3$), $12,12,20$ ($d=4$), none with $5$ or $9$ (`outputs/validation_two_frontpairs.txt`) | < 1 s each |
| `families16b.py` | The families of Lemma 16.3 built literally from their descriptions: caps (p alone / p with $q_2$), cubes (two kinds), (a), (b) with $r=q$ and $r=q'$, (c), and (iii) | used by the next script | |
| `analyze_threelow.py D L [i,j]` | Search against families, as sets, both layers $\lambda$; table of the pairs $(J,p')$ by $(n,E)$ and by the layer of $p'$; negative checks ($\lvert P\rvert=2$ with $4,4,4$ or $4,4,5$; $\lvert P\rvert=1$ with other sizes $\le7$) | identical sets: $d=3$: $48$, $24$, $24$; $d=4$: $876$, $320$, $312$; $d=5$: $4776$, $1256$, $1240$ (= the general formulas with $\nu_2=1584$); same on $L=8,12,16$ and for the three orientations of $q$ in $d=4$; all negative checks $0$ | 1 s ($d=3$), 5 s ($d=4$), 40 s ($d=5$) |
| `brute_pfirst.py D L n sizes lam` | A structurally different enumeration: all $P$ with a plaquette at a vertex of $q$, then fronts covering the links that need one, then `is_junction` | identical sets for all cases in $d=3$ (both layers) and for (ii), (iii), $\lvert P\rvert=2$ with $4,4,4$, and $\lvert P\rvert=1$ in $d=4$ | 2 s ($d=3$); 2–3 min per case ($d=4$) |
| `brute_pfirst_i.py D L lam` | The same for case (i) in $d=4$, with the extra (justified) pruning $\lvert T\rvert\le8$ | $1177513$ plaquette sets kept; $876$ junctions, identical set | 5 min |
| `variant_orders.py D L` | The deficiency search with the largest, or a random, deficient link | identical sets | 10 s |
| `nu2_general.py` | $\nu_2$ for $d=3,\dots,8$; the triple identity $\#\{(q,b,c)\}=2\nu_2$ is checked in `analyze_threelow.py` | $\nu_2=16,288,1584,5248,13200,27936$ $=8(d-2)^2(3d^2-14d+17)$ (fitted on $d=3..7$, confirmed at $d=8$) | 10 s |
| `rmax_gain.py d beta [maxonly]` (with `toe_rmax.py`, a copy of the authors' `../third_order_effective.py` with one added option `RMAX`) | $\delta'$ of Remark 16.5(iii) ($SU(3)$, settings of `../practical_beta.py`) with the factor $R$ of Lemma 16.4 (and the two-term factor of $M_4$) replaced by the maximum of its terms, which is also valid | $\delta'/\beta^4$: $d=4$, $\beta=1.7\cdot10^{-6}$: $1170$ (sum, as in the paper) $\to439$ (max); $d=3$, $\beta=1.8\cdot10^{-5}$: $107.5\to40.2$. With the max, $\beta^3/432$ is resolved up to about $2.4\cdot10^{-6}$ ($d=4$; at $2.5\cdot10^{-6}$: $949$ against $926$) and beyond $2.6\cdot10^{-5}$ ($d=3$; $83.9$ against $89.0$), instead of $1.7\cdot10^{-6}$ and $1.8\cdot10^{-5}$ (`outputs/rmax_gain_d*.txt`, `outputs/rmax_threshold_d*.txt`) | 1–5 min per $\beta$ |
| `exact_part_size.py` | The exact part $e^{c_J}(2N_3b^3a^2+N_2b^2a^3+N'_2b^2a^{7/2})$ of $M_5$ at the same $\beta$ ($c_J\le1.1$) | $\le2.2\cdot10^{-4}\beta^4$ ($d=4$) and $1.3\cdot10^{-4}\beta^4$ ($d=3$), i.e. about $10^{-7}$ of $\delta'$: the gain of Lemma 16.4 comes from the remainder starting at the corners $(4,2),(3,3),(2,4)$, not from the values of $N_3,N_2,N'_2$ | 1 s |

## Main findings

* Lemma 16.3 is confirmed family by family: the candidate junctions with exactly three front-pairs containing a given four-edge front-pair are exactly the
  caps and cube configurations of (i), the configurations (a), (b), (c) of (ii) and (a), (b) of (iii), in $d=3,4,5$, on the torus of side $8$ as well as
  on larger ones, in both layers. No family is missing and every listed configuration is a junction, so the bounds $N_3,N_2,N'_2$ are equalities.
* (F1), (F2) and the structure of six- and seven-edge fronts hold as stated (in $\mathbb Z^2,\mathbb Z^3,\mathbb Z^4$).
* Lemma 16.4: the pairs $(J,p')$ of these junctions have $(n,E)=(3,2)$ ($2N_3$ pairs), $(2,3)$ and $(2,4)$ ($N_2$ each), $(2,\frac72)$ and $(2,5)$ ($N'_2$ each),
  as the proof says. All pairs with $(n,E)=(2,3)$ or $(2,\frac72)$, and $2\nu_2+8(d-3)$ of the $2N_3$ pairs with $(3,2)$, have $p'$ in the layer of $p$.
* Outside the hypothesis: for sides $5$, $6$, $7$ the counts in $d=4$ are unchanged; with sides $4$ there are more (straight four-edge fronts)
  (`outputs/small_sides_diagnostic.txt`).

## Commands used

```
python3 fronts16b.py 2 16 9 ; python3 fronts16b.py 3 12 8
for a in "2 8" "2 12" "3 8" "3 12" "4 8"; do python3 check_F12.py $a; done
python3 analyze_threelow.py 2 8 ; python3 analyze_threelow.py 3 8 ; python3 analyze_threelow.py 4 8
python3 analyze_threelow.py 2 16 ; python3 analyze_threelow.py 3 12 ; python3 analyze_threelow.py 3 8 0,2 ; python3 analyze_threelow.py 3 8 1,2
for L in 7 6 5 4; do python3 analyze_threelow.py 3 $L; done
python3 brute_pfirst.py 2 8 <n> <sizes> <lam> ; python3 brute_pfirst.py 3 8 2 4,6 0 ; python3 brute_pfirst.py 3 8 2 4,7 0 ; ...
python3 brute_pfirst_i.py 3 8 0 ; python3 variant_orders.py 3 8 ; python3 nu2_general.py
python3 rmax_gain.py 4 1.7e-6 ; python3 rmax_gain.py 3 1.8e-5 ; python3 rmax_gain.py 4 2.5e-6 maxonly ; python3 rmax_gain.py 3 2.6e-5 maxonly
python3 exact_part_size.py
```

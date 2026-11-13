# Referee checks for Part IV, Section 16 (fifth pass): Lemma `mg:lem:uppertails`, Remark `mg:rem:upperexact`, Remark `mg:rem:nextorder`(viii)

Independent checks, made on 13 November, of these parts of `theory/mass-gap-strong-coupling.tex`:

* Lemma `mg:lem:uppertails` (the tails of the upper bound by exact counts; Theorem-level), with the background of Lemmas `mg:lem:five`, `mg:lem:twopoint`,
  `mg:lem:closed`, `mg:lem:char`(ii), `mg:lem:capelement`(iv) and Theorem `mg:thm:upper2`.
* Remark `mg:rem:upperexact` (i)–(iv) (Numerical): the counts of the admissible sets, the bound $E(\beta)$, the two-sided statement for $SU(3)$ in $d=3$, and $d=4$.
* Remark `mg:rem:nextorder`(viii), the new sentence of the abstract and the two new rows of the summary table.

These checks are a compass, not part of the proofs. The enumerations were written from the definition of "admissible" in Lemma `mg:lem:uppertails`, with two
different algorithms, and use none of the authors' code (`rel_closed_count.py` was read only after our counts were complete). The evaluation re-implements
$\mathcal T'_S$, Theorem `mg:thm:upper2` and $E(\beta)$ from the paper; the $SU(3)$ constants $u,\Upsilon_\pm,\kappa$ are computed by our own Weyl integration and
agree with the authors' `second_order_constants.su3` to 12 digits.

## Scripts

| Script | What it does | Result | Runtime |
|---|---|---|---|
| `admissible_partition.py d A\|Q nmax [--torus L] [--list]` | Method A. Connected admissible sets by a partition search: root = smallest plaquette of $\mathcal P_S\setminus\mathcal C$ in the set; branch on the first candidate (include it, exclude the earlier ones) through the smallest deficient link, or, at a closed state, through any plaquette sharing a link with the set; each set reached exactly once (asserted). Pruning: $\lceil D/4\rceil$ more plaquettes for $D$ deficient links. Then unions of link-disjoint connected ones. `--torus L`: the periodic lattice with all sides $L$ | see tables | 0.1–7 s |
| `admissible_redelmeier.c d A\|Q nmax [2]` | Method B. Redelmeier enumeration of all sets $Y$ (outside $\mathcal C$) such that $Y\cup\{o\}$ is connected, $o$ a virtual root joined to $\mathcal P_S\setminus\mathcal C$, i.e. conditions (1) and (3); condition (2) tested at every node; multi-component sets generated directly. Pruning $\lceil D/4\rceil$, and with `2` also a greedy set of deficient links pairwise in no common plaquette | see tables | 7 s – 18 min |
| `classify_sets.py d nmax` | shapes of the sets found by Method A | see below | 1–10 s |
| `cap_bruteforce.py` | the five-face cap: sum over all $\lvert G\rvert^8$ configurations of the eight inner links (no gauge fixing, random bottom links) for $\mathbb Z_3,\mathbb Z_4,S_3$; vertical links gauge fixed for the Frobenius group $F_{21}$ (complex 3-dim irreps with trivial determinant, like $SU(3)$); two values of $\beta$ each | $g=f^{*5}=\sum_{r\ne0}d_ru_r^5\chi_r$ to $10^{-16}$, $\max\lvert g\rvert=g(1)=\sigma_5=\sum d_r^2u_r^5$ | 6 s |
| `upper_eval.py` | $\mathcal T'_S$ of the Lemma; the exact error $D'_2+\log\Upsilon$ of Theorem `mg:thm:upper2`; $E(\beta)$ of Remark (ii); monotonicity; thresholds in $d=4$ with counts to 8, 9, 10 | see below | 40 s |
| `ablation.py` | $E(\beta_0(3))$ without each ingredient of the Lemma, and with the unit cubes also evaluated exactly | see below | 5 s |
| `crude_check.py` | the crude-tail error of Remark `mg:rem:nextorder`(viii) | thresholds $2.967\cdot10^{-8}$ ($d=3$), $6.706\cdot10^{-9}$ ($d=4$); coefficients $2.630\cdot10^{12}\beta^5$, $5.147\cdot10^{13}\beta^5$ | 5 s |

## Counts of admissible sets (non-cap sets by size; caps separately)

| lattice | $(S,\mathcal C)$ | $N^{\rm cap}$ | $n=6$ | 7 | 8 | 9 | 10 | 11 | 12 | methods |
|---|---|---|---|---|---|---|---|---|---|---|
| $\mathbb Z^3$ | $(S_A,\{p_0\})$ | 2 | 8 | 0 | 0 | 10 | 47 | 60 | 96 (16 two-comp.) | A to 12, B to 12 |
| $\mathbb Z^3$ | $(S_Q,\mathcal Q)$ | 6 | 12 | 0 | 0 | 30 | 93 (3 two-comp.) | 144 (48) | 162 (42) | A to 12, B to 12 |
| $\mathbb Z^4$ | $(S_A,\{p_0\})$ | 4 | 32 | 0 | 0 | 60 | 506 | | | A to 10, B to 10 |
| $\mathbb Z^4$ | $(S_Q,\mathcal Q)$ | 18 | 60 | 0 | 24 | 262 | 1143 (27 two-comp.) | | | A to 10, B to 10 |

No admissible set has fewer than five plaquettes (Lemma `mg:lem:five`), and every five-plaquette set is a cap. The paper's counts (Remark (i): $\mathbb Z^3$ to 10;
Remark (iv): $\mathbb Z^4$ to 8) are reproduced by both methods; the other entries are new, and every entry of the table is given by both methods. Shapes (`outputs/classify.txt`): the six-plaquette sets are unit cubes sharing a link with $S$ and containing no plaquette of
$\mathcal C$; the nine-plaquette sets in $\mathbb Z^3$ are boxes $2\times1\times1$ minus their unique face in $\mathcal C$; the two-component ten-plaquette sets are pairs
of caps on opposite faces of $Q$; in $\mathbb Z^4$ the 24 eight-plaquette sets for $(S_Q,\mathcal Q)$ are $\partial\Gamma_1\triangle\partial\Gamma_2$ minus $F_1,F_2$, for cubes
$\Gamma_i\supset F_i$ on two adjacent faces of $Q$ that share a face.

Tori (`outputs/torus_d3_10.txt`, Method A with `--torus L`, up to 10 plaquettes): the counts equal those of $\mathbb Z^3$ for every $L=4,\dots,13$; for $L=3$ they differ
($S_A$: one set of 8, the wrapped $3\times3$ plane through $p_0$ minus $p_0$; 14 of 9; 43 of 10). So the claim "sides $\ge12$" holds with room.

## The cap

`outputs/cap_bruteforce.txt`: for $\mathbb Z_3$, $\mathbb Z_4$, $S_3$ (all $\lvert G\rvert^8$ inner configurations) and $F_{21}$, at $\beta=0.3$ and $1.5$, the integral over the
eight inner links equals the five-fold convolution $f^{*5}=\sum_{r\neq0}d_ru_r^5\chi_r$ (difference $\le2\cdot10^{-16}$, independent of the non-gauge-fixed bottom
links), and $\max_V\lvert g(V)\rvert=g(1)=\sum_{r\ne0}d_r^2u_r^5=\sigma_5$.

## Evaluation, $SU(3)$ (`outputs/upper_eval.txt`)

* $d=3$, $n_0=10$, at $\beta_0(3)=4.0262\cdot10^{-4}$: $\mathcal T'_{S_Q}=1.232\cdot10^{-3}\beta^4$ (tail $1.106\cdot10^{-3}$, sets $1.25\cdot10^{-4}$, caps $2.3\cdot10^{-8}$),
  $\mathcal T'_{S_A}=5.05\cdot10^{-4}\beta^4$; $q=0.02628$, $\sigma_5/u^5=18.000$, $\epsilon/u=36.01$; exact error $D'_2+\log\Upsilon=0.40386\beta^4$, $E=0.40392\beta^4$
  (paper: $0.404$). $B$ is $99.7\%$ of $E$, and its term $(1+a'^2)e^{16\beta}\frac{N^4}{\tau_0^2}\mathcal T'_{S_Q}$ ($N^4/\tau_0^2=324$) is $99.7\%$ of $B$.
* Monotonicity on 91 points of $[4\cdot10^{-13},\beta_0(3)]$: $E/\beta^5$, $u/\beta$ and $\bar\sigma_5/\beta^5$ increase; $\Upsilon\ge1$ and $\kappa>18u^2$; exact error $\le E$.
* $\beta^3/432=5.7494\beta^4$ at $\beta_0(3)$ (the paper's "$\ge5.75\beta^4$" in Remark (iii) should be $\ge5.74\beta^4$).
* $d=3$, $n_0=12$ (our counts): $E(\beta_0)=0.0415\beta_0^4$ (then the cubes dominate).
* $d=4$, $n_0=8$: error $<\beta^3/432$ for $\beta\le8.261\cdot10^{-5}=\beta_0(4)/1.86$; $1.97\cdot10^{-3}\beta^4$ at $1.2\cdot10^{-5}$ (paper: $8.2\cdot10^{-5}$, $2\cdot10^{-3}$).
* $d=4$, $n_0=9$: $E(\beta_0(4))=10.6\beta^4<\beta^3/432=15.05\beta^4$, $E/\beta^5$ increasing on 91 points: **the upper side resolves the term on the whole range in $d=4$**;
  $n_0=10$: $E(\beta_0(4))=0.207\beta^4$.
* Crude tails (Remark `mg:rem:nextorder`(viii)): at $\beta_0(3)$, $\mathcal T_Q=1.381\cdot10^{-4}=5.25\cdot10^9\beta^4$ (text: $5.2\cdot10^9$); thresholds $3.0\cdot10^{-8}$, $6.7\cdot10^{-9}$
  and coefficients $2.63\cdot10^{12}$, $5.15\cdot10^{13}$ (text: $5.2\cdot10^{13}$).

## Which ingredients matter (`outputs/ablation.txt`, $SU(3)$, $d=3$, $n_0=10$, at $\beta_0(3)$; $\beta^3/432=5.749\beta^4$)

| variant | $E(\beta_0)/\beta_0^4$ |
|---|---|
| as in the paper ($\theta_5$, caps by characters) | $0.404$ |
| caps counted with $w^5$ instead of $\sigma_5$ | $25.9$ (the term is **not** resolved) |
| $\theta$ instead of $\theta_5$ | $2.3\cdot10^5$ |
| the six-plaquette sets (unit cubes meeting $S$ in one link) weighted by $\sigma_6$, $n_0=10$ | $0.363$ |
| the same with $n_0=12$ (our counts) | $2.7\cdot10^{-4}$ |

A unit cube $\Gamma$ of an admissible set meets $S$ in exactly one link $\ell$ (two links of $S$ in $\Gamma$ span a face of $\mathcal C$, or $\Gamma=Q$); its other eleven links lie in
no other plaquette of $X$, and integrating them gives a function of $U_\ell$ invariant under $U_\ell\mapsto hU_\ell k^{-1}$, i.e. the constant $\sigma_6$. This is a suggestion, not used
in the paper.

## Commands

```
S=<scratch directory>
python3 admissible_partition.py 3 A 12; python3 admissible_partition.py 3 Q 12        # and d = 4 up to 10
for L in 3 4 5 6 7 8 9 10 11 12 13; do python3 admissible_partition.py 3 Q 10 --torus $L; done
gcc -O2 -o $S/admissible_redelmeier admissible_redelmeier.c
$S/admissible_redelmeier 3 A 10; $S/admissible_redelmeier 3 Q 10; $S/admissible_redelmeier 4 A 9; $S/admissible_redelmeier 4 Q 9
$S/admissible_redelmeier 4 A 10 2; $S/admissible_redelmeier 4 Q 10 2; $S/admissible_redelmeier 3 A 12 2; $S/admissible_redelmeier 3 Q 12 2
python3 classify_sets.py 3 10; python3 classify_sets.py 4 9
python3 cap_bruteforce.py; python3 upper_eval.py; python3 ablation.py; python3 crude_check.py
```

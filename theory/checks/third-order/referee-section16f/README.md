# Referee checks for Part IV, Section 16 (sixth pass): Remark `mg:rem:nextorder`(ix), the lower bound in $d=4$

Independent checks, made on 13 November, of these parts of `theory/mass-gap-strong-coupling.tex`:

* Remark `mg:rem:nextorder`(ix) ("The lower bound in $d=4$"): the capacity prune of the junction search, the twelve classes of order 7 and the classes of order 8
  with pairs of at least four plaquettes in $d=4$, $S''_8$ and its corners, $Z_{8,4}$, $Z_{9,4}$, $N^v_{11}$, the exits $Z_{3,k}$, the closed polymers of $10$ and $11$
  plaquettes of $\mathbb Z^4$, the evaluation and the scaling below the threshold;
* the row of the summary table and the $d=4$ part of the abstract sentence on the two-sided third-order term.

The text changed during the review: first $S''_8=S(7)\cup\{n+E\le\frac{17}2,\ n\ge4\}$ and the threshold $4.3\cdot10^{-5}$ (`orders/d4_s8.json`), then the same set with the
exponent $(3,5)$ added (four classes with three plaquettes) and the threshold $5\cdot10^{-5}$ (`orders/d4_s8b.json`). Both versions are checked below.

These checks are a compass, not part of the proofs. The counts below were made with code written for this pass from Definition `mg:def:junction` and
Lemma `mg:lem:junction` (fronts, junctions, exits, closed sets); it uses neither the authors' nor the earlier referees' code or lists of fronts. The authors' search
(`../order_counts.py`) was run only to test its prune against itself, and the authors' evaluator (`../practical_beta.py`, `../third_order_effective.py`) only
to reproduce the evaluation and to measure which ingredients carry it.

## The capacity prune (`../order_counts.py`, `count_junctions`, `PRUNE`)

*Argument.* Let $X=(P_X,\Phi_X)$ be a state of the search contained in a junction $J$ that the search is to find ($|P_J|\le$ `NMAX`, further fronts of exactly the
class sizes, no link of coverage one). Coverage only grows along the search, so every link of coverage one in $X$ has coverage at least two in $J$ and lies in an
element of $J\setminus X$; an element contains at most $4$ (plaquette) or $|G|$ (front) links. Hence
$$\#\{\text{links of coverage one in }X\}\ \le\ 4\,(\texttt{NMAX}-|P_X|)+\sum(\text{sizes of the fronts of }J\setminus X)=\texttt{capacity}(P_X,\texttt{rem}).$$
The search only adds elements, so every state on the path to $J$ is contained in $J$. The three uses are this inequality for:
* the state itself (start of `rec`);
* $X\cup\{H\}$ for a front $H$ seeded at a closed state: its links of coverage one are exactly $H\setminus U$, its capacity `capf - len(H)`;
* $X\cup\{H\}$ for a front candidate $H$ at a deficient link: its links of coverage one are $(D_1\setminus H)\cup(H\setminus U)$, of size `len(D1) - len(H & D1) + len(H - U)`.

In the last two, $H$ is dropped from every option (with or without plaquettes $T$ and a second front); every junction containing such an option contains $X\cup\{H\}$,
which violates the inequality. Excluded plaquettes, fronts of one layer, a front in both layers, $p$ itself: none enters the argument, which uses only that
coverage never decreases and that $J\setminus X$ consists of at most `NMAX`$-|P_X|$ plaquettes and the remaining fronts. Closed states are never cut, so `seen`
and `finish` are unaffected. **The prune is a valid necessary condition in all three places.**

*Tests.*
* `prune_ab.py`: the authors' `count_junctions` with `PRUNE = True` and `False`, nothing else changed.
  * $d=3$ (`outputs/prune_ab_d3_a.txt`, `prune_ab_d3_b.txt`): **34 classes, identical counts in all**: the 22 multisets of order 8, each to its $n_{\max}$, and 12 runs with
    larger $n_{\max}$ (e.g. $4,4$ to $n=8$; $4,6$, $4,7$ to $6$; $4,4,4$ to $5$; $4,8$, $4,9$ to $5$; $4,10$, $4,11$, $4,4,6$, $4,4,7$, $6,6$, $6,7$ to $4$). Every count that the
    fourth referee made (`referee-section16d`, a third method) is the same; $4,4$ with $n\le8$ gives only $48$ at $n=3$, as Remark (iv) proves. Speed-up $1.1$ to $18$ ($4,13$: $68$ s against $1197$ s).
  * $d=4$ (`outputs/prune_ab_d4.txt`): **14 runs, identical in all**: $4,4,4$ ($n\le3$ and $\le4$), $6,6$, $6,7$, $4,4,6$, $4,4,7$, $4,8$, $4,9$ ($n\le2$), $4,6$, $4,7$
    ($n\le3$ and $\le4$), $4,4$ ($n\le5$), $4,8$ ($n\le3$); $7$ to $21$ times faster with the prune ($4,8$, $n\le3$: $51$ s against $1069$ s). The runs without the prune of the
    classes with fronts of $9$ to $11$ edges and three or more plaquettes take hours; they were not made (the referee's own search covers them).
* `junc_partition.c` (below), a different search, with and without its own implementation of the inequality (`noprune`): identical counts in all $d=4$ classes of
  the table and in $Z_{n,4}$, $n\le9$, with $20$ to $120$ times more nodes without it (`outputs/junc_partition_d4_noprune.txt`, `junc_partition_Zn4_d4.txt`).
* Above all, the independent counts below equal the authors' pruned counts in every $d=4$ class used: a prune that lost junctions would show there as a deficit.

## Scripts

| Script | What it does | Result | Runtime |
|---|---|---|---|
| `fronts_enum.c D K [dump]` | Fronts (connected, all degrees $\ge2$) of $\mathbb Z^D$ up to translation, by Redelmeier's algorithm on edges rooted at the lexicographically smallest vertex (dead-vertex and degree prunes); $N^v_k=\sum_{\rm classes}\lvert V\rvert$; dumps the classes | $\mathbb Z^3$: $12,132,108,2034,4344,43650,$ **$147888$** ($k=4,6,\dots,11$), and $N^v_{12}=1121083$, $N^v_{13}=4773186$; $\mathbb Z^2$: $4,12,12,70,160,634,1728,6369,18704,65228,199224,671188$ ($k\le16$) | 0.9 s ($\mathbb Z^3$, 11); 17 s (13) |
| `junc_partition.c DUMP D NMAX SIZES [noprune]` | Candidate junctions $\ni p=(\partial q,s-1)$ with conditions (1)–(5) of its header (Definition `mg:def:junction`, Lemma `mg:lem:junction`; front–front adjacency included), by a *partition* search: at the smallest link of coverage one, branch $i$ = "$c_i$ is the first element through it that $J$ contains" ($c_1,\dots,c_{i-1}$ excluded, in a LIFO hash set); at a closed state, $J=X$ or the first element of $J$ sharing a vertex with $X$. Each junction is reached once (fingerprints: $0$ duplicates in every run); every candidate is re-checked from scratch | see the tables | 0.04–74 s per class |
| `zexit_pfirst.py D NMAX` | $Z_{n,k}$ plaquettes first: $P$ with $\{q\}\cup P$ vertex-connected (derived in the header from $n_J\ge2$ and connectivity), $\varphi'=\mathrm{Def}(P)\cup T$ ($T$ any set of doubly covered links of $\partial q\cup E(P)$), both layers, all conditions from scratch | $d=4$: **$Z_{3,k}=288,264,862,760,3490,6596,14838,18444,23672,15296,10868$** ($k=6,\dots,16$); $Z_{2,k}=20,12,182,296,606,456,458$; $Z_{1,k}=12,12,20$; none in the layer of $p$. $d=3$: the paper's $Z_{1,k},Z_{2,k},Z_{3,k}$ | 12 s |
| `closed_count.py d NMAX` | link-connected closed sets through a plaquette, partition search | $\mathbb Z^4$ to 11: $\{6{:}4,\ 10{:}60,\ 11{:}66\}$, none with $7,8,9$; $\mathbb Z^3$ to 12: $\{6{:}2,\ 10{:}10,\ 11{:}11,\ 12{:}24\}$ | 3 s |
| `assemble_check_d4.py` | the classes of $S(7)$, $S''_8$ and $S''_8\cup\{(3,5)\}$ enumerated directly from the definition of $E$; the terms rebuilt from the referee's counts; the corners as minimal exponents outside the set | **equal to `orders/d4_m7.json` (21 terms), `orders/d4_s8.json` (28) and `orders/d4_s8b.json` (32)**; corners $(2,6),(3,5),(5,4),(6,3),(7,2)$, resp. $(2,6),(3,\frac{11}2),(4,5),(5,4),(6,3),(7,2)$ | 1 s |
| `prune_ab.py MODE D K:NMAX ...` | the authors' search with the prune on and off | identical in every class | see outputs |
| `ablation_d4.py [beta] [--terms F]` | the authors' evaluator with one ingredient removed or inflated | see below | 4 min |
| `corner_factors_d4.py [beta] [--scanA ..]` | the factor $(b/\bar b)^n(a/\bar a)^E$ of each corner at the chosen $(\bar b,\bar a)$, and two hypothetical larger sets | see below | 2 min |
| `same_layer_exits.py D N` | the exits with exactly $N$ plaquettes and $p'$ in the layer of $p$ (the enumeration of `zexit_pfirst.py`, restricted to $P\supseteq$ a cover of $\partial q$) | $d=4$, $N=4$: $8$, all six-edge | 39 s |
| `bonus_eval_d4.py beta ..` | the authors' evaluator with `d4_s8b.json` and, in addition, the referee's $N^v_{12},N^v_{13},Z_{10..12,4}$ (next step only) | see below | 1 min |

## Junction counts in $d=4$ (`outputs/junc_partition_d4*.txt`), referee = authors in every entry

| further fronts $K$ | $n=1$ | $n=2$ | $n=3$ | $n=4$ | $n=5$ | $n=6$ |
|---|---|---|---|---|---|---|
| $4,4$ | 0 | 0 | 876 | 42 | 912 | 1866 |
| $4,6$ | 0 | 320 | 680 | 18276 | 11744 | |
| $4,7$ | 0 | 312 | 344 | 17028 | 3440 | |
| $4,8$ | 36 | 992 | 9220 | 56960 | | |
| $4,9$ | 72 | 256 | 15576 | 50920 | | |
| $4,10$ | 276 | 2260 | 45700 | | | |
| $4,11$ | 240 | 4636 | 51640 | | | |
| $6,6$ | 0 | 64 | 3028 | | | |
| $6,7$ | 0 | 0 | 5712 | | | |
| $4,4,4$ | 0 | 4 | 16 | 24314 | | |
| $4,4,6$ | 0 | 0 | 9768 | | | |
| $4,4,7$ | 0 | 0 | 9600 | | | |
| $4,4,4,4$ | 0 | 0 | 0 | | | |

The $n=3$ entries of $4,10$, $6,6$, $6,7$, $4,4,6$, $4,4,7$ equal the authors' `orders/d4_m8_classes.json` (and the journal); $4,11$ with $n=3$ ($51640$) and $4,4,4,4$
were not yet in the authors' files. $Z_{n,4}$ (one further front of four edges): $1,32,0,12,26,302,288,$ **$1988$**, **$406$** for $n=1,\dots,9$ (also without prune), and
beyond: $Z_{10,4}=3193$, $Z_{11,4}=5911$, $Z_{12,4}=26956$.

## Evaluation (`outputs/repro_practical_d4*.txt`, `outputs/scaling_d4*_repro.txt`, `outputs/ablation_d4*.txt`)

* `practical_beta.py --terms orders/d4_s8.json --closed 10:60,11:66 --closed-from 12 4 4.3e-5 4.5e-5`: $D/\beta^4=48.17$ against $53.83$ (majorant $45.5$, one-plaquette
  junctions $1.99$, exits $0.332$, $\theta$-loss $0.0009$), and $57.61$ against $51.44$ at $4.5\cdot10^{-5}$.
* The same with `orders/d4_s8b.json` at $4.5,5,5.2\cdot10^{-5}$: $30.23<51.44$, **$45.57<46.30$**, $53.12>44.52$; at $5\cdot10^{-5}$ majorant $41.5$, one-plaquette $3.02$, exits $0.499$.
* `scaling_check.py --d4 orders/d4_s8.json 4.3e-5 1`: all hypotheses hold at the 17 values of $t$, $D/\beta^4\le(1+\kappa)^3D_*/\beta_*^4=48.181<53.833$.
  `scaling_check.py --d4 orders/d4_s8b.json 5e-5 1`: all hypotheses hold at the 17 values of $t$, $D/\beta^4\le(1+\kappa)^3D_*/\beta_*^4=45.572<46.296$ (the same numbers as the authors' `orders/d4_scaling.out`). The scaling paragraph of (vii) uses neither $d$ nor $S$; its monotonicity facts are those of
  $SU(3)$ and of $\theta/\beta^6$, $a/\beta$, $b/\beta$, which hold for $K=20$ and the $\theta$ of (ix) term by term (each term of $\theta/\beta^6$ is $\beta^{-6}$ times a power of $\epsilon$, $x$ or $\sigma_6$).
* Which ingredients carry the threshold ($D/\beta^4$; `ablation_d4.py`):

| variant | at $4.3\cdot10^{-5}$, `d4_s8.json` (target $53.83$) | at $5\cdot10^{-5}$, `d4_s8b.json` (target $46.30$) |
|---|---|---|
| as in the paper | 48.18 | 45.57 |
| $\theta$ crude from 10 plaquettes (no exact 10, 11) | 53.65 | **59.1** (not resolved) |
| $\theta$ with 10 exact, crude from 11 | 48.24 | 45.76 |
| $Z_{9,4}$ crude | 61.0 | 72.8 |
| $Z_{8,4},Z_{9,4}$ crude | 1255 | 2255 |
| $N^v_{11}$ crude | 107 | 132 |
| $Z_{3,k}$ crude | 62.5 | 68.1 |
| $S(7)$ instead | 128 | 233 |
| all exact pair counts $\times10$ | – | 45.62 |
| all exact pair counts $\times100$ | 48.68 | 46.15 |
| all exact pair counts $\times10^6$ | 5133 | 5959 |
| $Z_{8,4},Z_{9,4},N^v_{11},Z_{3,k}$ all $\times1000$ | 48.18 | 45.57 |

Each of $Z_{8,4}$, $Z_{9,4}$, $N^v_{11}$, $Z_{3,k}$ and the new pairs is needed, but none needs to be accurate: what matters is that they are exhaustive counts, small, and
that the set $S$ is complete (checked by `assemble_check_d4.py`). At $5\cdot10^{-5}$ the margin is $1.6\%$, and the closed polymers of $10$ and $11$ plaquettes of $\mathbb Z^4$ become load-bearing: counted crudely from $10$ plaquettes, $\theta$ alone costs $13.6\beta^4$. Their numbers ($60$, $66$) are reproduced by `closed_count.py`, and only their smallness matters ($\times1000$ changes nothing); the absence of closed sets of $8$, $9$ plaquettes (Remark (ii)) is essential.

* Corners (`outputs/corner_factors_d4.txt`, `corner_A_scan.txt`). At $4.3\cdot10^{-5}$ with `d4_s8.json` ($\bar b=2^{-9}$, $\bar a=2^{-13}$) the factors are $(3,5)$: $2.45\cdot10^{-13}$,
  $(7,2)$: $1.23\cdot10^{-13}$, $(2,6)$: $1.09\cdot10^{-13}$, $(6,3)$: $5.5\cdot10^{-14}$, $(5,4)$: $2.4\cdot10^{-14}$, so $(3,5)$ set the majorant, as the authors' revision assumes. Adding
  $(3,5)$ and also $(3,\frac{11}2)$ (with the referee's counts, $4,11$ included) gives $45.57$ at $5\cdot10^{-5}$, the same as with $(3,5)$ alone: after the move the corner
  $(2,6)$ dominates ($1.4\cdot10^{-15}$ against $9.3\cdot10^{-16}$ for $(7,2)$ at $(\bar b,\bar a)=(2^{-8},2^{-12})$), i.e. the classes with two plaquettes and fronts of $12$, $13$ edges,
  as the text says.

## For the next step (not used in (ix); not checked against the authors' code)

* **Exits in the layer of $p$.** `zexit_count_n.py` assumes (docstring: "By Lemma `mg:lem:twofront`(i), $p'$ lies in the other layer") that $p'$ is in the other layer, and counts
  each $(P,\varphi')$ once. The lemma allows $\lambda'=\lambda$. For $n\le3$ there are none (`zexit_pfirst.py` tries both layers), so the $Z_{3,k}$ of (ix) are right. With $n=4$ in $d=4$
  there are $8$, all with $|\varphi'|=6$ (`same_layer_exits.py`, `outputs/same_layer_exits_D3_n4.txt`): $q$ together with the three faces of a unit cube at a corner that shares
  exactly one vertex with $q$, $\varphi'$ the skew hexagon bounding those faces. Each such $(P,\varphi')$ gives two pairs $(J,p')$, one per layer. Both layers
  (`outputs/zexit_pfirst_D3_n4.txt`): $Z_{4,k}=152,40,3460,6164,19664,28668,84115,159020,356520,552112,828023,862000,843880,479680,265315$ for $k=6,\dots,20$; a count of the other
  layer only gives $144$ at $k=6$.
* **$N^v_k$ beyond 11 and $Z_{n,4}$ beyond 9** (`fronts_enum.c`, `junc_partition.c`): $N^v_{12}=1121083$, $N^v_{13}=4773186$, $N^v_{14}=31228116$, $N^v_{15}=149251700$ (`fronts_enum 3 15`, 6 min; `outputs/fronts_enum_D3_k15.txt`); $Z_{10,4}=3193$, $Z_{11,4}=5911$,
  $Z_{12,4}=26956$. The crude $N^v_k$, $k\ge12$, also enter the barred tree sums of the majorant, and they dominate it. With `d4_s8b.json` (`bonus_eval_d4.py`,
  `outputs/bonus_eval_d4.txt`):

| $\beta$ | $\beta^3/432$ | as in the paper | with $N^v_{12},N^v_{13}$, $Z_{10..12,4}$ | and $N^v_{14},N^v_{15}$ |
|---|---|---|---|---|
| $5\cdot10^{-5}$ | 46.30 | 45.57 | 14.70 | 9.47 |
| $5.5\cdot10^{-5}$ | 42.09 | 66.20 | 21.25 | 13.60 |
| $6\cdot10^{-5}$ | 38.58 | 93.26 | 29.84 | 19.00 |

  (the majorant drops from $41.5$ to $8.4\beta^4$ at $5\cdot10^{-5}$). The scaling would have to be re-checked at a new point; these counts have not been compared with the authors'.
* Classes still running in the authors' session: $4,11$ with three plaquettes: $51640$; $4,4,4,4$ with two or three plaquettes: none.

## Commands

```
S=<scratch directory>
gcc -O2 -o $S/fronts_enum fronts_enum.c; gcc -O2 -o $S/junc_partition junc_partition.c
$S/fronts_enum 3 11 $S/fronts_D3_k11.txt; $S/fronts_enum 3 13; $S/fronts_enum 2 13 $S/fronts_D2_k13.txt; $S/fronts_enum 2 16
$S/junc_partition $S/fronts_D3_k11.txt 3 6 4,4          # and the other classes of the table; append noprune for the run without the prune
$S/junc_partition $S/fronts_D3_k11.txt 3 12 4           # Z_(n,4), n <= 12
$S/junc_partition $S/fronts_D2_k13.txt 2 4 4,8          # d = 3 validation (as in referee-section16d)
python3 zexit_pfirst.py 3 3; python3 zexit_pfirst.py 2 3
python3 closed_count.py 4 11; python3 closed_count.py 3 12
python3 assemble_check_d4.py
python3 prune_ab.py both 2 4,4:6 4,6:5 ...; python3 prune_ab.py both 3 4,4,4:3 6,6:2 ...
cd ..; python3 practical_beta.py --terms orders/d4_s8.json --closed 10:60,11:66 --closed-from 12 4 4.3e-5 4.5e-5
python3 practical_beta.py --terms orders/d4_s8b.json --closed 10:60,11:66 --closed-from 12 4 4.5e-5 5e-5 5.2e-5
python3 scaling_check.py --d4 orders/d4_s8.json 4.3e-5 1; python3 scaling_check.py --d4 orders/d4_s8b.json 5e-5 1
cd referee-section16f; python3 ablation_d4.py 4.3e-5; python3 ablation_d4.py 5e-5 --terms ../orders/d4_s8b.json
python3 corner_factors_d4.py 4.3e-5; python3 corner_factors_d4.py --scanA 5e-5,5.5e-5,6e-5
python3 zexit_pfirst.py 3 4; python3 same_layer_exits.py 3 4; python3 bonus_eval_d4.py 5e-5 5.2e-5; $S/fronts_enum 3 15
```

Total CPU time of all runs: about two and a half hours, most of it the authors' search without its prune; the referee's C code needs seconds per class.

# Referee checks for Part IV, Section 16 (fourth pass): Remark `mg:rem:nextorder` (vi)–(vii)

Independent checks, made on 13 November, of these parts of `theory/mass-gap-strong-coupling.tex`:

* Remark `mg:rem:nextorder` (vi): the sets $S(m)$, the 22 classes and 34 kinds of pairs of $S(8)$, the exits $Z_{n,k}$ with $3\le n\le5$, the threshold $3.2\cdot10^{-4}$.
* Remark `mg:rem:nextorder` (vii): $S''_9$, its twelve new classes, their counts and its corners; the table; the interpolation between grid points; and the
  paragraph *Scaling*, added during the review (commit `b2f8186`), which now carries the whole range.

These checks are a compass, not part of the proofs.

## Methods

* **Junction counts by a third method.** `junc_faces.c` counts candidate junctions in $d=3$ (slice $\mathbb Z^2$) from Definition `mg:def:junction` and
  Lemma `mg:lem:junction`, by a method that differs from the deficiency-driven searches of the authors (`order_counts.py`, `next_order_count.py`) and of
  the third referee (`referee-section16c/junction_search.py`): the fronts are placed first, and the plaquettes are read off the faces of $W=A\cup B$.
  Two elementary facts organise the search (proofs in the header of the file):
  * (F1) in $\mathbb Z^2$ a link lies in exactly two plaquettes, so by $n_J\ge2$ a link outside $W$ lies in $0$ or $2$ plaquettes of $P$: $P$ is a union
    of bounded faces of $W$;
  * (F2) a cluster has a non-empty mod-2 boundary, which lies in $W$; so in a breadth-first order from $\partial q$ every further front has a vertex within
    $L^\infty$-distance $n$ of the vertices placed before it.

  Every candidate is then checked against all conditions of the definition from scratch ($n_J\ge2$ on $E(J)$, fronts of a layer vertex-disjoint,
  touching, connectivity), and duplicates are removed. Runs with the placement radius $R=n+1$ or $n+2$ instead of $R=n$ give the same counts for all order-9 classes (`outputs/junctions_order9_largerR.txt`). The fronts through
  the origin are dumped by `referee-section16c/fronts_redelmeier.c` (its $N^v_k$ agree with the three earlier enumerations; here also
  $N^v_{14,15,16}=65228,199224,671188$); `junc_faces.c` re-checks every template (connected, all degrees $\ge2$, through the origin).
* **Hypotheses of Theorem `mg:thm:lower5` at given parameters.** `interp_frozen.evaluate` recomputes every condition from the formulas of the paper
  (Definition `mg:def:jconst`, Lemmas `mg:lem:jsums`, `mg:lem:jsums4`, `mg:lem:jsums6`, `mg:lem:corners`, Theorems `mg:thm:lower2`, `mg:thm:lower4`,
  `mg:thm:lower5`): the three conditions on $\xi_F,\xi_C,h_J$ (for $k=4,\dots,59$) and $c_J$, the six inequalities, the barred constants at the fixed
  $(\bar b,\bar a)$, $\Sigma^{(4)}_4<\bar\Upsilon\delta'$, the $F$- and $F_6$-conditions, the third and fourth hypotheses. It reuses the authors' `su3`,
  `dressed`, `counts`, `tree_constants`, `marked` (whose outputs are checked against their definitions: the six inequalities in `evaluate`,
  `barred_check.py`), and the authors' `lower3e` only to find the parameters at the upper point. At the upper point it reproduces `lower3e`'s
  $\Sigma^{(4)}_4$ and $D$ to all printed digits.

## Scripts

| Script | What it does | Result | Runtime |
|---|---|---|---|
| `junc_faces.c` (`DUMP n R k1 k2 ...`) | candidate junctions with $n$ plaquettes and further fronts of sizes $k_i$, fronts first, $P$ from the faces | see below | 0.1 s – 6 min per class |
| `assemble_check.py` | the classes of $S(8)$, of $S''_9$ and of order 9 enumerated directly (not with `order_counts.classes`), the terms built from the referee's counts | $S(8)$: 22 multisets $K$ (43 classes $(n,K)$), 34 kinds of pairs, **equal to `orders/d3_m8.json`**; $S''_9$: 12 multisets with new pairs, 44 terms, **equal to `orders/d3_s9.json`**; 36 multisets of order 9 | 1 s |
| `corner_analysis.py` | corners of any finite $S$ from the definition; the size of each part of $\delta'$ at $\beta_0(3)$ | corners of $S''_9$ as in the paper; the exact terms of order 9 are $1.1\cdot10^{-17}\beta^4$ of $D$; the corner term $2.349\beta^4$ comes from $(8,2)$ | 30 s |
| `alt_sets.py` | $D$ at $\beta_0(3)$ (authors' evaluator) for smaller sets $S$ | $S''_9$: $2.5543$; $S(8)\cup\{(6,3),(6,\frac72),(7,2)\}$ (only the zero pairs of order 9): $4.6995$; adding $(5,4),(5,\frac92)$: $2.5543$; $S(8)$: $21.08$ (all $\times\beta^4$; $\beta^3/432=5.749\beta^4$) | 2 min |
| `interp_frozen.py b0 BETA..` | the interpolation of (vii): all parameters of $\beta_0(3)$ frozen ($c_J$ and the marked constants too), every condition at $\beta<\beta_0(3)$ | all hold at $4.02,4.0,3.95,3.9,3.85,3.8,3.0\cdot10^{-4}$; $D_{\rm frozen}\le D(\beta_0)$; at $\beta_0$ itself the decoration conditions fail by $10^{-32}$ (fixed point reached from below) | 3 min |
| `scaling_frozen.py` | the scaling of (vii) ($\varpi=t^p\varpi_*$, $F,F_6\propto(1+\kappa)t$, $\delta'\propto(1+\kappa)^2t^4$, $c_J$ recomputed) with the referee's evaluator, $t=1$ to $10^{-20}$ | all conditions hold for $t<1$ (at $t=1$ only the $10^{-32}$ rounding above); $D/\beta^4\le(1+\kappa)^3D_*/\beta_*^4=2.5574$, from $2.5564$ ($t=1$) to $2.5247$; smallest margins: $F$-condition $4.0\cdot10^{-4}$ ($t=1$) to $3.0\cdot10^{-2}$ | 4 min |
| `scaling_d4.py` | the same scaling with the authors' `lower3e(fixed=...)` in the other settings that the revised text now covers: $d=4$, Remark `mg:rem:practical`(iii) ($M_5$, $\theta_7$, $\beta_*=2.6\cdot10^{-6}$) and Remark `mg:rem:nextorder`(iii) ($S_6$, $\theta_{10}$, $\beta_*=1.2\cdot10^{-5}$) | all hypotheses hold at $t=1,0.9,0.5,0.1,10^{-3},10^{-6},10^{-12},10^{-20}$ in both settings; $D/\beta^4\le(1+\kappa)^3D_*/\beta_*^4$: $871.395\le871.397$ at $t=1$ ($M_5$; $D_*/\beta_*^4=871.39<890.3$), $187.949\le187.951$ ($S_6$; $<192.9$) | 5 min |
| `monotone_check.py` | the monotonicity facts of the scaling paragraph on 400 points of $[10^{-24},\beta_0(3)]$, 220 digits | $a/\beta$, $b/\beta$, $\epsilon/\beta$, $b_2/\beta^2$, $\sigma_6/\beta^6$, $\theta/\beta^6$, $\log(\bar\Upsilon/\Upsilon)/\beta^4$, $u/\beta$ strictly increasing; $1\le\bar\Upsilon\le1+\kappa$, $\bar\Upsilon(\beta_0)-1=\kappa/6$ | 3 min |
| `theta_sensitivity.py` | $D$ at $\beta_0(3)$ with the closed polymers counted exactly below $15$, $14$, $13$ plaquettes and crudely from there | $2.554$, $2.911$, $7.897$ ($\times\beta^4$; against $5.749$): the 13-plaquette count ($0$) is load-bearing, the 14-plaquette count ($70$) is not | 1 min |
| `barred_check.py` | the barred constants at $(\bar b,\bar a)=(2^{-7},2^{-11})$, $\beta_0(3)$, against their definition with 200-term series | conditions hold with equality up to $10^{-28}$ ($\bar\xi$) and $10^{-23}$ (six inequalities): truncation and fixed point from below; $\bar Y$ at $0.09$, $\bar A$ at $0.34$ of their radii | 30 s |

Junction counts (`outputs/junctions_*.txt`, `outputs/exits_Znk.txt`), each equal to the authors':

* Validation, orders $\le7$: $N_3=48$, $N_2=N'_2=24$, $N_{3,6}=N_{3,7}=24$, $N_{2,8}=48$, $N_{2,9}=16$, $N_4=N_{2,444}=0$; $S'_7$: $276,276,252,0,424,100,212$ and the zeros.
* Order 8 (all 22 multisets, every $n$ needed): $(4,4)$: $48$ at $n=3$, else $0$; $(4,6)$, $(4,7)$: $24,24,276,0$; $(4,8)$: $4,48,252,604$; $(4,9)$: $8,16,424,704$;
  $(4,10)$: $20,100,732$; $(4,11)$: $16,212,856$; $(4,12)$: $12,380$; $(4,13)$: $412$; $(6,6)$: $0,64$; $(6,7)$: $0,128$; $(6,8)$: $48$; $(6,9)$: $60$; $(7,7)$: $0$;
  $(4,4,4)$: $0,0,292$; $(4,4,6)$, $(4,4,7)$: $0,168$; $(4,4,8)$, $(4,4,9)$: $0$; $(4,6,6)$: $20$; $(4,6,7)$: $40$; $(4,4,4,4)$: $0$.
* Order 9 for $S''_9$: $(4,4)$ with $n=5,6,7$: $0$; $(4,6)$, $(4,7)$ with $n=5,6$: $0$; $(6,6)$: $64$, $(6,7)$: $128$ ($n=4$); $(4,4,4)$ with $n=5$: $4$;
  $(4,4,6)$, $(4,4,7)$: $168$ ($n=4$); $(4,8)$: $1152$, $(4,9)$: $2288$ ($n=5$); $(4,10)$: $2432$, $(4,11)$: $4556$ ($n=4$).
* Exits of two front-pairs: $Z_{2,k}=4,4,22,36,42,24,18$ ($k=6..12$); $Z_{3,k}=20,20,46,52,150,260,382,360,288,144,76$ ($k=6..16$);
  $Z_{4,k}=0,0,76,148,328,484,1143,2020,3114,3672,3651$ and $Z_{5,k}=0,0,12,0,376,880,2271,4096,8562,15532,25227$ ($k=6..16$; $k\ge17$ not checked).

Also reproduced with the authors' evaluator (`outputs/practical_m8_repro.txt`): with $S(8)$, $D/\beta^4=6.761$ at $3.2\cdot10^{-4}$ (against $7.234$) and $7.867$ at
$3.3\cdot10^{-4}$ (against $7.015$).

## Commands

```
S=<scratch directory>
gcc -O2 -o $S/fronts_redelmeier ../referee-section16c/fronts_redelmeier.c
$S/fronts_redelmeier 2 13 $S/fronts_D2_k13.txt 4; $S/fronts_redelmeier 2 16 $S/fronts_D2_k16.txt 4     # 3 s, 5 min
gcc -O2 -o $S/junc_faces junc_faces.c
$S/junc_faces $S/fronts_D2_k13.txt 4 4 4 10          # n = 4, R = 4, further fronts 4 and 10: 2432; etc.
$S/junc_faces $S/fronts_D2_k16.txt 3 3 16            # Z_{3,16} = 76
python3 assemble_check.py; python3 corner_analysis.py; python3 alt_sets.py
python3 interp_frozen.py b0 4.02e-4 4.0e-4 3.95e-4 3.9e-4 3.85e-4 3.8e-4 3.0e-4
python3 scaling_frozen.py b0; python3 scaling_d4.py; python3 monotone_check.py 400; python3 barred_check.py; python3 theta_sensitivity.py
```

Total CPU time of all runs: about one hour (the slowest single class, $(4,4,4)$ with five plaquettes and $R=6$, 5.7 min).

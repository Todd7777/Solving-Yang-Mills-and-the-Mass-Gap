# Referee report on Part IV, Section 16 (seventh pass): Lemma 16.10 and Remark 16.7(x), the lower bound in $d=4$ on the whole range

Independent checks, made on 13–14 November, of `theory/mass-gap-strong-coupling.tex`:

* Lemma `mg:lem:frontsdecision` ("Fronts counted by decisions", Lemma 16.10, lines 4843–4869), its constants, and its use in place of
  $(d-1)\binom{(4d-6)k}{k-1}$;
* Remark `mg:rem:nextorder`(x) ("The lower bound in $d=4$ on the whole range", lines 4744–4774) and what depends on it: the abstract (lines 132–136), the
  rows of the summary table for Lemma 16.10 and Remark 16.7(x) (lines 5025–5026), the row of Remark 16.9 (line 5023), Remark `mg:rem:upperexact`(iv).

The claim: for $SU(3)$ in $d=4$, in infinite volume, $-11.4\beta^4\le m(\mu_\beta)-4\log(1/u)+\log\Upsilon\le0.21\beta^4$ for all $\beta<\beta_0(4)$ (a
numerical evaluation; the upper half is Remark 16.9(iv), refereed before and not rechecked here).

These checks are a compass, not part of the proofs. The counts were redone with code written for this pass (fronts, closed polymers, the set
$S^\times_{10}$, a fronts-first count of the classes with a front of $12$ or $13$ edges) and with the sixth referee's partition search
(`../referee-section16f/junc_partition.c`, independent of the authors' code), run on the referee's own list of fronts. The authors' evaluator
(`../practical_beta.py`, `../third_order_effective.py`, `../scaling_check.py`, `../monotone_check.py`) was run to reproduce the evaluation and, with
ingredients removed or inflated, to find what carries it. The authors' search programs (`junctions.c`, `zexits.c`, `closed_sets.c`, `fronts_decision.c`)
were read, not run.

## Recommendation

**Minor revisions.** No error found.
* Lemma 16.10 and its proof are correct, and the constants are right.
* Every count of (x) that was redone agrees exactly with the authors':
  * $N^v_k$ to $k=14$, $Z_{n,4}$ to $n=12$, and the closed polymers of $\mathbb Z^4$ to 16 plaquettes;
  * 30 class runs, by two methods (table below), among them every class that the result needs;
  * the exits with three and four plaquettes (the sixth referee's code, rerun);
  * the set $S^\times_{10}$, with its 73 terms and 8 corners.
* The evaluation at $\beta_0(4)$ ($D/\beta^4=11.383<15.054$) reproduces, as do the quoted grid points, the monotonicity facts in $d=4$ and the scaling
  from $\beta_0(4)$ ($D/\beta^4\le11.388$ at all 17 values of $t$).
* The result is far more robust than the text says. Beyond $S(8)$, only three exponents of orders 9 and 10 are needed: $(6,3)$, $(7,2)$, $(8,2)$. Their
  counts may be wrong by a factor $10^9$, and the scaling with this smaller set passes too ($D/\beta^4\le11.840$). See issue 4.

The issues below concern wording and presentation; none changes a number.

## Issues

1. *(Imprecision; abstract, lines 132–133.)* "with orders $8$ and $9$ and part of order $10$": order 9 is counted only in part (the pairs with $n\ge4$ and the
   exponent $(3,6)$); $(2,7)$ and $(3,\frac{13}2)$ are corners. Say "order $8$ and parts of orders $9$ and $10$", as the summary table (line 5026) in effect does.
2. *(Imprecision; abstract, lines 135–136.)* "…$\le0.21\beta^4$ for all $\beta<\beta_0(4)$, and the same holds (numerical evaluations)": "the same holds" has no
   clear antecedent. Say "so the third-order term is determined from both sides on the whole range there as well".
3. *(Imprecision; (x), first item, line 4746.)* The evaluation does not use $N^v_k\le0.734\cdot7.355^k$ (the lemma at $s=0.194$). For each $k\ge14$ it takes
   the minimum, over $s\in\{\frac1{400},\dots,\frac{399}{400}\}$, of the bound of Lemma 16.10 (`third_order_effective.decision_bound`).
   * This is legitimate: the lemma holds for every $s>0$, and the bound is used only as an upper bound.
   * It makes no difference: with $s=0.194$ fixed, $D/\beta^4=11.383$ as well.
   * One clause would align the text with the code.
4. *(Suggestion, the main one; (x).)* State which steps carry the result, as (ix) does. The referee's ablation at $\beta_0(4)$ (`ablation_x*.out`,
   `essential_x.out`, `scaling_min_set.out`) gives the following.
   * **Needed:**
     * Lemma 16.10: $46.5$ with the binomial bound beyond $13$;
     * $Z_{n,4}$, $n=10..12$: $310$ when counted crudely;
     * the exits with four plaquettes: $39.2$ when crude;
     * the closed polymers of 12 plaquettes: $35.9$ with $\theta$ crude from $12$, while crude from $14$ is enough ($11.43$);
     * of the exponents of orders 9 and 10, only $(6,3)$, $(7,2)$ and $(8,2)$, i.e. the class $4,4,6$ at $n=6$ ($41848$) and the class $4,4,4$ at $n=7,8$
       ($21552$, $22038$); of $S(8)$, at least $(4,4)$, $(5,3)$, $(5,\frac72)$, $(6,2)$ (already in the set of (ix)). Removing $(8,2)$ gives $106$, removing $(7,2)$ gives $1305$, removing $(6,3)$ gives $75.7$.
   * **Sufficient:** $S(8)\cup\{(6,3),(7,2),(8,2)\}$ (45 terms) gives $D/\beta^4=11.835$. The scaling from $\beta_0(4)$ passes at all 17 values of $t$ with it
     ($D/\beta^4\le11.840<15.054$).
   * **Not needed, each alone** ($D/\beta^4\le11.84$):
     * $N^v_{12}$ and $N^v_{13}$ exact: $11.385$ with the decision bound from $k=12$;
     * the exponent $(3,6)$: $11.383$;
     * the order-9 pairs with four plaquettes: $11.475$ without all of them;
     * every other single new exponent (list below).
   * **The values of the counts hardly matter.**
     * All order-9 and order-10 counts $\times10^9$: $11.383$.
     * All pair counts above the lowest order ($n+E\ge\frac{11}2$) $\times10^6$: $11.976$.
     * $Z_{8..12,4}$, $Z_{4,k}$ and the closed-polymer counts $\times100$: unchanged.
     * Only the lowest-order counts $N_3$, $N_2$, $N'_2$ of Lemma 16.3 (a Theorem) carry weight: all pair counts $\times10^3$ gives $29.6$.

   So (x) rests on the enumerations of orders 9 and 10 only through their bookkeeping: every class with a pair at an exponent of $S$ must have *some* count.
   That bookkeeping was checked independently (below). This robustness is worth stating, and the evidence for the large C enumerations then matters much
   less.
5. *(Suggestion.)* At $\beta_0(4)$ the majorant ($8.14\beta^4$ of $D=11.38\beta^4$) is set by the corner $(9,2)$. At $(\bar b,\bar a)=(2^{-8},2^{-11})$ its factor is
   $3.6\cdot10^{-14}$, and the next corner, $(8,3)$, has $7.9\cdot10^{-15}$. The next largest part of $D$ is the crude count of the exits with five plaquettes,
   $1.11\beta^4$ of the exits' $1.78\beta^4$ (`parts_x.out`). These are the places to gain margin. Naming the binding corner would help the reader.
6. *(Imprecision; proof of Lemma 16.10, lines 4860–4868.)* Two small points are used tacitly:
   * $y\ge1$: an edge of $T$ to an open vertex of degree $\ge2$ multiplies $\Phi$ by $1\le y$. This holds since $y\ge(1+s)^{D-1}$.
   * For $N\ge1$, a state with no open vertex has $Z_N=1=\Phi$.

   One clause each.
7. *(Suggestion; Lemma 16.10.)* The proof uses only that the graph is simple with all degrees $\le2D$. So the bound holds verbatim on periodic lattices
   with all sides $\ge3$, where $N^v_k$ is defined (Section 7). This is not needed for (x), which is in infinite volume, but is worth a sentence.
8. *(Imprecision.)* Two places were not updated for (x):
   * the title/date line (line 63: "in $d=3$ to the whole range, from both sides");
   * the introduction of Section 16 (line 4213: "in $d=3$ until the whole range is covered").
9. *(Imprecision; checks.)* (x) cites `monotone_check.py --d4`, but `orders/monotone_check.out` holds only the $d=3$ run; the $d=4$ output is not stored.
   * The referee's run is in `monotone_check_d4.out`. All facts hold with the $\theta$ of (x), and $\bar\Upsilon(\beta_0(4))-1=0.1667\,\beta_0(4)\le\kappa/5$.
   * Cite `orders/d4_x_scaling.out` and `orders/d4_whole_range.out` by name in (x).
10. *(Suggestion; journal 2026-11-13, Part 17.)* Its caveat says "the classes of order 10 have not been counted by a second method". They now have
    been (below), and all counts are equal.

## Lemma 16.10 (fronts counted by decisions)

*Proof.* Checked line by line.
* **The decision process reaches each front through $v$ exactly once.**
  * The state determines the next vertex.
  * The decisions dictated by a front $\varphi$ keep $S\subseteq\varphi$ and $X\cap\varphi=\emptyset$.
  * A decided vertex has all its edges decided. So at termination $S$ contains every edge of $\varphi$ at a vertex of $S$, and $S=\varphi$ by
    connectedness.
  * Conversely, every terminal $S$ is connected with all degrees $\ge2$.
  * Infinite branches are not completions.
* **Potential.** Deciding $w$ removes the factor $y^{\min(\deg_Sw,2)}$. Each edge of $T$ goes to a distinct undecided vertex, since a decided vertex has no
  undecided edge. That vertex is of one of three kinds:
  * new: factor $y$;
  * open of degree one: $y\to y^2$;
  * open of degree $\ge2$: factor $1\le y$.

  The closing edges to open vertices are the second and third kinds, so they are covered. Then:
  * degree one: $T\ne\emptyset$, with at most $2D-1$ undecided edges, so $\sum_{T\ne\emptyset}s^{|T|}\le(1+s)^{2D-1}-1\le y$;
  * degree $\ge2$: at most $2D-2$ edges, so $\sum_Ts^{|T|}\le(1+s)^{2D-2}\le y^2$.

  Checked by brute force for every degree and every number of undecided edges, $D=2,3,4$ (`decision_lemma_check.py`).
* **First step.** Afterwards $\Phi=y^{|T|}$, which gives $\sum_{|T|\ge2}s^{|T|}=(1+s)^{2D}-1-2Ds$.
* **Limit.** $Z_N$ is monotone in $N$, and every completion is finite.

*Constants* (`decision_lemma_check.out`).
* $D=3$, $s=0.194$: $C=0.733517\le0.734$ and $y/s=7.354283\le7.355$. The optimum over $s$ is $7.35284$, at $s=0.1939$.
* $D=2$, $s=0.325$: $C=0.782219\le0.783$ and $y/s=4.080625\le4.081$.
* The bound lies above every exact count: in $\mathbb Z^3$ to $k=14$ the ratio exact/bound is at most $5.6\cdot10^{-3}$ ($3.1\cdot10^{-5}$ at $k=14$); in
  $\mathbb Z^2$ to $k=16$ likewise.
* $(d-1)\binom{(4d-6)k}{k-1}$ grows like $(10^{10}/9^9)^k=25.81^k$, as stated.

*Use.* $N^v_k$ enters, always as an upper bound for the number of fronts with $k$ edges through a vertex and always with non-negative coefficients:
* Definition `mg:def:jconst` ($A_j$, $B_j$);
* Lemma `mg:lem:twofrontcount` ($X_{n,k}$);
* the one-plaquette bound ($\mathsf F,\mathsf A$: "at most $8N^v_k$" fronts through a vertex of $c$, in two layers);
* the barred tree sums.

So any upper bound is admissible, and Lemma 16.10 gives one. In the code (`third_order_effective.lower3e`, `NV_DECISION`) it replaces the binomial for
$k\ge14$ only; $k\le13$ are exact. At $\bar a=2^{-11}$ the barred series converge with ratio about $7.355\sqrt{\bar a}=0.16$ per edge, against $0.57$ with
the binomial. So the truncation at $k<60$ is negligible beside the relative $10^{-15}$ inflation.

## Counts reproduced

Classes are named as in the paper, by the sizes of all their fronts, $\partial q$ included (so $4,4,4$ has the further fronts $4,4$); the counts are by $n$, the number of plaquettes.

| Count | Referee's method | Result (referee = authors in every entry) |
|---|---|---|
| $N^v_k$, $\mathbb Z^3$, $k=4..14$ | `fronts_redel.c`: Redelmeier on edges, each translation class once with its smallest *edge* at the origin (the sixth referee rooted at the smallest vertex; the authors use a decision search); 3.5 min | $12,132,108,2034,4344,43650,147888,$ **$1121083$, $4773186$**, $31228116$; $N_l(k)$ also equal |
| $Z_{n,4}$, $n\le12$ | partition search on the referee's list of fronts | $1,32,0,12,26,302,288,1988,406,$ **$3193,5911,26956$** |
| closed polymers of $\mathbb Z^4$, to 16 plaquettes | `closed_partition.c`: partition search (exclusions never cleared, most constrained link first, each set reached once, 0 duplicates; 3 s) | $6{:}4$, $10{:}60$, $11{:}66$, **$12{:}320$, $13{:}0$, $14{:}1316$, $15{:}2520$, $16{:}12552$**; none with $7$–$9$ |
| closed polymers of $\mathbb Z^3$, to 14 | same | $6{:}2$, $10{:}10$, $11{:}11$, $12{:}24$, $14{:}70$ |
| classes of order 10 with $n\ge5$, and order 9 with $n\ge4$ | partition search | $4,4,4$: $21552$, **$22038$** ($n=7,8$). $4,4,6$: **$41848$**, $60272$ ($n=6,7$). $4,4,7$: $24992$, $80632$ ($n=6,7$). $4,4,4,4$: $1928$, $43742$ ($n=5,6$). $4,4,8$: $258944$, $336752$ ($n=5,6$). $4,4,9$: $459880$, $415216$ ($n=5,6$). $4,4,10$: $266876$, $1616240$ ($n=4,5$). $4,4,11$: $536716$ ($n=4$). $4,6,6$: $7524$, $95900$ ($n=4,5$). $4,4,4,6$: $23352$, $831600$ ($n=4,5$). $4,6,7$: $11208$, $172856$ ($n=4,5$). $4,4,4,7$: $10800$, $778488$ ($n=4,5$). |
| classes with a pair at $(3,6)$; order 8 with $n\le3$ | partition search | $4,6,8$: $1872,20616$. $4,6,9$: $2076$. $4,7,7$: $0$. $4,4,6,6$: $820,6112$. $4,4,6,7$: $1608$. $4,4,4,8$: $0,29664$. $4,4,4,9$: $0$. $4,4,4,4,4$: $0,0$. $4,4,11$: $240,4636,51640$. $4,6,7$: $5712$. $4,4,4,6$: $9768$. $4,4,4,7$: $9600$. |
| classes $4,4,k$ with $n\le2$, $k=9..13$ | `junc_bigfront.py`, fronts first: the big front through $\partial q$, then $h$, then $P$; layers fixed by $\lvert A\triangle B\rvert\le4n$ | $k=9$: $72,256$. $k=10$: $276,2260$. $k=11$: $240,4636$. $k=12$: **$300,14096$**. $k=13$: $0$, **$19644$**. |

Every completed partition-search run (25 runs, `junc_partition_counts.out`, compared by `compare_counts.py`) equals the authors' counts at every $n$, with 0
duplicates and 0 candidates cut by the box.

The partition search could not count $4,4,12$ and $4,4,13$. At closed states its candidate sets hold every placement of a $12$- or $13$-edge front
through a vertex, and the hash sets overflow; with $2^{27}$ slots the $4,4,13$ run was still going after 17 minutes. Hence the fronts-first count. The
class $4,4,12$ with $n=3$ ($133640$, at $(3,6)$, not needed) was not redone.

The exits with three and four plaquettes, $Z_{3,k}$ and $Z_{4,k}$, were not counted by new code:
* the sixth referee's `zexit_pfirst.py` (plaquettes first, $\varphi'=\mathrm{Def}(P)\cup T$, both layers), rerun for this pass (`zexit_pfirst_D3_n4_rerun.out`,
  8.5 min), gives exactly the authors' lists `ZEXITN_D4` ($Z_{4,k}=152,40,3460,\dots,862000,\dots,265315$ for $k=6..20$, of which $8$ in the layer of $p$);
* the evaluation does not depend on their values ($\times100$: $11.383$).

## The set $S^\times_{10}$ (`assemble_check_x.py`)

* **Rebuilt independently.** The set was rebuilt from the definition of $E$ alone, without `order_counts.classes` or `orders/assemble_set.py`: every multiset
  $K$ of further front sizes $4,6,7,\dots,24$ with $2\le|K|\le7$, and every $n\le12$.
* **Bookkeeping.** 22 classes have a pair in $S^\times_{10}$. Each is counted in the merged caches at least as far as needed, and the caches agree where they
  overlap.
* **Terms and corners.** The 73 terms and the corners $(2,7),(3,\frac{13}2),(4,6),(5,\frac{11}2),(6,5),(7,4),(8,3),(9,2)$ equal `orders/d4_s10x.json` exactly.
* **Excluding $(5,\frac{11}2)$ is legitimate.** Lemma `mg:lem:corners` allows any finite $S\subseteq\mathcal E$. Its proof uses only that every exponent
  outside $S$ dominates a minimal one, which is true in $\mathbb N\times\frac12\mathbb N$. So $(5,\frac{11}2)$ simply becomes a corner.
* **Unused counts.** The classes $4,6,7$ and $4,4,4,7$ at $n=5$ are counted, but their only pairs of order $\le10.5$ lie at $(5,\frac{11}2)$. They are
  correctly not used.

## Code reading

* **`junctions.c`** follows `order_counts.count_junctions` step for step:
  * the smallest link of coverage one;
  * the subsets $T$ of the plaquettes through it (the others excluded), times at most one new front per free layer;
  * closed states memoised by their full key, with the exclusions cleared below them;
  * seeds at closed states: plaquettes at a vertex, and fronts of a remaining size through a vertex, in both layers.

  Checks:
  * The keys are exact. The plaquettes are sorted, and fronts of one layer are vertex-disjoint, so (first link, layer) orders them uniquely.
  * The per-size shortcuts are valid. At a seed, $|H\setminus U|\ge k-|U|$. At a deficient link,
    $|D_1|-|H\cap D_1|+|H\setminus U|\ge|D_1|-\min(k,|D_1|)+\max(0,k-|U|)$.
  * A pair of new fronts is pruned one front at a time. This is weaker, but valid.
  * One harmless inefficiency: in the seed loop the recursion can overwrite the de-duplication stamp of the vertex set, so a front may be seeded twice.
    The memo `seen` and the set `found` absorb this.
* **`zexits.c`**:
  * $P$ by Redelmeier rooted at $V(q)$, each set once;
  * fronts inside $G=\partial q\cup E(P)$, each from its smallest edge, by a restricted decision search;
  * all conditions checked;
  * a front vertex-disjoint from $\partial q$ counts in both layers. The coverage formula is the same in both, and in the layer of $p$ the disjointness is
    what Lemma `mg:lem:junction` requires.

  Correct, and it agrees with the sixth referee.
* **`closed_sets.c`**:
  * Only link-sharing extensions: at a deficient link, plaquettes through it; at a closed set, plaquettes sharing a link with it. So every set produced is
    link-connected, and every link-connected closed set is reached.
  * The memo of closed sets is valid, because the expansion below a closed set ignores the exclusions.
  * The prune $\#\{\text{coverage one}\}\le4(n_{\max}-|Y|)$ is valid.
  * The box ($R=n_{\max}+4$) suffices.
* **`fronts_decision.c`**: the decision search with the prune $|S|+\lceil\#\{\deg1\}/2\rceil\le K$. It agrees with `fronts_redel.c` to $k=14$.

## Evaluation, scaling, monotonicity

* **The evaluation at $\beta_0(4)$.** `practical_beta.py --terms orders/d4_s10x.json --closed 10:60,…,16:12552 --closed-from 17 4 1.53766448e-4`
  (`repro_practical_beta0.out`) gives $D/\beta^4=11.38$ against $15.05$, with $(\bar b,\bar a)=(2^{-8},2^{-11})$. The point used is below
  $\beta_0(4)=1.5376644804\cdot10^{-4}$ by $4\cdot10^{-14}$.
* **The quoted grid points.** At $1.25\cdot10^{-4}$ and $1.4\cdot10^{-4}$: $4.02$ and $6.842$ (`repro_grid_points.out`), as quoted.
* **Parts of $D/\beta^4$ at $\beta_0(4)$** (`parts_x.out`):
  * majorant $8.14$;
  * exits $1.78$: one plaquette $0.61$, nearly all $K_Saf(6)$; $Z_{2,k}$ $0.011$; $Z_{3,k}$ and $Z_{4,k}$ below $10^{-4}$; five plaquettes, crude, $1.11$;
  * one-plaquette junctions $1.06$;
  * shifts and jumps $0.25$;
  * two four-edge front-pairs $0.145$, mostly $Z_{5,4}=26$;
  * exact terms of $M_S$ $0.018$, all from $n+E=5$ (Lemma 16.3);
  * $\theta$-loss $4.8\cdot10^{-6}$.
* **Monotonicity** (`monotone_check.py --d4`, `monotone_check_d4.out`; 400 points of $[10^{-24},\beta_0(4)]$, 220 digits).
  * Increasing: $a/\beta$, $b/\beta$, $\epsilon/\beta$, $\gamma/\beta$, $b_2/\beta^2$, $\sigma_5/\beta^5$, $\sigma_6/\beta^6$, $\theta/\beta^6$ (closed polymers 10–16 exact,
    crude from 17), $\log(\bar\Upsilon/\Upsilon)/\beta^4$ and $\bar\Upsilon-a$.
  * Also $\bar\Upsilon-a\ge1$ and $\Upsilon\ge1$.
  * $\bar\Upsilon(\beta_0)-1=0.16667\beta_0\le\kappa/5$.
* **The scaling of (vii) in this setting.** Its argument uses $d$ and $S$ nowhere.
  * Lemma 16.10 enters only through the constant coefficients $N^v_k$, which do not depend on $\beta$.
  * The new $\theta$ enters only through the monotonicity of $\theta/\beta^6$. This holds term by term ($\sigma_6/\beta^6$, $(e\epsilon)^n/\beta^6$ and $x^{17}/\beta^6$
    increase) and is checked above.
  * Orders in $t$: the corner term is $O(t^{n+E})$ with $n+E\ge11$; the exact terms are $O(t^{n+E})$ with $n+E\ge5$; the exits with four plaquettes are
    $O(t^{k+1})$, $k\ge6$.
  * The barred tree sums are affine in $F,F_6$ with non-negative coefficients, so the factor $1+\kappa$ enters linearly, as the argument requires.

  `scaling_check.py --d4 orders/d4_s10x.json --closed … --closed-from 17 1.53766448e-4 1` (`scaling_d4_x_repro.out`): all hypotheses hold at all 17
  values of $t$, from $1$ to $10^{-20}$, with $D/\beta^4\le(1+\kappa)^3D_*/\beta_*^4=11.3879<15.054$. This is identical to `orders/d4_x_scaling.out`. With the
  smaller set $S(8)\cup\{(6,3),(7,2),(8,2)\}$ (`s8_min_set.json`, `scaling_min_set.out`), the bound is $11.8403<15.054$.
* **The grid.** `orders/d4_whole_range.out` was not rerun in full; three of its points were recomputed (above). It is a confirmation only: the
  whole-range claim rests on the scaling.

## Ablation at $\beta_0(4)$ (`ablation_x.py`: `ablation_x_beta0.out`, `ablation_x_extra.out`; `essential_x.py`: `essential_x.out`)

$D/\beta^4$ against $1/(432\beta_0)=15.054$:

| variant | $D/\beta^4$ |
|---|---|
| as in the paper | 11.383 |
| binomial bound beyond 13 (no Lemma 16.10) | 46.45 |
| $N^v_k$ exact to 11, binomial beyond (the setting of (ix)) | 213.8 |
| decision bound from 12 ($N^v_{12},N^v_{13}$ not used) | 11.385 |
| decision bound at fixed $s=0.194$ | 11.383 |
| $N^v_{12},N^v_{13}\times100$ | 11.383 |
| $Z_{10..12,4}$ crude | 309.8 |
| $Z_{10..12,4}\times100$; $Z_{8..12,4}\times100$ | 11.383 |
| $Z_{4,k}$ crude / $\times10$ / $\times100$ | 39.16 / 11.383 / 11.383 |
| $Z_{3,k}\times100$ | 11.387 |
| $\theta$ crude from 12 / 14 / 15 / 16 | 35.87 / 11.433 / 11.385 / 11.383 |
| closed counts 12–16 $\times100$; 10–16 $\times10^4$ | 11.383 |
| $S''_9=S(8)\cup\{n+E\le\frac{19}2,\,n\ge4\}$ (`orders/d4_s9.json`); the same with $(3,6)$ | 105.99; 105.99 |
| $S^\times_{10}$ without $(3,6)$ | 11.383 |
| $S(8)\cup\{n+E\le\frac{21}2,\,n\ge5\}\cup\{(3,6)\}$ (no order-9 pairs with $n=4$) | 11.475 |
| $S^\times_{10}$ without $(8,2)$; without $(7,3),(8,2)$ | 105.99 |
| $S(8)\cup\{(7,2),(8,2)\}$ | 75.65 |
| $S(8)\cup\{(6,3),(7,2),(8,2)\}$ | **11.835** |
| all exact pair counts $\times10$ / $\times10^3$ / $\times10^6$ | 11.546 / 29.56 / 18210 |
| exact pairs with $n+E\ge\frac{11}2$ $\times10^6$ | 11.976 |
| order-10 classes $\times10^3$; order-9 classes $\times10^3$; both $\times10^9$ | 11.383; 11.383; 11.383 |
| order-$\le8$ classes $\times10^3$ | 29.56 |

Single exponents removed from $S^\times_{10}$ (those with $n+E\ge8$; `essential_x.out`):
* **needed:** $(4,4)$ $28.8$; $(5,3)$ $923$; $(5,\frac72)$ $46.3$; $(6,2)$ $16537$; $(6,3)$ $75.7$; $(7,2)$ $1305$; $(8,2)$ $106$;
* **not needed** (each alone, at most $11.84$): $(2,6)$, $(2,\frac{13}2)$, $(3,5)$, $(3,\frac{11}2)$, $(3,6)$, $(4,\frac92)$, $(4,5)$, $(4,\frac{11}2)$, $(5,4)$, $(5,\frac92)$,
  $(5,5)$, $(6,\frac72)$, $(6,4)$, $(6,\frac92)$, $(7,3)$, $(7,\frac72)$.

The binding pairs are those with small $E$ and many plaquettes, since $a/\bar a=0.0175\ll b/\bar b=0.079$ at $(2^{-8},2^{-11})$. Removing exponents is done
by Lemma `mg:lem:corners` (the corners recomputed); the evaluator then reoptimises $(\bar b,\bar a)$.

## What I could not check

* The class $4,4,12$ with three plaquettes ($133640$, at $(3,6)$) by a second method. It is not needed: $(3,6)$ can be dropped.
* The exits with four plaquettes by a method of my own (the authors' and the sixth referee's agree, and the result does not depend on their values).
* The full 20-point grid `whole_range_d4.py`; three of its points were recomputed.
* The upper half ($0.21\beta^4$, Remark 16.9(iv)), refereed before.
* A proof of the monotonicity facts. Like the authors, I checked them on a grid of 400 points, which is what (vii) claims.

## Files

| File | What |
|---|---|
| `fronts_redel.c` | $N^v_k$ and $N_l(k)$ by Redelmeier on edges, with an optional dump of the classes (the input of the junction searches); `fronts_redel_D3_k14.out` |
| `closed_partition.c` | closed link-connected polymers through a plaquette, by a partition search; `closed_partition_d4_16.out` ($\mathbb Z^4$ to 16, $\mathbb Z^3$ to 14) |
| `junc_partition_big.c` | the sixth referee's `junc_partition.c`, unchanged except for larger hash sets (credit: sixth referee); `junc_partition_counts.out` |
| `junc_bigfront.py` | the classes $4,4,k$ ($k\ge9$) with $n\le2$, fronts first; validated on $k=9,10,11$ (and $k=12$ by both of its versions); `bigfront.out` |
| `zexit_pfirst_D3_n4_rerun.out` | the sixth referee's exit counter, rerun: $Z_{1..4,k}$ in $d=4$ |
| `compare_counts.py` | the partition-search counts against the authors' caches |
| `assemble_check_x.py` | $S^\times_{10}$, its corners and terms, from the definition of $E$; the class bookkeeping |
| `decision_lemma_check.py` | the constants of Lemma 16.10, the one-step potential inequality, the bound against exact counts; `decision_lemma_check.out` |
| `ablation_x.py`, `essential_x.py`, `parts_x.py` | ablations and the parts of $\delta'$, with the authors' evaluator; `s8_min_set.json`, the smaller sufficient set |
| `repro_practical_beta0.out`, `repro_grid_points.out`, `monotone_check_d4.out`, `scaling_d4_x_repro.out`, `scaling_min_set.out` | the authors' evaluator, rerun |

Commands (`S` is a scratch directory):

```
gcc -O2 -o $S/fronts_redel fronts_redel.c; $S/fronts_redel 3 14; $S/fronts_redel 3 13 $S/fronts_D3_k13.txt
gcc -O2 -o $S/closed_partition closed_partition.c; $S/closed_partition 4 16; $S/closed_partition 3 14
gcc -O2 -o $S/junc_partition_big junc_partition_big.c
$S/junc_partition_big $S/fronts_D3_k13.txt 3 12 4         # Z_(n,4); and e.g. 3 8 4,4   3 7 4,6   3 5 4,10   3 4 4,11   3 3 4,6,6 ...
python3 junc_bigfront.py $S/fronts_D3_k13.txt 12           # and 9, 10, 11, 13
python3 compare_counts.py junc_partition_counts.out
python3 assemble_check_x.py; python3 decision_lemma_check.py
REFSCR=$S python3 ablation_x.py; REFSCR=$S python3 essential_x.py; python3 parts_x.py
cd ..; python3 practical_beta.py --terms orders/d4_s10x.json --closed 10:60,11:66,12:320,13:0,14:1316,15:2520,16:12552 --closed-from 17 4 1.53766448e-4 1.25e-4 1.4e-4
python3 monotone_check.py --d4
python3 scaling_check.py --d4 orders/d4_s10x.json --closed 10:60,11:66,12:320,13:0,14:1316,15:2520,16:12552 --closed-from 17 1.53766448e-4 1
python3 scaling_check.py --d4 referee-section16g/s8_min_set.json --closed 10:60,11:66,12:320,13:0,14:1316,15:2520,16:12552 --closed-from 17 1.53766448e-4 1
```

CPU time of all runs: about two and a half hours, on at most two cores (the slowest: the Python fronts-first counts, 20 min, the ablations, 25 min, and a stopped partition-search run of $4,4,13$, 17 min).

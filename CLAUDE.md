# CLAUDE.md

Guidance for Claude Code sessions in this repository, local or cloud.

## What this is

A long-horizon research program on Yang–Mills existence and the mass gap: the precise statement and its
lattice reduction, a survey, a research program, rigorous strong-coupling monographs in `theory/`, and a
tested SU(2)/SU(3) lattice Monte Carlo toolkit in `lattice/`. The problem is open, and this repository does
not solve it.

## Ground rules

- Follow the README's **Standards** section: every claim labelled Theorem, Cited, Conjecture, Heuristic or
  Numerical, the U(1) test, the strong-coupling test, numerics as a compass and never as evidence, and no
  citation without verification.
- Do not claim or imply that the problem is solved.
- End every session with an entry in `journal/`: what was established, what failed, what is next.
- Commit verified work in small steps and push before the session ends. Work left only inside a cloud
  session can be lost with it.

## Setup and checks

```bash
pip install -e "lattice[test,plot]"
cd lattice && pytest -q      # 2026-11-09: 41 pass, 4 slow tests skipped (add --runslow to run them)
```

To build the monographs, first install the TeX packages listed in `.github/workflows/ci.yml` (including
`lmodern`, which the cloud container lacked on 2026-11-09) if `latexmk` or `lmodern.sty` is missing, then:

```bash
cd theory && for f in lattice-foundations infinite-volume-strong-coupling area-law-strong-coupling mass-gap-strong-coupling; do latexmk -pdf -interaction=nonstopmode $f.tex; done
```

The PDFs in `theory/` are committed. Rebuild and commit them whenever their `.tex` source changes.

## State on 2026-11-14

- `theory/` holds four monographs: the foundations [F] (`lattice-foundations.tex`, errata 2026-11-09); Part II
  (`infinite-volume-strong-coupling.tex`, the infinite-volume theory, milestone T1.1; refereed and revised 2026-11-09);
  Part III (`area-law-strong-coupling.tex`, the area law; two referees, revised 2026-11-09); Part IV
  (`mass-gap-strong-coupling.tex`, the strong-coupling gap is $4\log(1/\beta)+O(1)$: $4\log(1/q)\le m\le4\log(1/u)+O(\beta)$ for all
  $\beta<\beta_0(d)$ when $\dim\rho\le4$ ($\le5$ if self-conjugate); new on 2026-11-09; the upper bound was reproved by transfer-matrix positivity on 2026-11-10, refereed the same day, with checks in
  `theory/checks/upper-bound/`; the lower bound sharpened the same day to $m=4\log(1/u)+O(\beta)$ by transfer-kernel norms and counting by fronts, Sections 6–7, checks in
  `theory/checks/lower-bound/`; Section 7 refereed, then revised to be effective on the whole range (for $SU(2)$,
  $-1.38\le m-4\log(1/u)\le0.0014$ in $d=4$), revision refereed (minor revisions, applied); the asymptotics is classical (Schor), see Part IV's Literature paragraph).
  Section 8 (junctions, 10 November) gives $m\ge4\log(1/u)-2\beta+O(\beta^2)$; refereed (minor revisions, applied). Section 9 (the cap by
  characters, 10 November) gives the first-order term: $m=4\log(1/u)-\kappa(\beta)+O(\beta^2)$, $\kappa=\tau_0^{-1}\int(\mathrm{Re}\,\chi_\rho)^2f_\beta$,
  which is $O(\beta^2)$ for $SU(2)$, $U(1)$, $SU(N\ne3)$ and $\beta/6+O(\beta^2)$ for $SU(3)$ (three fundamental fluxes can meet at a link); for $SU(2)$
  $-1.3\cdot10^{-4}\le m-4\log(1/u)\le4.8\cdot10^{-4}$ on the whole range in $d=4$; checks in `theory/checks/first-order/`; refereed (minor revisions, applied). Section 10 (11 November) gives the
  second-order term: $m=4\log(1/u)-\log\Upsilon(\beta)+O(\beta^3)$, $\Upsilon$ the largest variance of a real unit plaquette function in
  $\mathrm{span}(\chi_\rho,\bar\chi_\rho)$ under $w_\beta/c_0$; for $SU(2)$ $m=4\log(1/u)+\beta^2/8+O(\beta^3)$ (the shifts and jumps of the tube vanish on
  the swap-odd sector of the duplicated system); checks in `theory/checks/second-order/`; refereed: minor revisions, applied. Section 11 (11 November) shows that nothing new enters at
  third order: $m=4\log(1/u)-\log\Upsilon+O(\beta^4)$, so the $\beta^3$ term is that of $-\log\Upsilon$ ($\beta^3/432$ for $SU(3)$, $-\beta^3/6000$ for
  $SU(5)$, $0$ for $U(1)$, $SU(2)$, $SU(4)$, $SU(N\ge6)$); asymptotic only (constant $\sim10^{11}$); checks in `theory/checks/third-order/`; refereed:
  minor revisions, applied. Section 12 (12 November) on the fourth order: the plaquette state moves (bulges and hops through six-edge fronts,
  turns and jumps through unit cubes), with a hopping matrix whose top is $h_+=26d-70$ ($h_-=10d-22$ for $C$-odd). Proved: the upper bound
  $\limsup\beta^{-4}(m-4\log(1/u)+\log\Upsilon)\le-hu_1^4$ (e.g. $-\frac{17}{128}$ for $SU(2)$, $d=4$). Conjectured: equality,
  $m=4\log(1/u)-\log\Upsilon-hu^4+O(\beta^5)$. Evidence: exact $\mathbb Z_2$/$\mathbb Z_3$ spectra, and the classical series (Münster, Seo, as quoted by
  Langelage–Münster–Philipsen) at order $u^4$ in four of five channels ($-\frac{98}3u^4$ for the $SU(2)$ scalar; the $SU(3)$ $T_1^{+-}$ differs by $u^4/4$,
  unresolved). Checks in `theory/checks/fourth-order/`; refereed: minor revisions, applied. Section 13 (12 November) identifies the junctions of order
  $\beta^4$: the shifts and jumps, the turns, jumps and U-turns through unit cubes ($Z_{4,4}\le12(d-3)$, proved via Lemmas 5.4 and 12.5), and the widenings
  to six-edge fronts, which end at leading order in one of two narrowings (Lemma 13.3). With $\varpi=b^{5/2}$ the shifts and jumps drop out, and the U-turns
  (downward exits, cost $\eta=O(\beta^{3/2})$) too. So the third-order lower bound holds with $\delta'=(\delta'_4+o(1))\beta^4$, $\delta'_4=160(d-3)+8K_Su_1^2$: the moves of the
  hopping matrix ($h_+=2K_S+10(d-3)$) at the sup norm. That is $166$ for $SU(2)$, $d=4$, instead of $8.4\cdot10^{10}$, and $2$ in $d=3$. With Corollary 12.9 the
  $\liminf$ and $\limsup$ of $\beta^{-4}(m-4\log(1/u)+\log\Upsilon)$ lie in $[-166,-\frac{17}{128}]$ (Theorem). Still asymptotic at practical $\beta$: the crude $\beta^5$
  terms dominate. Checks in `theory/checks/third-order/` (with the referee's independent code in `referee-section13/`); two referees: minor revisions, applied.
  Section 14 (12 November) evaluates the moves by characters. Their junctions are closed surfaces (each link in exactly two faces), on which the transfer
  operator is diagonal in the representations with norm exactly $u_*$ per plaquette (Lemma 14.1: $u_*^6$ widening/narrowing, $u_*^8$ turn/jump). So $b\to a$ for them,
  and $\delta'_4=h_+u_1^4$ (Theorem 14.2). The recursion follows the moving tube because the induction's bound is uniform over all front-pairs; each move
  is counted in absolute value, which gives the row sums of $H_+$. Corollary 14.3, **the fourth-order term**: $m=4\log(1/u)-\log\Upsilon-(26d-70)u^4+o(\beta^4)$:
  * for $SU(N)$, $2\le N\le5$, and for every self-conjugate $\rho$, in every $d\ge3$;
  * $-8u^4$ for every $\rho$ in $d=3$.

  This is Conjecture 12.4 there, with $o(\beta^4)$ in place of $O(\beta^5)$; for $SU(2)$, $d=4$, the coefficient is $-\frac{98}3u^4$, the classical one. For $U(1)$, $d\ge4$, the
  coefficient lies between $-(26d-70)/16$ and $-(10d-22)/16$, because norms lose the $C$-odd signs. Checks: `theory/checks/fourth-order/moves_by_characters.py`
  ($\mathbb Z_3$, $S_3$, exact), and the referee's own code in `referee-section14/` (five finite groups; the full operator from the definitions for $\mathbb Z_2$,
  $\mathbb Z_3$). Refereed: minor revisions, applied.
  Section 15 (12 November), the sectors of charge conjugation. With a charge conjugation ($U(1)$, $U(N)$, $SU(N)$, abelian), $m=\min(m_+,m_-)$, and each sector
  rate is $4\log(1/u)-\log\Upsilon_C+O(\beta^4)$, the even one to $o(\beta^4)$; for $SU(3)$, $m_--m_+=\beta/3+O(\beta^2)$. Lemma 15.6: $\Upsilon_+-\Upsilon_-=\beta^4/248832+O(\beta^5)$
  for $SU(6)$ and $O(\beta^5)$ for $SU(N\ge7)$, so the fourth-order term holds for $SU(N)$, every $N$, every $d\ge3$. The odd sector at fourth order ($U(1)$, $d\ge4$)
  is Conjecture 15.11. The hop/turn cancellation is proved as an operator identity ($Y^0_{\rm turn}=Y_{\rm hop}\mathsf C$, Lemma 15.8) with a partner lemma, but the
  grouping of sets along the path is open. A first draft's proof had major gaps (referee: non-closed grouping, unequal $r_X^2$ factors), which is why it is a
  conjecture. Checks in `theory/checks/fourth-order/` (`c_odd_pairing.py`, `upsilon_split.py`, `referee-section15/`). Second referee of the revision: minor revisions
  (a gap in the passage to $\mu_\beta$, repaired by defining the sector gaps spectrally via the ground-state chain), applied (journal 2026-11-12, Part 11). Next: the odd sector (routes in Remark 15.12), the $u^6$ term, effective constants, near $\beta_0$.
  Section 16 (12 November), closed polymers and the third-order term at practical $\beta$. Closed sets of plaquettes have 6 (a unit cube) or $\ge8$ plaquettes
  (Lemma 16.1), so the ratio bound holds with $\theta_7=2(d-2)e^6\sigma_6+x^8/(1-x)=O(\beta^6)$ in every section (Lemma 16.2). The junctions with three front-pairs of
  lowest order are classified (three caps or a unit cube; a widening with a cap, or a cube with a bent domino) and counted exactly in the junction sum (Lemmas 16.3–16.4).
  With exact $N^v_k$, $Z_{n,4}$, $Z_{2,k}$ (enumerations), the lower bound resolves the $SU(3)$ third-order term $\beta^3/432$ for $\beta\le2.6\cdot10^{-6}$ ($d=4$; before
  $2\cdot10^{-9}$) and $\beta\le2.3\cdot10^{-5}$ ($d=3$) (Remark 16.5, Numerical; the criterion is the whole deficit $8(K+1)\theta+\log(\bar\Upsilon/\Upsilon)+\log(1+\delta')$, and an
  earlier $2.6\cdot10^{-5}$ for $d=3$ had left out $\theta_7$). Lemma 16.6 (13 November): junction sums by corners, any order; Remark 16.7: the next order by exhaustive
  search, and $\theta_{10}$ (no closed sets of 8, 9 plaquettes), give $\beta\le1.2\cdot10^{-5}$ ($d=4$) and $\beta\le9\cdot10^{-5}$ ($d=3$, about $\beta_0(3)/4.5$); one order more in $d=3$, with the closed sets of 10–12 plaquettes
  counted exactly and $Z_{n,4}=0$ for $n\ge3$ (a planar argument), $\beta\le1.8\cdot10^{-4}$; two orders more (classes generated automatically, `order_counts.py`; exact exits
  up to 5 plaquettes; closed sets up to 14), $\beta\le3.2\cdot10^{-4}\approx\beta_0(3)/1.26$ (Numerical). Third referee (13 November; Lemma 16.6, Remark 16.7(i)–(v)): minor revisions, applied. With order 9 in part
  ($S''_9=S(8)\cup\{n+E\le\frac{19}2,\ n\ge4\}$; the twelve classes with such pairs counted, `orders/assemble_d3_s9.py`), **the whole range in $d=3$**: the lower bound resolves the $SU(3)$
  third-order term for all $0<\beta<\beta_0(3)$ (Remark 16.7(vii), Numerical: $D\le2.56\beta^4$ at $\beta_0(3)$, and a scaling of the parameters in $\beta$ ($\varpi\propto t^p$, $F,F_6\propto t$,
  $\delta'\propto t^4$) gives $D(\beta)/\beta^4\le(1+\beta_*)^3D(\beta_*)/\beta_*^4$ below any computed $\beta_*$, `scaling_check.py`; a 22-point grid confirms it, `whole_range_d3.py`). The
  scaling also justifies the monotonicity of $D/\beta^4$ that the earlier thresholds used. Fourth referee (Remark 16.7(vi)–(vii) with the scaling): minor revisions, no error, applied (`referee-section16d/`). Checks in `theory/checks/third-order/`
  (`three_front_count.py`, `practical_beta.py`, `order_counts.py`, `whole_range_d3.py`, ...). Two referees (Lemmas 16.1–16.2 and Remark 16.5; Lemmas 16.3–16.4, with
  independent searches in $d=3,4,5$): minor revisions, applied. The upper bound (Theorem 10.2) saw the term only for $\beta\le3\cdot10^{-8}$ (crude tails $\mathcal T_Q,\mathcal T_p$ of Lemma `mg:lem:twopoint`,
  Remark 16.7(viii)); Lemma 16.8 (Theorem: $\theta_5$ in the tails, exact counts of the admissible sets, the five-face caps by characters, $\sigma_5=\sum d_r^2u_r^5$) and Remark 16.9
  (Numerical: counts by `rel_closed_count.py`, evaluation by `upper_exact.py`) give **for $SU(3)$ in $d=3$: $-2.56\beta^4\le m(\mu_\beta)-4\log(1/u)+\log\Upsilon\le0.042\beta^4$ for all
  $\beta<\beta_0(3)$**, the third-order term from both sides on the whole range; in $d=4$ the upper side on the whole range ($\le0.21\beta^4$, counts to 10 plaquettes) and both
  sides for $\beta\le5\cdot10^{-5}\approx\beta_0(4)/3.1$ (Remark 16.7(ix): order 7 complete and order 8 for $n\ge4$ and at $(3,5)$ in $d=4$, exact $Z_{8,4}=1988$, $Z_{9,4}=406$, $N^v_{11}=147888$, exits with three
  plaquettes; the junction search now has a capacity prune; sixth referee: minor revisions, no error, applied, `referee-section16f/`). **Then both sides on the whole range in $d=4$ too**
  (Remark 16.7(x), 13 November; seventh referee: minor revisions, no error, applied, `referee-section16g/`): Lemma 16.10 (Theorem, fronts counted by decisions: $N^v_k\le0.734\cdot7.355^k$ in $\mathbb Z^3$, against $25.8^k$; the crude $N^v_k$ had
  dominated the majorant, as the sixth referee found), searches in C (`junctions.c`, `zexits.c`, `closed_sets.c`, `fronts_decision.c`; `junctions_validate.py` reruns 37 Python counts),
  $S(8)$, order 9 with $n\ge4$, $(3,6)$, order 10 with $n\ge5$ (`orders/d4_s10x.json`, `orders/assemble_set.py`), $Z_{n,4}$ to 12, exits with four plaquettes, closed polymers of $\mathbb Z^4$ to 16
  plaquettes ($320,0,1316,2520,12552$ for $12..16$); $D\le11.39\beta^4$ for all $\beta\le\beta_0(4)$ by the scaling from $\beta_0(4)$: **$-11.4\beta^4\le m-4\log(1/u)+\log\Upsilon\le0.21\beta^4$, $SU(3)$, $d=4$**
  (Numerical). Robust (the referee's ablation, rechecked by `ablation_x_recheck.py`): beyond $S(8)$ only $(6,3),(7,2),(8,2)$ are needed ($S(8)\cup\{(6,3),(7,2),(8,2)\}$: $11.84\beta^4$,
  the scaling passes), and the order-9/10 counts $\times10^9$ change nothing; needed are Lemma 16.10, $Z_{10..12,4}$, the four-plaquette exits, the 12-plaquette closed polymers. `practical_beta.py --ix` reproduces (ix). Fifth referee (Lemma 16.8, Remark 16.9): minor revisions, no error, applied (`referee-section16e/`). Remark 16.7(iv) now also proves (planar faces) that in $d=3$ pairs with $E=2$ have $n=3$ and pairs with $E\in\{3,\frac72\}$ have $n\le4$.
  **$SU(5)$** (Lemma 16.11 and Remark 16.12, 14 November; eighth referee: minor revisions, applied, incl. an error in Lemma 16.11 for $N=4$ ($\Lambda^2$ self-conjugate), fixed, no effect on $SU(5)$; `referee-section16h/`): the only other $SU(N)$ with a third-order term ($-\beta^3/6000$, Corollary 11.6). Group data from exact
  Haar moments (`sun_group.py`, Schur–Weyl; agrees with the $SU(2)$, $SU(3)$ code to 40 digits); Lemma 16.11 (Theorem) bounds $u_2$ by the exact order-2 multiplicities,
  $u_2\le\beta^2/2000+R/100$ for $SU(5)$ instead of $e^\beta-1-\beta$ (needed in $d=3$: otherwise $\log(\bar\Upsilon/\Upsilon)=0.254\beta^4$). Lower bound: $d=3$ with $S''_9$, $D\le0.3454\beta^4$
  by the scaling; $d=4$ with $S^\times_{10}\cup\{(9,2)\}$ (class $4,4,4$ with nine plaquettes: $148894$, two methods), $D\le0.519\beta^4$. Upper bound: $E=0.318\beta^4$ ($d=3$); $0.2511\beta^4$ ($d=4$), with
  the admissible sets of 11 plaquettes in $\mathbb Z^4$ ($760$ for $S_A$, $2616$ for $S_Q$; `admissible_split.c`, the fifth referee's enumeration in resumable pieces,
  `run_admissible_pieces.sh`, about 1 hour on 4 cores; `orders/admissible_d4_n11.out`). **For $SU(5)$: $-0.346\beta^4\le m-4\log(1/u)+\log\Upsilon\le0.318\beta^4$ ($d=3$), $-0.519\beta^4\le\dots\le0.252\beta^4$ ($d=4$), whole range**
  (Numerical); the `--group 5` option of `practical_beta.py`, `scaling_check.py`, `monotone_check.py`, `upper_exact.py`. The corner $(9,2)$ and the sets of 11 plaquettes also give $SU(3)$, $d=4$: $-3.94\beta^4\le m-4\log(1/u)+\log\Upsilon\le0.033\beta^4$ (was $-11.4$, $0.21$).
  The $d=4$ upper half does not rest on the enumeration: the binomial bound of Lemma 5.5's proof as $N_{11}$ already gives $0.360\beta^4$ (the referee).
  **Computer-assisted proof** (Lemmas 16.13–16.17, Theorem 16.18, Remark 16.19, 14 November): with the group data replaced by proven bounds ($\bar u$, $u_{\rm low}$, $u_2$ of Lemma 16.11,
  $e^{-2\beta}\le\Upsilon\le e^\beta$, $\bar\Upsilon\le\Upsilon+\Delta$; Lemma 16.14 "bound mode", `rigorous_bounds.py`, `--bounds`), every monotonicity fact of the scaling holds by inspection
  (Lemma 16.15, no grid checks); the planar argument in $d=3$ is Lemma 16.16, the upper bound with exact counts Lemma 16.17; only counts made by two programs (`--certified`);
  $\Upsilon\ge1$ for $SU(5)$ by an exact series (Lemma 16.13(ii)). **Theorem 16.18: for $SU(3)$ and $SU(5)$, $d=3,4$, $-c_-\beta^4\le m-4\log(1/u)+\log\Upsilon\le c_+\beta^4$ for all
  $\beta<\beta_0(d)$, $(c_-,c_+)=(2.558,0.042)$, $(3.937,0.047)$, $(0.346,0.318)$, $(0.520,0.361)$**: the third-order term from both sides on the whole range, proved (computer-assisted).
  The evaluation at $\beta_*=\beta_0(d)(1+10^{-30})$ is done twice: ours (mpmath, 60 digits, tails $\le1.4\cdot10^{-21}$ against room $\ge2.8\cdot10^{-16}$, `error_budget.py`) and the ninth referee's in
  Arb ball arithmetic with rigorous remainders (`referee-section16i/arb_eval.py`, `upper_arb.py`). Ninth referee: minor revisions, applied; one error (the marked constants were all raised by one
  factor, leaving the two quadratic inequalities of Definition 8.9 short by $\sim10^{-16}$; fixed in `third_order_constants.marked`, now checked explicitly, no constant changes); the arguments
  that sat in Remarks became Lemmas 16.14–16.17; `run_cap.sh` reruns everything. The summary table of Part IV is now a longtable (its last part had been cut off in the PDF).
  **The two plaquette glueballs** (Corollary 16.20, 14 November; computer-assisted; tenth referee: minor revisions, no error affecting a result, applied, with every constant
  certified in ball arithmetic, `referee-section16j/`): the same bounds for each sector of charge conjugation, for $SU(3)$, $SU(4)$, $SU(5)$, $d=3,4$:
  $-c'_-\beta^4\le m_C-4\log(1/u)+\log\Upsilon_C\le c'_+\beta^4$ for $C=\pm$, all $\beta<\beta_0(d)$, $(c'_-,c'_+)=(0.600,0.042)$, $(3.937,0.047)$ [$SU(3)$], $(0.134,0.131)$, $(1.083,0.148)$ [$SU(4)$],
  $(0.073,0.319)$, $(0.520,0.361)$ [$SU(5)$]. So the odd glueball of $SU(3)$ has the third-order term $\beta^3/1296$, resolved on the whole range; the splitting $m_--m_+$ is
  $2(\beta/2N)^{N-2}/(N-2)!$ to leading order ($\frac\beta3-\frac{\beta^3}{648}$, $\frac{\beta^2}{64}$, $\frac{\beta^3}{3000}$), determined on the whole range, with the even glueball the lighter one; and in
  $d=3$ the empty kinds of Lemma 16.16(ii) lower $c_-$ of Theorem 16.18 to $0.600$ ($SU(3)$) and $0.073$ ($SU(5)$). Checks: `run_cap.sh 2 corollary`, `sector_upper.py`, `rigorous_bounds.py`
  (`sector_series`, `split_lower`), `orders/d3_s9z.json`.
  Theorem 15.5(iii) now states the fixed-$\beta$ sector bound, and the proof of Theorem 10.2 notes that any eigenvector of $\mathsf V$ may be used. For $SU(N)$, $N\ge6$, $d\ge4$, the glueballs split
  at order $\beta^4$ through the hopping too ($h_+\ne h_-$), so "order $\beta^{N-2}$" is claimed only for $N\le5$.
  Next: $d\ge5$ (a crude probe gives $D\approx3\cdot10^{11}\beta^4$ at $\beta_0(5)$: the whole counting programme again), or the $C$-odd sector at fourth order (Conjecture 15.11).
- Glueball runs: raw `.npz` files and checkpoints in `results/` are gitignored and do not survive the container; only
  the analyses (`results/*.md`) are committed. `glueball.py` resumes from its checkpoint when rerun with the same
  arguments.
- The newest journal entry is `journal/2026-11-14.md` (Part 18: the seventh referee, on Lemma 16.10 and Remark 16.7(x); Part 19: the third-order term of $SU(5)$; Part 20: the computer-assisted proof; Part 21: the ninth referee; Part 22: the two plaquette glueballs; Part 23: the tenth referee). Before it, `journal/2026-11-13.md` (Part 13: the next order of the third-order term at practical $\beta$; Part 14: further orders in $d=3$, the whole range; Part 15: the upper bound; Part 16: the lower bound in $d=4$; Part 17: the whole range in $d=4$). Before it, `journal/2026-11-12.md` (session 3, Part 8: the fourth-order term; Part 9: the third-order term with an effective constant; Part 11: the sectors of charge conjugation;
  Part 12: the third-order term at practical $\beta$; Part 10: the lower bound to fourth order); `journal/2026-11-11.md` has Parts 6–7 (the second- and
  third-order terms), `journal/2026-11-10.md` has Parts 1–5,
  `journal/2026-11-09b.md` is the first part of session 3 and `journal/2026-11-09.md` is session 2.

## Where the next steps are

`docs/03-research-program.md` and the newest entry in `journal/`.

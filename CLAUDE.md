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

## State on 2026-11-13

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
  third-order term for all $10^{-12}\le\beta<\beta_0(3)$ (Remark 16.7(vii), Numerical: 22 grid points and an interpolation argument with the parameters of the upper neighbour
  fixed, `whole_range_d3.py`; below $10^{-12}$ only asymptotically, by Theorem 14.2). Orders 8–9 are not yet refereed. Checks in `theory/checks/third-order/`
  (`three_front_count.py`, `practical_beta.py`, `order_counts.py`, `whole_range_d3.py`, ...). Two referees (Lemmas 16.1–16.2 and Remark 16.5; Lemmas 16.3–16.4, with
  independent searches in $d=3,4,5$): minor revisions, applied. Next: a referee for Remark 16.7(vi)–(vii); an effective argument below $10^{-12}$; $d=4$ (order 7 is being
  counted; then the junctions with two four-edge front-pairs and $\ge8$ plaquettes).
- Glueball runs: raw `.npz` files and checkpoints in `results/` are gitignored and do not survive the container; only
  the analyses (`results/*.md`) are committed. `glueball.py` resumes from its checkpoint when rerun with the same
  arguments.
- The newest journal entry is `journal/2026-11-13.md` (Part 13: the next order of the third-order term at practical $\beta$; Part 14: further orders in $d=3$, the whole range). Before it, `journal/2026-11-12.md` (session 3, Part 8: the fourth-order term; Part 9: the third-order term with an effective constant; Part 11: the sectors of charge conjugation;
  Part 12: the third-order term at practical $\beta$; Part 10: the lower bound to fourth order); `journal/2026-11-11.md` has Parts 6–7 (the second- and
  third-order terms), `journal/2026-11-10.md` has Parts 1–5,
  `journal/2026-11-09b.md` is the first part of session 3 and `journal/2026-11-09.md` is session 2.

## Where the next steps are

`docs/03-research-program.md` and the newest entry in `journal/`.

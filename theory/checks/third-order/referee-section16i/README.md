# Referee report (ninth pass on Section 16): Lemma 16.13, Theorem 16.14, Remark 16.15

Independent referee checks, made on 14 November 2026, of the items added at the end of Section 16 of `theory/mass-gap-strong-coupling.tex`
(working tree of commit `ceebc44`, branch `claude/nice-lamport-ciczuy`; I ran no git command and changed no file outside this folder):

* Lemma `mg:lem:upsbounds` (Lemma 16.13, lines 4993–5010);
* Theorem `mg:thm:thirdcap` (Theorem 16.14, lines 5012–5075), the computer-assisted theorem;
* Remark `mg:rem:thirdcap` (Remark 16.15, lines 5077–5090);
* what depends on them: the abstract sentence (lines 141–143), the date block (lines 64–65), the end of the Section 16 introduction (lines 4222–4223), the
  paragraph before the lemma (lines 4987–4991), the summary-table rows (lines 5244, 5247, 5249–5251) and the new rows 50–53 of `../README.md`.

Nothing here bears on the Yang–Mills problem itself, which is open; the items refereed are statements about lattice Yang–Mills theory at strong
coupling ($\beta<\beta_0(d)\approx4\cdot10^{-4}$).

## Recommendation

**Minor revisions.** The result stands. I re-derived every hypothesis of Theorem 14.2 (`mg:thm:lower5`) in bound mode from the formulas of the paper,
in my own code, and evaluated it at $\beta_*=\beta_0(d)(1+10^{-30})$ in **Arb ball arithmetic** (python-flint, 400 bits) with **rigorous remainders for
every infinite series**: all hypotheses hold in all four cases ($SU(3)$, $SU(5)$; $d=3,4$), and $(1+\kappa)^3D(\beta_*)/\beta_*^4=2.557639$, $3.936887$,
$0.345985$, $0.519411$, below $c_-=2.558$, $3.937$, $0.346$, $0.520$. The upper constants $c_+$ are confirmed in Arb as well, and Lemma 16.13(ii) in exact
rational arithmetic by an independent method. All four bound-mode scalings, `error_budget.py` and `certify_upper.py` reproduce to the last printed digit.

There is **one error**, in the verification bookkeeping, not in the numbers: the marked constants that the code uses violate inequalities (5) and (6)
of Definition `mg:def:jconst` by a relative $\approx5\cdot10^{-16}$, and the paper's argument that they are safe (Theorem 16.14(c)) is wrong (Issue 1).
Two repairs are available, each of which I checked, and neither changes any constant. There is one **gap of presentation** that matters for the label
"Theorem": the proof rests on arguments that sit in Remarks, which the paper's Convention calls commentary (Issue 2). The rest are a formal gap in the
coverage of `error_budget.py` (numerically immaterial), the rounding argument (now settled by my ball-arithmetic run), and imprecisions of wording,
ranges and the summary table.

## Issues

Line numbers refer to `theory/mass-gap-strong-coupling.tex`.

1. **Error (verification of the marked constants; no effect on the constants).** Theorem 16.14(c), line 5054: "the marked constants solve monotone
   $2\times2$ systems with determinant $\ge0.57$, so they move by much less". `third_order_constants.marked` (lines 45–49 of `third_order_constants.py`) solves the $\kappa$-system with
   the *un-raised* $g$'s and then multiplies all six constants by `INFLATE` $=1+\iota$, $\iota=10^{-15}$. The bilinear sources $g_Fg'_F$, $g_Cg'_C$ of
   inequalities (5), (6) (lines 1703–1704) then grow by $(1+\iota)^2$ while $\kappa_F,\kappa_C$ grow by $1+\iota$, so (5) and (6) fail by
   $\iota(1+\iota)g_Fg'_F$, $\iota(1+\iota)g_Cg'_C$. Measured (`my_eval.py`, `outputs/my_eval.out`), as $1-\text{lhs}/\text{rhs}$:
   * actual constants, all four cases: (5) $-1.08\cdot10^{-16}$, (6) $-4.93\cdot10^{-16}$;
   * barred constants: (5) $-0.55$ to $-0.81\cdot10^{-16}$, (6) $-2.99$ to $-4.07\cdot10^{-16}$.

   The values are the same with the authors' own truncated sums, so this is not a tail effect; and the determinant argument does not apply, because the
   $\kappa$-system is not linear in the raised $g$'s. The same holds at every $t$ of `scaling_check.py`, which calls the same routine. Either repair suffices:
   * (i) Lemma `mg:lem:marked` (line 1672) needs only $G_1G_2$, i.e. the left sides of (1)·(3) and (2)·(4), not $g_Fg'_F$, $g_Cg'_C$. In that form the authors'
     constants satisfy (5), (6) with room $\ge1.0\cdot10^{-16}$ (actual) and $\ge5.5\cdot10^{-19}$ (barred, $d=4$; $\ge3.2\cdot10^{-17}$ in $d=3$), with the
     extended sums included. Smallest: $5.55\cdot10^{-19}$, $SU(3)$, $d=4$, barred (5).
   * (ii) Solve the $\kappa$'s with the raised $g$'s, then raise them. This is my construction: all six inequalities hold, with room $\ge2.8\cdot10^{-31}$ at
     $\iota=10^{-30}$. Every hypothesis and $D$ are then unchanged to $\ge8$ digits (`my_eval.py`), and certified in Arb (`arb_eval.py`).

   The sentence in (c), and either the code or Definition `mg:def:jconst`/Lemma `mg:lem:jsums` (state the weaker form), should be corrected. The six
   inequalities should also be checked explicitly, as the decoration inequalities already are.

2. **Gap (presentation; matters for the label "Theorem").** By the Convention (line 191, "Remarks are commentary"), a Theorem's proof may not rest on
   Remarks. The proof of Theorem 16.14 rests on four:
   * the scaling argument of Remark 16.7(vii), lines 4690–4724, which proof (b) only adapts;
   * the bound $E(\beta)$ and $\Upsilon\ge1$ for $SU(3)$, Remark 16.9(ii), lines 4851–4862, used in (e);
   * $Z_{n,4}=0$ for $n\ge3$ in $d=3$, Remark 16.7(iv), lines 4613–4628, which `load_terms` uses for all $n$ in $d=3$;
   * the bookkeeping that the sets $S''_9$ and $S^\times_{10}\cup\{(9,2)\}$ have the stated corners, Remarks 16.7(vii), (x), 16.12(iii).

   The mathematics in them looks right to me (§D, §H below), but they should be promoted to lemmas or written into the proof. The scaling lemma should be
   stated in bound mode, with its list of monotonicity facts; the list in §D below could serve.

3. **Gap (formal; numerically immaterial).** Theorem 16.14(c) (lines 5049–5053) and row 52 of `../README.md` ("every truncated series") overstate what
   `error_budget.py` does. It does not bound the tails of:
   * the full sum $\Sigma_4$ (both $n$-sums; used in $\eta$ and in the third hypothesis);
   * the first $n$-sum of $\Sigma^{(6)}$ (the four-edge fronts);
   * the one-plaquette sums $\mathsf A,\mathsf A_8,\mathsf F,\mathsf F_8$;
   * $Y_1,Y_2$ at the actual $b$, and $A_1,A_2$ at the actual $a$;
   * the series for $e^{\mathsf A}-1-\mathsf A$, cut at 59 terms.

   I computed all of them (`tails_check.py`, `outputs/tails_check.out`). The largest is $4.6\cdot10^{-52}$ relative (the crude $Z_{n,4}$ in $\Sigma_4$, $d=4$),
   so nothing changes, but the claim should match the script, or the script should list them. `arb_eval.py` bounds every remainder rigorously.

4. **Gap, closed here (rounding).** Remark 16.15, lines 5081–5083. The argument that mpmath's rounding is negligible is informal: only $+,-,\times,/$
   are correctly rounded; $\exp$, $\log$ and the binomials for large arguments (via $\Gamma$) are accurate but carry no stated error bound; and no
   propagation bound is given.
   * Measured: the 50- and 100-digit evaluations of every quantity that enters a hypothesis agree to $\sim10^{-50}$ relative (`precision_check.py`),
     against rooms $\ge6.7\cdot10^{-16}$.
   * Proved: `arb_eval.py` redoes the evaluation at $\beta_*$ in ball arithmetic with rigorous remainders and certifies every hypothesis
     (`outputs/arb_eval.out`). `upper_arb.py` certifies $E(\beta)/\beta^4<c_+$ at an exact point above $\beta_0(d)$.

   I suggest adopting such a ball-arithmetic evaluation (python-flint is installed in the container) and citing it in Remark 16.15.

5. **Imprecision.** Lemma 16.13(ii), line 4996, states $0<\beta\le4.03\cdot10^{-4}$. The proof (line 5009) and `rigorous_bounds.py` (which calls
   `prove_upsilon_ge_1(5, beta0(3))`) cover only $\beta\le\beta_0(3)(1+10^{-20})=4.0262\cdot10^{-4}$. The statement is true: the authors' routine
   passes with `beta_max = 4.03e-4` (margin $0.0208333$), and so does my exact-rational check. Align the text. Also, line 5008, "Bound the rest of $P$
   by these majorants": the degrees $4,\dots,30$ are bounded by their exact $|p_k|$, which gives the $5.6\cdot10^{-11}$; only the degrees $>30$ use the
   majorants.

6. **Imprecision.** Theorem 16.14(a), line 5032: "each hypothesis is monotone in each of them" is not literally true, though what the code does is
   right.
   * $u$ enters $\Sigma^{(4)}_4$ increasingly (through $a=\lambda u$) and decreasingly ($b_2\propto u^{-4}$, the shifts and jumps).
   * $\bar\Upsilon$ enters $\eta$ and the third hypothesis decreasingly (through $C_\varpi$) and increasingly (through $2\bar\Upsilon(1+\delta')$).

   What is used is that every *term* is monotone. It also uses that Theorem 14.2 holds with $C_\varpi$ replaced by any larger number, and with
   $\bar\Upsilon$ replaced by any $\bar\Upsilon'\ge\norm{P_-P^qD_qP^qP_-}/u^4$; here $\bar\Upsilon'=\Upsilon+\Delta\in[e^{-2\beta},e^\beta+\Delta]$. This is
   legitimate, since the proof uses them only through $\nu(\pi)\le C_\varpi w_\varpi(\pi)$ and $\nu^-(D_q)\le u^4\bar\Upsilon$. It should be said, because
   it is what justifies using different values of $\bar\Upsilon$ in different hypotheses.

7. **Imprecision.** Remark 16.15, line 5084: "no step uses the exact group data: only $N$, $\iota$ and the proven bounds of (a)". This is not accurate.
   * The $SU(5)$ upper bound uses Lemma 16.13(ii), an exact computation with the Haar moments of $SU(5)$.
   * $u_*=u$ (Corollary `mg:cor:window`) is a group fact too.
   * `certify_upper.py` calls the exact group data, although only for an assertion.

   "No step uses floating-point values of the group data" would be accurate.

8. **Imprecision.** Lines 5049 and 5081, "It runs at 50 digits":
   * the evaluation at $\beta_*$ in `scaling_check.py` runs at 60 digits (the module default of `third_order_effective.py`); `error_budget.py` runs at 50;
   * `error_budget.py` and `certify_upper.py` evaluate at $\beta_0(d)$, not at $\beta_*$.

   Both are immaterial: my runs at $\beta_*$ agree to every printed digit.

9. **Imprecision.** Line 5069, "the exits with four and five plaquettes and exits of at least $17$ edges". The code (`practical_beta.load_terms`, lines
   76–84) makes crude only the exits with $n\in\{4,5\}$ *and* $k\ge17$. Write "with four or five plaquettes and at least $17$ edges". The list in (d)
   should also name the $d=3$ exits that are used: $n=2,3$, and $n=4,5$ with $k\le16$, all made by two programs. Effect on $D$: none to 12 digits
   (`certified_effect.py`).

10. **Imprecision.** Summary table, lines 5244 and 5249: "(proved in Thm. thirdcap)" is attached to constants the theorem does not prove.
    * $0.033\beta^4$ ($SU(3)$, $d=4$, upper): the theorem gives $0.047$.
    * $-0.519\beta^4$ and $0.252\beta^4$ ($SU(5)$, $d=4$): the theorem gives $0.520$ and $0.361$.

    Write "(proved, with the constants of Thm. 16.14)". The abstract (lines 141–143) says "these statements become a computer-assisted theorem"; adding "with
    slightly weaker constants in $d=4$" would be exact.

11. **Suggestion (checks README).** The header of `../README.md` still says "A compass, not part of the proofs", but `rigorous_bounds.py`, the
    `--bounds --certified` mode, `error_budget.py` and `certify_upper.py` are now part of a computer-assisted proof. Say so. Also, the four scaling commands
    (with the $\beta_*$ strings) are written down nowhere, and the `cap_*_scaling.out` files carry no command line. I had to reconstruct them; they are in
    `run_reproductions.sh`, and with them the outputs are identical.

12. **Suggestion (code).** `third_order_effective.py` line 203 computes $\log(1+\delta')$ as `log(1+dl)`. With $\delta'\sim10^{-14}$ this loses 14 digits
    ($D$ at 50 digits agrees with $D$ at 100 digits only to $10^{-38}$). Use `log1p`. Harmless at the 60 and 130 digits used.

13. **Suggestion (labels).** The Standards of the README list Theorem, Cited, Conjecture, Heuristic and Numerical. In the Status column of the summary table,
    "Theorem (computer-assisted)" for Lemma 16.13 and Theorem 16.14 would make the dependence on enumerations and machine evaluation visible where readers
    look for the status.

## What I checked, and how

### A. Lemma 16.13(i) (item 1 of the brief)

* **The definition.** Section 10, display before Lemma `mg:lem:Upsilon` (line 2193): $\Upsilon$ is the maximum of the variance under $w_\beta/c_0$ over
  $\psi\in\mathrm{span}_\R\{\mathrm{Re}\,\chi_\rho,\mathrm{Im}\,\chi_\rho\}$ with $\int\psi^2\dU=1$. The real elements of the complex span of
  $\chi_\rho,\overline{\chi_\rho}$ are exactly this real span: for $\iota=0$ write $\psi=a\chi+\bar a\bar\chi$; for $\iota=1$, $\psi=c\chi$ with $c$ real. So
  the lemma's reformulation is the paper's definition, and $\int\psi\dU=0$ since $\rho$ is non-trivial.
* **The weight.** $w_\beta=\exp(\frac\beta N\mathrm{Re}\,\chi_\rho)=e^{t(\chi_\rho+\overline{\chi_\rho})}$, $t=\beta/(2N)$ (as in the proof of Lemma
  `mg:lem:u2sun`), so $w_\beta\in[e^{-\beta},e^\beta]$ since $|\mathrm{Re}\,\chi_\rho|\le N$. And $c_0=\int w_\beta\ge e^{\int(\beta/N)\mathrm{Re}\chi}=1$ (Jensen), $c_0\le e^\beta$.
* **The two bounds.**
  * Upper: $\Var_\mu\psi\le\int\psi^2w_\beta/c_0\le e^\beta$.
  * Lower: $\Var_\mu\psi=\min_c\int(\psi-c)^2d\mu\ge e^{-\beta}e^{-\beta}\min_c\int(\psi-c)^2\dU=e^{-2\beta}$.

  Both are correct line by line.
* **Numerical sanity check** (`ups_bounds_check.py`, by Weyl integration). $c_0\in[1,e^\beta]$ and $e^{-2\beta}\le\Upsilon\le e^\beta$ hold for $U(1)$ and
  $SU(2)$ at $\beta=0.01,0.3,1,3,10$, and for $SU(3)$ at $\beta=0.01,0.3,1,3$; e.g. $SU(3)$: $\Upsilon=1.0016667$ at $0.01$ (that is $1+\beta/6$).

### B. Lemma 16.13(ii) (item 2)

* **Independent series** (`ups_ge1_check.py`). I recomputed the series of $Z(t)=\int_{SU(5)}e^{t\,\mathrm{Re}\,\chi}\dU$ from the Bessel-determinant formula
  $Z=\sum_{k\in\Z}\det[I_{k+i-j}(t)]_{i,j\le5}$, not from Schur–Weyl moments, in exact rational arithmetic. Then $Z_1=Z'$, $Z_2=Z''$ and
  $P=2Z_2Z-2Z_1^2-Z^2$.
  * Haar moments $m_2,\dots,m_{10}=\frac12,0,\frac34,\frac1{16},\frac{15}8,\frac{35}{64},\frac{105}{16},\frac{261}{64},\frac{7581}{256}$.
  * $p_0=p_1=p_2=0$, $p_3=\frac1{48}$, $p_4=0$, $p_5=\frac{11}{1280}$, $p_6=0$, $p_7=\frac{47}{26880}$, $p_8=\frac1{737280}$.
  * Equal to the authors' $P$ for every degree $\le30$.
* **The majorants.** The coefficient of $t^k$ in $Z$, $Z_1$, $Z_2$ is $\frac1{k!}\int(\mathrm{Re}\,\chi)^{k+j}\in[0,5^{k+j}/k!]$ for $j=0,1,2$, since the
  moments of $\chi,\bar\chi$ are multiplicities, hence $\ge0$, and $|\mathrm{Re}\,\chi|\le5$. So the majorants $5^k/k!$, $5^{k+1}/k!$, $5^{k+2}/k!$ are right.
* **The tail.** `prove_upsilon_ge_1` bounds the part of degree $>30$ of the three products by the degree-$>30$ part of the products of exact truncations
  plus majorant tails, all with non-negative coefficients. So it is free of cancellation, and divided by $t^3$ it is a series with non-negative
  coefficients, hence increasing. The degrees $4..30$ enter as $\sum|p_k|t^{k-3}$, also increasing. One evaluation at $t_{\max}$ therefore suffices.
  Its sum to $K+59$ times $(1+10^{-20})$ covers the rest of the exponential tail by a factor $>10^{100}$.
* **My bound.** $|p_n|\le(4N^2+1)(2N)^n/n!$, evaluated exactly in rationals, gives margin $0.0208333332775>0$ on $(0,4.03\cdot10^{-4}]$; the authors'
  routine gives the same with `beta_max=4.03e-4`. $U_+-1=\beta^3/6000+O(\beta^5)$ as stated. And $\Upsilon\ge U_+$ by Lemma `mg:lem:Upsilon`(i). (Range: Issue 5.)

### C. Bound mode, place by place (item 3, the main one)

I read `lower3e` with `g['bounds']` against Theorems `mg:thm:lower2`, `mg:thm:lower4`, `mg:thm:lower5`, Definition `mg:def:jconst`, Lemmas
`mg:lem:jsums`, `mg:lem:jsums4`, `mg:lem:jsums6`, `mg:lem:corners`, `mg:lem:dressed`, `mg:lem:theta7`. Every group-dependent quantity enters as follows
(✓ = safe direction; "line" = line of `third_order_effective.py`):

| Where | Quantity used | True quantity | Direction |
|---|---|---|---|
| $\theta$ (line 60) | $2(d-2)e^6\bar\sigma_6+\dots$ | $\theta_7$ with $\sigma_6$ | any $\theta\ge\theta_7$ is allowed (proof of Lemma `mg:lem:theta7`); ✓ |
| $a=\lambda u$ (line 64) | $\lambda\bar u$ | $\lambda u$ | ✓ (see below) |
| $b=\lambda\gamma$ | $\lambda\gamma_2$, group-free | $\lambda\gamma_g$, $\gamma_g\le\gamma_2$ | ✓ |
| $C_\varpi$ (line 70) | $1+\max(1/\varpi,(1+\gamma_2)/e^{-2\beta})$ | $1+\max(1/\varpi,(1+\gamma_g)/\bar\Upsilon')$ | larger $C$ allowed (Issue 6); ✓ |
| $\omega$, $c_J$, $h_J$-condition | with $\bar C$, $\bar a$ | — | increasing in $C$, $a$; ✓ |
| $b_2$ (line 71) | $\lambda^2(2\bar\sigma_6+\bar\sigma_6^2)/u_{\rm low}^4$ | $\lambda^2(2\sigma_6+\sigma_6^2)/u^4$ | ✓ (the only place where $1/u$ enters) |
| $A_j$, $\mathsf A$, $f(k)$, exits, one-plaquette, $\Sigma^{(6)}$, $\Sigma_4$ | $\bar a$ | $a$ | all terms increasing in $a$; ✓ |
| $f(k)/a^{k-4}$ in the exact corner terms | $a$ cancels | — | ✓ |
| `chars` (Theorem 14.2) | $10\bar a^4$, $K_S\bar af(6)$, $2\bar a$ | $u_*=u$ exactly (Lemma `mg:lem:moves`, Cor. `mg:cor:window`) | ✓ |
| majorant $(b/\bar b)^n(a/\bar a)^E$, grid `ab >= a` | $\bar a$ | $a$ | ✓ |
| right sides of the $\delta'$-, $F$-, $F_6$-conditions | $e^{-2\beta}(1+\delta')-\bar a$, $e^{-2\beta}\delta'$ | $\bar\Upsilon'(\cdot)-a$ | ✓ |
| $\eta$, third hypothesis, $a^4\bar\Upsilon(1+\delta')\le\varpi$, $\le1$ | $e^\beta+\Delta$ | $\bar\Upsilon'\le e^\beta+\Delta$ | ✓ |
| deficit | $8(K+1)\bar\theta+e^{2\beta}\Delta+\log(1+\delta')$ | $8(K+1)\theta-\log\bar\Upsilon'+\log\Upsilon+\log(1+\delta')$ | $\log(1+\Delta/\Upsilon)\le e^{2\beta}\Delta$; ✓ |

* **No unbounded group data.** `g['u']` (that is $\bar u$) is used only in $a$, which appears only in increasing positions on the left sides and as
  $-a$ on the right sides. `Ups` is never used in bound mode. `dressed()` is not called.
* **$\Delta$.** It is right: $\bar\Upsilon-\Upsilon=n_B^2/(\Upsilon-c_2)+2e^{2\beta}\sigma_6/(\upsilon_1-u^4(1+\gamma_g))$, and $\Upsilon,\upsilon_1\ge e^{-2\beta}$,
  $n_B\le\bar n_B$, $c_2\le\bar c_2$, $\sigma_6\le\bar\sigma_6$. The denominators are asserted positive.
* **My own implementation.** `my_eval.py` implements all of this from the paper's formulas (not from `lower3e`), with sums to 260 edges and 130
  plaquettes. With the authors' parameters at $\beta_*$:
  * $D/\beta^4$ agrees with theirs to 10 digits in all four cases ($2.554552302$, $0.345567128$, $3.93507159$, $0.5191715676$);
  * every part of $\Sigma^{(4)}_4$ agrees to 8 digits, and $\Sigma^{(6)}/b$ to 12;
  * rooms: decoration $\ge9.7\cdot10^{-16}$, barred $\ge6.70\cdot10^{-16}$, $\delta'$-, $F$-, $F_6$-conditions $\ge1.0\cdot10^{-9}$, the others $\approx1$.

### D. Monotonicity and the scaling (item 4)

The scaling of Remark 16.7(vii) uses the following facts. In bound mode each holds by inspection, as checked here (and confirmed on 2000 points of
$[10^{-24},\beta_*]$ by `monotone_bounds_check.py`, which is only a sanity check):

1. $\bar a/\beta=\lambda\,((1+\iota)+\beta e^\beta)/(2N^2)$ increases, since $\bar\theta$, hence $\lambda=e^{2(K+1)\bar\theta}$, increases.
2. $b/\beta=\lambda(\epsilon/\beta+\beta e^\beta)$ and $\epsilon/\beta$ increase.
3. $\bar u_2/\beta^2$ increases: a maximum of $\beta^{-2}\cdot$(order-2 term $+R/d_r$)$/d_r$ and of $R/\beta^2$, with $R/\beta^2=\beta\sum\beta^{j}/(j+3)!$.
4. $\bar\sigma_6/\beta^6=(2-\iota)N^2(\bar u/\beta)^6+(\bar u_2/\beta^2)^5\beta^4(e^\beta-1)$ increases.
5. $\bar\theta/\beta^6$ increases: every term is $\beta^n$ times a series with non-negative coefficients, $n\ge6$.
6. $\bar b_2/\beta^2=\lambda^2(2\bar\sigma_6+\bar\sigma_6^2)\beta^{-6}\cdot(\beta/u_{\rm low})^4$ increases, with $\beta/u_{\rm low}=2N^2e^\beta/(1+\iota)$.
7. $e^{2\beta}\Delta/\beta^4=\beta^2e^{2\beta}\Delta/\beta^6$ increases (it replaces $\log(\bar\Upsilon/\Upsilon)/\beta^4$):
   * $\bar n_B/\beta^3$ and $\bar c_2/\beta^4$ increase, as $1/u_{\rm low}$ enters only as $e^\beta/\beta$;
   * the denominators $e^{-2\beta}-\bar c_2$ and $e^{-2\beta}-\bar u^4(1+\gamma_2)$ decrease and stay positive.
8. "$\bar\Upsilon-a\ge1$, increasing" is replaced: $e^{-2\beta}-\bar a$ *decreases* in $\beta$. So at $\beta=t\beta_*$ the right sides of the $F$- and
   $F_6$-conditions are $\ge(1+\kappa)t\,F_*(e^{-2\beta_*}-\bar a_*)$ (resp. with $F_{6*}$).
   * The left sides are $\le t\cdot(\text{value at }\beta_*)\cdot(1+O(\kappa\beta^3))$: the $f(k)$, $k\ge6$, carry $1+\kappa$, which reaches $g_F,g_C$ only
     through $B_0-N^v_4e^{4\xi_F}$.
   * So $1+\kappa$ has to absorb only $e^{-2\beta_*}\delta'_*$ (plus that $O(\kappa\beta^3)$), as the theorem says. Here $\delta'_*\approx10^{-14}$ against
     $\kappa=\beta_*\approx10^{-4}$.
9. "$\bar\Upsilon(\beta_*)-1\le\kappa/5$" is no longer needed: the right side $e^{-2\beta}\delta'$ of the $\delta'$-condition is
   $\ge e^{-2\beta_*}(1+\kappa)^2t^4\delta'_*$. Of the $(1+\kappa)^2$, one factor absorbs the $1+\kappa$ that the $f(k)$ carry into $\Sigma^{(4)}_4$:
   * linearly into the exits, the one-plaquette term and the exact corner terms;
   * linearly into the barred marked constants, since $\kappa$ is bilinear only in $(g,g')$ and $g'$ does not depend on $f$.
10. $\bar C_\varpi=1+\max(t^{-p}/\varpi_*,(1+\gamma_2)e^{2\beta})\le t^{-p}\bar C_{\varpi*}$. Hence:
    * $\omega\le t^{4-p}\omega_*$, so $c_J\le c_{J*}$;
    * the $h_J$-condition scales as $t^{1+k/2-p}\le t^{1/2}$ ($p=\frac52$ in all four cases);
    * $\bar a^4\bar\Upsilon_{\rm hi}(1+\delta')/\varpi=O(t^{3/2})$, as $\bar\Upsilon_{\rm hi}=e^\beta+\Delta$ increases.
11. The terms of $\Sigma^{(4)}_4$ have the orders $t^{p+2}$ (shifts and jumps), $t^4$ ($10\bar a^4$, $\bar af(6)$), $t^5$ (one plaquette), and $t^{n+E}$
    with $n+E\ge5$ (exact pairs, one more power of $t$ when $|p'|>4$; the corners have $n+E\ge9$). $\eta\le(1+\kappa)t^{4-p}\eta_*$.
12. $D(\beta)/\beta^4\le(1+\kappa)^3D_*/\beta_*^4$: $\bar\theta$ and $e^{2\beta}\Delta$ scale at least as $t^6$, and
    $\log(1+\delta')\le(1+\kappa)^2t^4\delta'_*\le(1+\kappa)^3t^4\log(1+\delta'_*)$, as $\delta'_*\le\kappa$.

So (b) is right. With my evaluator and correctly built marked constants, the scaled parameters are admissible at $t=1,0.9,0.5,0.1,10^{-3},10^{-6},10^{-12},
10^{-20}$ in all four cases, with $D/\beta^4\le(1+\kappa)^3D_*/\beta_*^4$ (`scaled_check.py`, `outputs/scaled_check.out`).

### E. Error budget (item 5)

* **The majorants are correct.**
  * Clusters: $\binom{K_Sm}{m-1}\le(eK_Sm/(m-1))^{m-1}\le e(eK_S)^{m-1}$, with ratio bound $(1+1/M)^peK_Sy$.
  * Fronts: Lemma 16.10. I rechecked the constants: $(1+s)^6-1-6s=0.733517$ and $y/s=7.35428$ at $s=0.194$; $0.782219$ and $4.08063$ at $s=0.325$ in $\Z^2$.
  * Polymers: $\binom{c+Tn}n\le(e(T+c/n))^n$.
  * Exits: $\sum_{k\ge5}\binom Lkk^2z^k\le L^2\binom L5z^5(1+z)^{L-5}$, by $\binom L{5+j}\le\binom L5\binom{L-5}j$ and $k\le L$. The ratio
    $\binom{L+4}5/\binom L5\le((L+4)/(L-1))^5$ holds for $L\ge16$, factor by factor.
  * The $k^4$-weighted relative tail dominates those of $B_0,B_1$, since the weights increase in $k$.
* **Reproduced:** largest relative tail $1.40\cdot10^{-21}$ ($A_2$ at $\bar a$, $SU(3)$, $d=4$; my extended sums give $1.39\cdot10^{-21}$); rooms with the
  tails $\ge6.70\cdot10^{-16}$.
* **Missing families**, all $\le4.6\cdot10^{-52}$ relative: Issue 3.
* **The marked-constant claim is wrong:** Issue 1.
* **"Room $\ge10^{-12}$" for the other hypotheses** is right: they are $\ge1.0\cdot10^{-9}$, because $\delta'$ gets a final extra factor $1+10^{-9}$.
* **mpmath vs interval arithmetic:** Issue 4. The Arb evaluation settles it for the four cases at hand.

### F. Provenance of the counts (item 6)

Every count that `--certified` keeps has a second, independent program, and the tables in `practical_beta.py` equal the referees' published numbers:

* $N^v_k$:
  * $\Z^2$ to 12: `referee-section16c/fronts_redelmeier.c`, `referee-section16f/fronts_enum.c` ($1728$, $6369$);
  * $\Z^3$ to 13: `fronts_enum.c`, `referee-section16g/fronts_redel.c`.
* $Z_{n,4}$, $d=4$, to 12: `referee-section16f/junc_partition.c` ($3193$, $5911$, $26956$). In $d=3$, $Z_{n,4}=0$ for all $n\ge3$ by the planar argument of
  Remark 16.7(iv) (Issue 2).
* Exits:
  * $Z_{2,k}$: `referee-section16c/zexit_check.py`;
  * $d=4$: $Z_{3,k}$ and $Z_{4,k}$ for $k=6..20$ (the authors' table, identical): `referee-section16f/zexit_pfirst.py`, rerun in `referee-section16g/`;
  * $d=3$: $Z_{3,k}$, and $Z_{4,k},Z_{5,k}$ for $k\le16$: `referee-section16d/junc_faces.c` ("$k\ge17$ not checked", hence the crude counts there).

  Same-layer exits do not arise in $d=3$ for $n\le5$: by the face argument of `referee-section16d`, $P$ would contain the annulus between $\partial q$ and an
  enclosing vertex-disjoint cycle, at least 8 plaquettes. So the "other layer only" counts there are complete.
* Junction classes:
  * $S''_9$ ($d=3$): `junc_faces.c`, all 44 terms;
  * $S^\times_{10}$: the seventh referee's partition search;
  * the $(9,2)$ class, $148894$: `junctions.c` and `referee-section16g/junc_partition_big.c`, rerun by the eighth referee.
* Closed polymers:
  * $\Z^3$ to 14: `ref_closed.py`, `closed_search.py`;
  * $\Z^4$ to 16: `referee-section16g/closed_partition.c`.
* $N_3,N_2,N'_2$: `referee-section16b/`.
* Admissible sets: $\Z^3$ 11–12 and $\Z^4$ 9–10 by the fifth referee's two methods. These are "two programs" but not "ours and a referee's"; (e) words it
  correctly.

### G. Reproductions (item 7)

All with `run_reproductions.sh` (2 jobs). Outputs in `outputs/repro_*.out`. The scaling outputs, the error budget and the upper-bound output are identical to the files in `../orders/`,
apart from timing lines and header comments.

| Run | Result | Time |
|---|---|---|
| `scaling_check.py --bounds --certified 0.000402621433522819744012853968950678326... 2` | $D_*/\beta_*^4=2.55455$; $(1+\kappa)^3D_*/\beta_*^4=2.55764<5.7494$; 17/17 $t$ pass | 38 s |
| the same with `--group 5` | $0.345567$; $0.345985<0.41395$; 17/17 | 41 s |
| `--d4 orders/d4_s10x_92.json --closed 10:60,…,16:12552 --closed-from 17 0.000153766448042548273791723769079994… 2` | $3.93507$; $3.93689<15.054$; 17/17 | 40 s |
| the same with `--group 5` | $0.519172$; $0.519411<1.0839$; 17/17 | 41 s |
| `certify_upper.py` | $E(\beta_0)/\beta_0^4=0.041462$, $0.0469392$, $0.317825$, $0.360293$; $\Upsilon\ge U_+\ge1$ for $SU(5)$ (margin $0.0208333$) | 2.6 min |
| `error_budget.py` | $\tau\le1.40\cdot10^{-21}$, room $\ge6.70\cdot10^{-16}$, $\det\ge0.5696$ | 1.3 min |
| `rigorous_bounds.py` | as `orders/`; with `beta_max=4.03e-4` also True | 0.3 s |

My independent checks:

| Script | What | Result |
|---|---|---|
| `ups_ge1_check.py` | Lemma 16.13(ii) by Bessel determinants, exact rationals | §B; margin $>0$ on $(0,4.03\cdot10^{-4}]$ |
| `ups_bounds_check.py` | Lemma 16.13(i) by Weyl integration, $U(1)$, $SU(2)$, $SU(3)$ | all within $[e^{-2\beta},e^\beta]$ |
| `my_eval.py G d` | bound-mode evaluator from the paper; checks of the authors' constants | §C; Issue 1 |
| `tails_check.py G d` | tail of every truncated series | Issue 3 |
| `scaled_check.py G d` | scaling with correct marked constants | §D |
| `arb_eval.py G d` | ball arithmetic at $\beta_*$, rigorous remainders, own marked constants, $\beta_*\ge\beta_0(d)$ checked | every hypothesis certified, 4/4; $(1+\kappa)^3D/\beta_*^4\le2.557639$, $0.345985$, $3.936887$, $0.519411$ |
| `upper_check.py`, `upper_arb.py` | $E(\beta)$ from the text of Remark 16.9(ii), Lemma 16.8; Arb at a point above $\beta_0(d)$ | $0.041462$, $0.046939$, $0.317825$, $0.360293$, all $<c_+$; $\binom{268}{11}=1.042\cdot10^{19}$, $\binom{240}{11}=3.02\cdot10^{18}$ |
| `monotone_bounds_check.py` | the bound-mode monotonicity facts on 2000 points | all increasing |
| `certified_effect.py` | `--certified` vs all exits, $d=3$ | $D$ unchanged to 12 digits |
| `precision_check.py G d` | authors' `lower3e` at 50 vs 100 digits | parts agree to $\sim10^{-50}$; $D$ to $10^{-38}$ (Issue 12) |

Commands, from `theory/checks/third-order`:

```
bash referee-section16i/run_reproductions.sh 2
cd referee-section16i
python3 ups_ge1_check.py 5 34; python3 ups_bounds_check.py
for c in "3 3" "5 3" "3 4" "5 4"; do python3 my_eval.py $c; python3 tails_check.py $c; python3 scaled_check.py $c; python3 arb_eval.py $c; done
python3 upper_check.py; python3 upper_arb.py; python3 monotone_bounds_check.py; python3 certified_effect.py
python3 precision_check.py 3 3; python3 precision_check.py 3 4
```

Total CPU time of all runs: about half an hour, on at most 3 cores.

### H. The theorem and remark text (item 8)

* **$\beta_*$.** $\beta_*=\beta_0(d)(1+10^{-30})\ge\beta_0(d)$ is right, and needed only as a starting point of the scaling: the hypotheses are numerical
  inequalities there, and the theorem is applied at $\beta=t\beta_*<\beta_0(d)$. `arb_eval.py` checks $\beta_*>\beta_0(d)$ in balls.
* **The constants are rounded up correctly.**
  * $c_-$: $2.557639\to2.558$, $3.936887\to3.937$, $0.345985\to0.346$, $0.519411\to0.520$.
  * $c_+$: $0.041462\to0.042$, $0.046939\to0.047$, $0.317825\to0.318$, $0.360293\to0.361$.
* **The comparisons.** For $\beta<\beta_0(d)$, $\beta^3/432>\beta^4/(432\beta_0)$: $5.7494$ ($d=3$) and $15.054$ ($d=4$). Likewise $\beta^3/6000$: $0.41395$ and
  $1.0839$. So "$\ge5.74$, $15.05$, $0.413$, $1.083$" are right, and all $c_\pm$ lie below them. The closest margin is $SU(5)$, $d=3$:
  $c_-=0.346$ against $0.413$.
* **Interpretation.** "Determined from both sides" means, as in Remarks 16.7, 16.9 and 16.12: the residual $|m-4\log(1/u)+\log\Upsilon|$ is smaller than the
  term at every $\beta$ of the range.
* **The upper bound (e).** It is correct.
  * Lemma 16.8 holds with any $N_n\ge$ the count, and $\binom{J'+Kn}n$ is the rooted count of the proof of Lemma `mg:lem:twopoint`, with $J'=48$, $20$ and $K=20$.
  * $E/\beta^5$ is a series with non-negative coefficients, so $E(\beta)\le(\beta/\beta_0)^5E(\beta_0)\le E(\beta_0)\beta^4/\beta_0^4$.
  * $\Upsilon\ge1$ is used three times in the derivation of $E$ ($\log(V_2/\Upsilon)$, $(\Upsilon-x)^2\ge\Upsilon^2(1-2x)$, and the $\mathcal T_Q$ term). I rederived the
    derivation, $V_2\le(1+\mathcal T)\Upsilon+\mathcal T(1+3a+3a^2)$ and $L_2\ge\Upsilon^2(1-B)$, and it is right.
  * The $SU(3)$ argument $\kappa-18u^2>0$ is right.
* **"Computer-assisted".** The description is honest about the enumerations and the evaluation. See Issue 4 (rounding), Issue 7 (exact group data) and
  Issue 2 (Remarks).

## What I could not check

* **The enumerations.** I did not recount any of them: no new enumeration code. I checked only their provenance, and that the tables used equal the
  numbers published by the referees.
* **The earlier lemmas.** The correctness of Theorems 14.2/13.4 and of the lemmas they use was refereed earlier. I checked only how the bound-mode data enter
  their hypotheses, and that my reading of the formulas reproduces the authors' sums to 8–12 digits.
* **The C search programs.** I did not read them.
* **The upper bound beyond its formula.** I did not re-derive Theorem 10.2 (`mg:thm:upper2`) and Lemma 16.8 beyond checking the formula for $E(\beta)$ and
  its use of $\Upsilon\ge1$.

## Files

`run_reproductions.sh`; `ups_ge1_check.py`; `ups_bounds_check.py`; `my_eval.py`; `tails_check.py`; `scaled_check.py`; `arb_eval.py`; `upper_check.py`;
`upper_arb.py`; `monotone_bounds_check.py`; `certified_effect.py`; `precision_check.py`; `outputs/` (all outputs named above).

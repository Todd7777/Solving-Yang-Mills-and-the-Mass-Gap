# Referee report (eleventh pass on Section 16): Corollary 16.21, $SU(2)$ and $U(1)$ on the whole range

Independent referee checks, made on 14 November 2026, of the items added in commit `a5a212a` to `theory/mass-gap-strong-coupling.tex`. I ran only read-only git
commands and wrote no file outside this folder. Bytecode caches were suppressed (`PYTHONDONTWRITEBYTECODE=1`), and the authors' `run_cap.sh 2 classic` commands
were rerun with their outputs redirected here.

The items refereed:
* Corollary `mg:cor:classic` (Corollary 16.21, lines 5340–5358) and its proof (lines 5360–5384);
* the subsection "$SU(2)$ and $U(1)$" (lines 5333–5338);
* the new sentence of the abstract (lines 148–150);
* the new row of the summary table (line 5526).

I checked Theorem 16.18, Corollary 16.20 and Lemmas 16.13–16.17 only where the corollary uses or extends them. That means $\iota=1$ ($SU(2)$) and $N=d_\rho=1$
($U(1)$), $\bar u_2=e^\beta-1-\beta$, $\varpi=b^2$, $\bar a=2^{-8}$ and the larger $F$ raise.

Nothing here bears on the Yang–Mills problem itself, which is open. The corollary concerns lattice gauge theory at strong coupling, $\beta<\beta_0(d)$, with
$\beta_0(3)\approx4.0\cdot10^{-4}$ and $\beta_0(4)\approx1.5\cdot10^{-4}$. For $U(1)$ a gap at strong coupling is expected (the massless phase of Guth and
Fröhlich–Spencer lies at weak coupling), so the README's $U(1)$ test is respected. Nothing in the refereed text claims or implies more.

## Summary

The corollary is correct, and every constant in it is certified here in ball arithmetic. No error affects a stated result. Each claim was checked by an
independent argument or program:

* **Applicability (item 1).** Lemma 16.14, Lemma 16.15 and the sector upper bound of Corollary 16.20 hold for $(N,\iota)=(2,1)$ and $(1,0)$. Every input is
  general in $(N,\iota)$:
  * $u_*=u$ (Corollary 6.5, $N\le49$; the condition $\beta e^{2\beta}\le(1+\iota)/N^2$ holds);
  * $u_2\le e^\beta-1-\beta$ and $\sigma_6\le(2-\iota)N^2u^6+u_2^5(e^\beta-1)$ (Lemma 9.1(iv));
  * $u_{\rm low}\le u\le\bar u$ (Lemma 5.1(iii): $m_\rho(0,1)=\iota$, $m_\rho(n,m)\le N^{n+m-1}$);
  * $e^{-2\beta}\le\Upsilon,\Upsilon_\pm\le e^\beta$ (Lemma 16.13(i));
  * Theorem 15.5(iii) for $U(1)$ ($\iota=0$, $\mathsf c(U)=\bar U$), and Lemmas 10.6, 15.3;
  * $\sqrt{2-\iota}\,N$ and $N^4/\tau_0^2$ in Lemma 10.6, Lemma 16.8 and Theorem 10.2.

  Lemma 16.17 is rightly not used: $\Upsilon<1$ for $SU(2)$, and $\Upsilon_\pm<1$ for $U(1)$. A numerical sanity check with the exact Bessel-function group data
  confirms every bound-mode input, the dressed-tube bound $\bar\Upsilon_{(C)}\le\Upsilon_{(C)}+\Delta$ and $a_1=2u,\sqrt2u,0$. It also shows that bound mode
  costs only 0.06–0.33% of $D$ (§A).
* **Code (item 2).**
  * `group_data` with `--group 1` and `--group 2` calls `bound_group(1, 0, sharp_u2=False)` and `bound_group(2, 1, sharp_u2=False)`.
  * `lower3e` honours `F_RAISE` and `KMAX_DEFAULT`, and its `strong` dict computes exactly the stronger form of Lemma 16.15.
  * `sector_upper.py` uses $\sqrt{2-\iota}$ and $\tau_0=(1+\iota)/2$, with $C=+$ only for $SU(2)$.

  All ten `classic` outputs of `run_cap.sh` reproduce identically, apart from the timing lines I added (§G). They take about 3 minutes on this machine, not an hour.
* **The two technical points (item 3).** Both are right, and both are needed.
  * With the default raise $10^{-9}$, the stronger $F$-condition for $U(1)$, $d=3$, fails (room $-3.65\cdot10^{-9}$). With $10^{-7}$ it holds: room
    $9.535\cdot10^{-8}$ in Arb, against $\delta'_*=4.649\cdot10^{-9}$.
  * With fronts to 60 edges in $d=4$, the barred tail is $1.09\cdot10^{-13}$, and four barred inequalities fail (down to $-1.09\cdot10^{-14}$). With 100 edges
    the tail is $1.77\cdot10^{-25}$ and the room $\ge2.503\cdot10^{-16}$. My Arb evaluation sums to 300 edges with rigorous remainders.
* **The series (item 4).** I recomputed them from closed forms in modified Bessel functions, not from Haar moments:
  $Z_{SU(2)}=2I_1(\beta)/\beta$, $\Upsilon=4(\log Z)''$, $\Upsilon_+=2(\log I_0)''$ and $\Upsilon_-=2I_1/(\beta I_0)$.
  * They are exact to $\beta^{10}$, and every printed coefficient is right.
  * High-precision numerics confirm the $\beta^{10}$ coefficients.
  * The factor $e^{-3\beta}$ in the splitting bound is right: $I_0\le e^\beta$ and $\Upsilon_-\le e^\beta$.

  Moreover the coefficients of $I_1^2-I_0I_2$ have a positive closed form. So $I_1^2-I_0I_2\ge(\beta/2)^2/2$ holds for every $\beta>0$, and that step needs no
  computer at all (Issue 5).
* **The numbers (item 5).** All 33 quoted numbers are correct and on the safe side, each the certified value rounded up at its quoted precision
  (`outputs/numbers_check_classic.out`).
* **Ball arithmetic (item 6).** I adapted the ninth and tenth referees' evaluators to general $(N,\iota)$, $\bar u_2=e^\beta-1-\beta$, the authors' $\varpi=b^p$,
  KMAX and $F$ raise.
  * Every hypothesis at $\beta_*=\beta_0(d)(1+10^{-30})$ is certified in all four cases, including the stronger form of Lemma 16.15 and $\varpi_*=\hat b_*^p$.
  * $(1+\kappa)^3D(\beta_*)/\beta_*^4$ is below the table's $c_-$ in all four cases (table in §B).
  * $E_C/\beta^4$ at an exact point above $\beta_0(d)$ is below $c_+$ in all six cases (three per dimension).
* **Labels (item 7).** "Computer-assisted" is justified once this certification is cited (Issue 4). "Determined" is used in the paper's sense (residual smaller
  than the term). That sense has been used before with residuals up to 84% (Theorem 16.18, $SU(5)$, $d=3$), so 23% for $U(1)$ in $d=3$ is honest. But it should
  be visible in the statement, and the subsection's "far below" overstates it (Issue 2).
* **A finding beyond the brief (Issue 1).** In $d=3$ the large lower constants come mainly from the crude count of fronts with more than 12 edges. The larger
  $a$ of $SU(2)$ and $U(1)$ exposes that count; it is not the powers $a^E$ alone, as the text says. $d=3$ uses $2\binom{6k}{k-1}\sim14.9^k$ there, although
  Lemma 16.10 proves $N^v_k\le0.783\cdot4.081^k$ for $\mathbb Z^2$ too. With that bound the constants fall (for $SU(3)$ only by a factor 3):
  * $SU(2)$, $d=3$: $c_-=1968$ becomes $12.7$;
  * $U(1)$, $d=3$: $c'_-=1.772\cdot10^5$ becomes $181$.

  Both are certified here in Arb at $\beta_*$, and the scaling check passes. The $U(1)$, $d=3$ residual then drops from 23% to 0.023% of $\beta^2/8$.

## Recommendation

**Minor revisions.** The corollary stands with the stated constants. The "Computer-assisted" label is justified once the ball-arithmetic evaluations of this
folder (or equivalent) are cited, as was done for Theorem 16.18 and Corollary 16.20. I strongly recommend Issue 1: it is a small change with a large effect
on the $d=3$ constants.

## Issues

Line numbers refer to `theory/mass-gap-strong-coupling.tex` as of commit `a5a212a`.

1. **Minor (the explanation of the large constants is incomplete in $d=3$; a proven bound is left unused).** Lines 5336–5337: "Its constants are much larger
   than for $SU(3)$, since $u=\beta/4$ and $\beta/2$ here, against $\beta/18$, enter the junction sums through powers $a^E$, $E\ge2$." The same explanation is in
   the journal (Part 24).
   * **In $d=4$ it is right.** The one-plaquette term $4P_vb(\mathsf F_8\mathsf A+\mathsf F\mathsf A_8+\mathsf F(e^{\mathsf A}-1-\mathsf A))$, of order
     $ba^4$, is 78% of $\delta'_*$ for $SU(2)$ and 88% for $U(1)$.
   * **In $d=3$ the main cause is the crude count.** `practical_beta.load_terms` keeps `NV_DECISION = False` in $d=3$ ("as refereed"). So beyond the exact
     counts ($k\le12$) the fronts of $\mathbb Z^2$ are counted by Definition 8.9's $(d-1)\binom{(4d-6)k}{k-1}$. That bound inflates the one-plaquette sums, and
     the barred sums force $\bar a=2^{-10}$, so the corner term $(a/\bar a)^6$ at $(3,6)$ dominates: 1586 of $1963\beta^4$ for $SU(2)$, $1.41\cdot10^5$ of
     $1.77\cdot10^5$ for $U(1)$.
   * **Lemma 16.10 already covers $D=2$**: $N^v_k\le0.783\cdot4.081^k$. With it in $d=3$ (and fronts to 100 edges, as $\bar a=2^{-6}$ then):

     | case | $D(\beta_*)/\beta_*^4$ | $(1+\kappa)^3D/\beta_*^4$, Arb, $\le$ | stated $c_-$ |
     |---|---|---|---|
     | $SU(2)$, $d=3$ | $1964.71\to12.6726$ | $12.68791$ | 1968 |
     | $U(1)$, $d=3$ | $176920\to180.613$ | $180.8309$ | $1.772\cdot10^5$ |

     * Every hypothesis holds, the stronger form included (`outputs/arb_eval_decision_d3_g*.out`), and the 17-point scaling check passes
       (`outputs/scaling_decision_d3_g*.out`).
     * The $U(1)$, $d=3$ residual at $\beta_0(3)$ becomes 0.023% of $\beta^2/8$ instead of 23%; for $SU(2)$, 0.0016% instead of 0.26%.
     * The $F$ raise $10^{-7}$ is then no longer needed, as $\delta'_*=4.7\cdot10^{-12}$.
     * The same change gives $0.1845$ for $SU(3)$, $d=3$, instead of $0.5993$ (mpmath only, outside this report's scope).
   * **Recommendation.** Adopt it, with an Arb evaluation as here. At the least, correct the explanation.

2. **Minor (honesty of "far below" and "determined").**
   * Line 5338: "They are still far below the second-order terms on the whole range." For $U(1)$ in $d=3$ the residual reaches $0.23\cdot\beta^2/8$ at
     $\beta_0(3)$: below the term, not far below. (It is 11.5% of the splitting $\beta^2/4$.)
   * Lines 5356–5357: "are determined from both sides on the whole range". The sense (residual smaller than the term) and the ratios are only in the proof
     (lines 5382–5383). Remark 16.9(iii) set the better precedent: "45% of the term at $\beta_0(3)$, 11% at $10^{-4}$". State the ratios in the corollary, or
     in a sentence after it. Issue 1 would make this moot.

3. **Minor (overstatement).** Line 5335: "so the gap has no terms of odd order beyond $4\log(1/u)$".
   * Only orders 1 and 3 are established: Section 9's $\kappa=O(\beta^2)$, and $-\log\Upsilon_{(C)}$ even with the $O(\beta^4)$ bound. The expansion of $m$ is
     proved only to $o(\beta^4)$, and its analyticity is listed as open (summary table).
   * Write "no first- or third-order term".
   * The reason for the evenness is $\rho(-1)=-1$ with $-1$ central, not centrality alone. Lemma 11.1(i) proves it: cite that lemma instead of the
     parenthesis "($-1$ is central)".

4. **Minor (what the evaluation covers; cite the certification).**
   * Lines 5367 and 5372 cite only `run_cap.sh` and `sector_upper.py`. That is, an mpmath evaluation at $\beta_*$ for the lower bound and at $\beta_0(d)$ (not
     above it) for the upper bound.
   * `error_budget.py` bounds only the families with the least room, as the ninth and tenth referees noted, and it evaluates at $\beta_0(d)$, not at $\beta_*$.
   * The gaps are closed by this folder's Arb runs (§B, §C).
   * Cite them, and add a "What is computer-assisted" paragraph to the proof, as in Corollary 16.20 (line 5326). The new computer steps are the four evaluations
     at $\beta_*$, the six upper bounds, the exact series and the bound on $I_1^2-I_0I_2$.

5. **Minor (the splitting bound: a small gap in the code, and a proof without computer).** Lines 5379–5380, `rigorous_bounds.u1_split_lower`.
   * The "tail bounded by the products of the non-negative majorant series" is computed by `tailprod`, which truncates each factor at $K+60$ terms and bounds
     nothing beyond. The neglected part is below $10^{-940}$ relative (its largest term is $3\cdot10^{-944}$), but it is not bounded. (`split_lower` for $SU(N)$ at least carries a factor
     $1+10^{-20}$.)
   * Better, the step needs no computer. By the Cauchy product and Vandermonde's identity,
     $I_\mu I_\nu=\sum_k(\beta/2)^{\mu+\nu+2k}(\mu+\nu+2k)!/(k!(\mu+\nu+k)!(\mu+k)!(\nu+k)!)$. Hence
     $$I_1^2-I_0I_2=\sum_{k\ge0}\frac{(2k+2)!}{k!\,(k+2)!\,((k+1)!)^2\,(k+2)}\Bigl(\frac\beta2\Bigr)^{2k+2},$$
     all terms positive (checked against the exact coefficients to $k=60$, `outputs/series_check_classic.out`).
   * So $I_1^2-I_0I_2\ge(\beta/2)^2/2$ and $\log(\Upsilon_-/\Upsilon_+)\ge e^{-3\beta}\beta^2/4$ for every $\beta>0$, without the $(1-10^{-7})$ and without
     the restriction to $(0,4.03\cdot10^{-4}]$.

6. **Minor (statement form).** Line 5354 chains an asymptotic expansion and an inequality: "$\log(\Upsilon_-/\Upsilon_+)=\frac{\beta^2}4-\frac{\beta^4}{48}+O(\beta^6)\ge
   e^{-3\beta}(1-10^{-7})\frac{\beta^2}4$". The inequality is a statement for $0<\beta\le4.03\cdot10^{-4}$ (or every $\beta>0$, by Issue 5). Split the two.

7. **Minor (consistency of $B_C$).**
   * Line 5371: "The bound $E_C$ of that proof, with $a=\sqrt{2-\iota}\,N\bar u$". The displayed $B_C$ of Corollary 16.20 (line 5308) has $\sqrt2N\bar u$ in its third term.
   * `sector_upper.py` uses $\sqrt{2-\iota}$ there too, which for $SU(2)$ is $1$ (it is $\norm{P_\rho f_\beta}/(Nu)$, Theorem 10.2).
   * Both are valid, and the difference is $6\cdot10^{-9}$ relative (`outputs/upper_arb_classic.out`). Say "$\sqrt{2-\iota}$ in place of $\sqrt2$ throughout".

8. **Minor (wording).** Line 5364: "$F_*,F_{6*}$ are raised by $10^{-7}$ instead of $10^{-9}$". The raise is relative, $F_*=F_{\min}(1+10^{-7})$ in
   `lower3e`. Write "by the factor $1+10^{-7}$". Optionally add that with $1+10^{-9}$ the stronger form fails (room $-3.65\cdot10^{-9}$).

9. **Minor (the comparison at $\beta_0(d)$).** Lines 5381–5382: "which exceeds $(c'_-+c'_+)\beta^4$ on the whole range: at $\beta_0(d)$ ...".
   * Say why the check at $\beta_0(d)$ suffices: $(c'_-+c'_+)\beta^2$ increases and $e^{-3\beta}(1-10^{-7})/4$ decreases. Corollary 16.20 says this for $SU(5)$.
   * "$0.2497$" is the $d=3$ value $0.249698$, rounded to nearest on the side the comparison must win. This is harmless: the margin is a factor 8.7.

10. **Minor (abstract).** Lines 148–150.
    * The old window "$-9.4\cdot10^{-5}\le m-4\log(1/u)\le4.8\cdot10^{-4}$" was labelled "a numerical evaluation of the explicit constants" in the old abstract,
      and is Numerical in the table (lines 5481, 5510). The new sentence drops the label: add "(numerical)".
    * "the same evaluation gives" could say that the new bounds are computer-assisted.

11. **Minor (table row, line 5526).**
    * It writes $(c_-,c_+)$ for $U(1)$, where the corollary has $c'_\pm$.
    * "$=O(\beta^4)$ on the whole range" would read better as "between $-c_-\beta^4$ and $c_+\beta^4$ on the whole range".
    * Optionally, mark the old $SU(2)$ whole-range row (line 5481) as superseded.

12. **Cosmetic (code).**
    * `rigorous_bounds.py` prints no bound-mode data ($\bar u/\beta$, $u_{\rm low}/\beta$, $\Delta/\beta^4$) for $N=1,2$, so `orders/rigorous_bounds.out`
      does not document the inputs of the new cases. (The values are in `outputs/applicability_check.out`.)
    * `bound_group(1, 0)` names itself "SU(1) bounds".
    * The header of `run_cap.sh` gives "about 15 + 15 + 25 minutes"; the classic part takes about 3 minutes here (4 cores).

## What I checked, and how

### A. Applicability (item 1)

* **The lemmas, read for $\iota=1$ and $N=1$.** Each input below is stated and proved for general irreducible non-trivial $\rho$, and nothing assumes $N\ge3$:
  * Lemma 5.1(iii): $m_\rho(0,1)=\iota$; for $U(1)$, $m_\rho(n,m)=[n-m=1]\le1=N^{n+m-1}$;
  * Corollary 6.5: $u_r\le e^\beta-1-\beta$ for $r\notin\{0,\rho,\bar\rho\}$;
  * Lemma 9.1(iv);
  * Lemma 10.1(i): $\iota=1$ gives the eigenvector $\chi_\rho$, $\Upsilon=1+\kappa-N^2u^2$;
  * Lemma 10.6: $\norm{P_\rho v}=\sqrt{2-\iota}\,Nu^3$; for $\iota=1$, $P_\rho f_\beta=Nu\chi_\rho$;
  * Theorem 10.2 ($a_1=|\langle f_\beta\psi\rangle|$, $\norm{P_\rho f_\beta}=\sqrt{2-\iota}\,Nu$, $N^4/\tau_0^2$), with Lemma 5.1(iv);
  * Lemma 16.8 ($\sigma_5$), Lemma 16.13(i) ($|\mathrm{Re}\,\chi_\rho|\le N$);
  * Lemma 16.14 ($\bar u_2$ "in general");
  * Lemma 16.15: $\bar u_2/\beta^2=\sum_j\beta^j/(j+2)!$ increases; $1/u_{\rm low}$ enters as $e^\beta/\beta$;
  * Theorem 15.5(iii) and Lemma 15.3 ($\iota=0$ and $\mathsf c(U)=\bar U$ for $U(1)$).

  Lemma 14.1 and Theorem 14.2 (the moves by characters, `chars=True`) are group-free given $u_*=u$.
* **Numerical sanity check** (`applicability_check.py`, Bessel functions, at $\beta_0(3)$, $\beta_0(4)$, $10^{-4}$, $10^{-6}$). Exact data:
  * $SU(2)$: $u_j=I_{2j+1}/I_1$, $\kappa=3I_3/I_1$;
  * $U(1)$: $u_n=I_n/I_0$, $\alpha_\rho=0$, $\gamma_\rho=I_2/I_0$.

  Results, all as expected:
  * $u_*=u$ and $u_{\rm low}\le u\le\bar u$;
  * $u_2\le e^\beta-1-\beta$ (true $u_2/\beta^2=1/24$, $1/8$, against $0.5$);
  * $\sigma_6\le\bar\sigma_6$ and $\sigma_5\le\bar\sigma_5$;
  * $e^{-2\beta}\le\Upsilon_{(C)}<1$;
  * Lemma 10.1(i) agrees with the direct variance to $10^{-30}$;
  * $\Upsilon_{(C)}\le\bar\Upsilon_{(C)}\le\Upsilon_{(C)}+\Delta$ with Lemma 10.6/15.3 at the exact data;
  * $n_B\le\bar n_B$ and $c_2\le\bar c_2$;
  * $a_1=2u$, $\sqrt2u$, $0$;
  * $\Upsilon_->\Upsilon_+$.
* **Exact data against bound mode** (`exact_vs_bound.py`). $D(\beta_*)$ with the exact $u,\Upsilon,\bar\Upsilon$ is lower by 0.20% ($SU(2)$, $d=3$), 0.064%
  ($d=4$), 0.33% ($U(1)$, $d=3$) and 0.094% ($d=4$). So the proven bounds cost almost nothing, and the size of $c_-$ is not an artefact of bound mode.

### B. The lower bounds in ball arithmetic (items 3, 6)

`arb_eval_classic.py` is the tenth referee's `arb_eval_cor.py` (itself the ninth referee's `arb_eval.py`), with the changes marked `16k`:
* $(N,\iota)=(2,1)$ and $(1,0)$; the earlier code hard-coded $\iota=0$ and $N=$ GROUP;
* $\bar u_2=e^\beta-1-\beta$;
* the authors' $p$ in $\varpi=\hat b_*^p$, rechecked in Arb;
* their $(\bar b,\bar a)$, KMAX 100 in $d=4$ and $F$ raise $10^{-7}$ for $U(1)$, $d=3$;
* the barred tree constants at the authors' KMAX, not a fixed 60;
* fronts to 300 edges.

As before: 400-bit balls, clusters and polymers to 110 plaquettes, rigorous remainders, the referee's own marked constants (six inequalities verified in Arb),
and $\beta_*\ge\beta_0(d)$ checked in balls.

| case | $\varpi$, $(\bar b,\bar a)$ | every hypothesis | stronger form ($F$ room) | $\delta'_*$ | $D(\beta_*)/\beta_*^4$ (= authors' to 12 digits) | $(1+\kappa)^3D/\beta_*^4\le$ | $c_-$ |
|---|---|---|---|---|---|---|---|
| $SU(2)$, $d=3$ | $b^{5/2}$, $(2^{-7},2^{-10})$ | holds | holds ($9.48\cdot10^{-10}$) | $5.16\cdot10^{-11}$ | $1964.70961701$ | $1967.083675$ | 1968 |
| $SU(2)$, $d=4$ | $b^{5/2}$, $(2^{-8},2^{-8})$ | holds | holds ($1.00\cdot10^{-9}$) | $3.18\cdot10^{-13}$ | $568.907134026$ | $569.1696109$ | 570 |
| $U(1)$, $d=3$ | $b^2$, $(2^{-8},2^{-10})$ | holds | holds ($9.54\cdot10^{-8}$) | $4.65\cdot10^{-9}$ | $176919.852420$ | $177133.6336$ | $1.772\cdot10^5$ |
| $U(1)$, $d=4$ | $b^2$, $(2^{-8},2^{-8})$ | holds | holds ($9.96\cdot10^{-10}$) | $4.49\cdot10^{-12}$ | $8038.46664120$ | $8042.175351$ | 8043 |

The other rooms:
* decoration conditions $\ge9.7\cdot10^{-16}$, barred $\ge6.5\cdot10^{-16}$;
* conditions on $\delta'$, $F$, $F_6$ $\ge10^{-9}$;
* third and fourth hypotheses $\approx1$, also with $1/(1+\kappa)$ and $\varpi_*/(1+\kappa)$;
* $\delta'_*\le\kappa$.

The parts of $\Sigma^{(4)}_4/\beta^4$:

| case | two | four | exits | one plaquette | multi (corner) |
|---|---|---|---|---|---|
| $SU(2)$, $d=3$ | 0.456 | 0 | 2.70 | 373.6 | 1586.4 $(3,6)$ |
| $SU(2)$, $d=4$ | 1.12 | 0.18 | 79.6 | 445.9 | 42.0 $(8,3)$ |
| $U(1)$, $d=3$ | 32.2 | 0 | 21.8 | 35468 | 141255 $(3,6)$ |
| $U(1)$, $d=4$ | 128.3 | 0.77 | 542.9 | 7035 | 329.1 $(8,3)$ |

The two technical points, tested with the authors' `error_budget.py` (both must fail with the defaults):
* $d=4$ with `--kmax 60`: the barred tail is $1.09\cdot10^{-13}$, and jc_bar3–6 fail (rooms $-5.8\cdot10^{-16}$ to $-1.09\cdot10^{-14}$). With
  `--kmax 100`: tail $1.77\cdot10^{-25}$, room $\ge2.503\cdot10^{-16}$. So the proof's "tails at most $2\cdot10^{-25}$, room at least $2.5\cdot10^{-16}$" is right.
* $U(1)$, $d=3$ with the default raise: strong_F $=-3.65\cdot10^{-9}$. With `--f-raise 1e-7`: $9.54\cdot10^{-8}$.

In my Arb run at $d=4$, the barred $\xi_F$ condition alone would still hold at 60 edges (room $6.29\cdot10^{-16}$). It is the $k$-weighted sums
$\bar A_1,\bar A_2$ of the six barred inequalities that need the 100 edges.

### C. The upper bounds in ball arithmetic (items 1, 6)

`upper_arb_classic.py` adapts the tenth referee's `upper_arb_sectors.py` to general $\iota$; its power-series check had hard-coded $\iota=0$. It also contains my
own derivation for $SU(2)$, with $\psi=\chi_\rho$ and $a_1=Nu$:
* $V_2\le(1+\mathcal T)\Upsilon+\mathcal T(1+2a_1+3a_1^2)$ in both cases of the case analysis;
* $L_2/\Upsilon^2\ge1-B$ with $B$ as in `sector_upper.py` ($\sqrt{2-\iota}=1$, $\tau_0=1$, $e^{20\beta}$ from $\Upsilon\ge e^{-2\beta}$).

$E_C/\beta^4$, certified at an exact number above $\beta_0(d)$:

| case | $E_C/\beta^4$ (ball upper end) | $c_+$ |
|---|---|---|
| $SU(2)$, $d=3$ | $0.0025082343$ | 0.0026 |
| $SU(2)$, $d=4$ | $0.0027671106$ | 0.0028 |
| $U(1)$, $d=3$, $C=+$ / $-$ | $0.0016790406$ / $0.0016788876$ | 0.0017 |
| $U(1)$, $d=4$, $C=+$ / $-$ | $0.0017947730$ / $0.0017947287$ | 0.0018 |

* All agree with `sector_upper.E_sector` to $5\cdot10^{-60}$ relative.
* $B_C\le6.2\cdot10^{-17}$, so $L_2>0$.
* $E_C/\beta^5$ is a power series with non-negative coefficients, by the structural argument of Lemma 16.17's proof. My exact expansion to $\beta^{14}$
  confirms it (the coefficients of $\beta^0..\beta^4$ vanish). So the value above $\beta_0(d)$ bounds the whole range.

### D. The series and the splitting (item 4)

* **The series** (`series_check_classic.py`, Fractions from the Bessel series, my own power-series log and derivatives). Results, to $\beta^{10}$:
  * $-\log\Upsilon_{SU(2)}=\frac{\beta^2}8-\frac{\beta^4}{192}+\frac{11\beta^6}{46080}-\frac{17\beta^8}{1474560}+\frac{37\beta^{10}}{61931520}$;
  * $-\log\Upsilon_+=\frac38\beta^2-\frac{13}{384}\beta^4+\dots$;
  * $-\log\Upsilon_-=\frac18\beta^2-\frac5{384}\beta^4+\dots$;
  * $\log(\Upsilon_-/\Upsilon_+)=\frac{\beta^2}4-\frac{\beta^4}{48}+\frac{\beta^6}{512}-\dots$

  All odd coefficients vanish, and the series are identical to `classic_series` to $\beta^8$. The two closed forms of $\Upsilon_\pm$ agree exactly. At
  $\beta=0.1,\dots,10^{-4}$ (80-digit mpmath), (exact $-$ series to $\beta^8$)$/\beta^{10}$ converges to the exact $\beta^{10}$ coefficient (6 digits
  for $\beta\le4\cdot10^{-4}$).
* **The Turán determinant.** $p_k$ for $k\le60$ equals the closed form of Issue 5, and all are positive. The authors' route (subtracting
  $\sum_{k\ge1}|p_k|x^k$) gives $(1/2)(1-2.707\cdot10^{-8})$, matching their $3\cdot10^{-8}$. My Arb tail bound is the geometric majorant of
  $2(2^k/k!)^2x^k$.
* **The chain** $\log(\Upsilon_-/\Upsilon_+)\ge(\Upsilon_--\Upsilon_+)/\Upsilon_-\ge2e^{-3\beta}P\ge e^{-3\beta}(1-10^{-7})\beta^2/4$ holds:
  * $\log x\ge1-1/x$;
  * $I_0\le e^\beta$ and $\Upsilon_-\le e^\beta$ (Lemma 16.13(i)) give $e^{-2\beta}\cdot e^{-\beta}$.

  The chain is also confirmed numerically at 2000 points of $(0,4.03\cdot10^{-4}]$; there $\min\log(\Upsilon_-/\Upsilon_+)/(\beta^2/4)=1-1.35\cdot10^{-8}$.
* **The consequences.** $m_+-m_-\ge\log(\Upsilon_-/\Upsilon_+)-(c'_-+c'_+)\beta^4>0$, so $m=m_-$. The comparison numbers are $0.028725$ ($d=3$) and
  $1.9017\cdot10^{-4}$ ($d=4$), against $0.249698$ and $0.249885$.

### E. The labels and the claims (item 7)

* **The corollary rests on three kinds of computer step.**
  * The enumerations inherited from Theorem 16.18 and Corollary 16.20, each by two programs. They are unchanged here: the same `d3_s9z.json`,
    `d4_s10x_92.json` and closed-polymer counts.
  * Exact rational computations (the series, the splitting), now with a rigorous tail or none (Issue 5).
  * One evaluation per case, now certified in ball arithmetic (§B, §C).

  No grid or numerical monotonicity is used. So "(Computer-assisted.)" and "Theorem (computer-assisted)" are justified, with Issue 4.
* **"Determined"** has the paper's meaning. At $\beta_0(d)$ the residual-to-term ratios are 0.26%, 0.011% ($SU(2)$, $d=3,4$) and 23%, 0.15% ($U(1)$, of
  $\beta^2/8$), and they fall like $\beta^2$. Honest, but see Issue 2.
* **"The third-order terms are zero"** is an asymptotic statement, correct ($-\log\Upsilon_{(C)}$ even, residual $O(\beta^4)$).
* **No overclaim.** Nothing claims or implies that the Yang–Mills mass gap problem is solved. The abstract sentence, the subsection, the corollary and the table
  row concern lattice theory at $\beta<\beta_0(d)$.
* **Unlabelled claims.** The subsection's "no terms of odd order" (Issue 3), and the abstract's quotation of a Numerical window without its label (Issue 10).

### F. The numbers (item 5)

`numbers_check_classic.py` checks 33 numbers against the outputs above: the constants, the bounds on $D$, the tails and rooms, the $E_C$, $3\cdot10^{-8}$, the
comparisons, the residual ratios, the abstract's $570$ and $0.0028$, and the subsection's $u_1$. All 33 are correct.
* Every bound is the certified value rounded up at its quoted precision (for example $569.17\ge569.16961$ and $1.6\cdot10^{-3}\ge1.521\cdot10^{-3}$).
* "$\delta'_*\approx5\cdot10^{-9}$" is $4.65\cdot10^{-9}$ to one digit.
* "$0.2497$" is the only rounding on the side a comparison must win (Issue 9).

### G. Reproductions (item 5)

`run_reproductions.sh` reruns the `classic` part of `run_cap.sh`, with the outputs in `outputs/repro_*.out`:
* the four scalings, the four error budgets, `sector_upper.py` and `rigorous_bounds.py`;
* all ten are identical to the files in `../orders/`, apart from the timing lines;
* about 3 minutes with 2 jobs.

## What I could not check

* **The enumerations.** The counts inherited from Theorem 16.18 and Corollary 16.20 (pairs, corners, exits, $Z_{n,4}$, $N^v_k$, closed polymers, admissible
  sets) I did not recount. Earlier referees checked each of them by two programs, and the corollary adds no count.
* **The general theorems.** Sections 10–15 and Lemmas 16.14–16.16 I read only as far as the corollary uses them, checking where $N$, $\iota$, $\tau_0$ and
  $\bar u_2$ enter. I did not re-derive them.
* **The literature.** I did not compare the $U(1)$ second-order terms with the classical strong-coupling literature: no source was re-read here.

## Files

Run from this folder: `bash run_checks.sh` (about 3.5 min on 4 cores). From `theory/checks/third-order`: `bash referee-section16k/run_reproductions.sh 2`
(about 3 min).

| file | what |
|---|---|
| `run_reproductions.sh` | the `classic` part of `run_cap.sh`, outputs to `outputs/repro_*.out` (identical to `../orders/`) |
| `run_checks.sh` | all my checks below, and the two failure tests of the technical points |
| `arb_eval_classic.py` | §B: the lower bound in Arb, four cases, stronger form of Lemma 16.15; `--decision-d3` for Issue 1; GROUP 3, 4, 5 for the addendum (Lemma 16.11 for $\bar u_2$) |
| `upper_arb_classic.py` | §C: $E_C$ in Arb above $\beta_0(d)$, six cases; my $2a$ variant and the $\sqrt2$ variant; against `sector_upper.py`; power-series structure for general $\iota$ |
| `series_check_classic.py` | §D: the series from Bessel closed forms (exact), numerics, the Turán coefficients and their closed form, the chain |
| `applicability_check.py` | §A: the bound-mode inputs and the dressed tube with the exact $SU(2)$, $U(1)$ data (Numerical) |
| `exact_vs_bound.py` | §A: $D(\beta_*)$ with exact group data against bound mode (Numerical) |
| `nv_decision_d3_probe.py` | Issue 1: $D(\beta_*)$ in $d=3$ with Lemma 16.10 for $\mathbb Z^2$ (authors' pipeline, mpmath), $SU(2)$, $U(1)$, $SU(3)$ |
| `scaling_decision_d3.py` | Issue 1: the authors' `scaling_check.py` in that setting (a wrapper; nothing in the repository is modified) |
| `numbers_check_classic.py` | §F: the 33 quoted numbers against the outputs |
| `outputs/` | all outputs, each with its command line |

## Addendum (follow-up request): Lemma 16.10 in $d=3$ for $SU(3)$, $SU(4)$, $SU(5)$ (Corollary 16.20)

The coordinator adopted Issue 1 everywhere in $d=3$, with a new option `--nv-decision-d3` in `practical_beta.py`, `scaling_check.py` and `error_budget.py`.
That option sets `NV_DECISION` in `load_terms`, which is what `--decision-d3` does here. Theorem 16.18's own setting, the default, is unchanged. I was asked to
certify the three $SU(N)$ cases of Corollary 16.20 in $d=3$ in ball arithmetic.

**The extension.** `arb_eval_classic.py` now takes GROUP $=3,4,5$, with $\iota=0$ and $\bar u_2$ of Lemma 16.11 ($c_\Lambda=2$ for $N=4$, no $\Lambda^2$ entry for
$N=3$), exactly as in `referee-section16j/arb_eval_cor.py`. I reread the lemma and `sun_group.u2_bound` for this. As a check of the extension, without the option
it reproduces the tenth referee's certified values for `d3_s9z.json` to all printed digits: $0.5993437096$, $0.1332401769$, $0.07288809269$
(`outputs/arb_eval_validate_su*_d3_s9z.out`).

**The certification.** The setting is that of the authors' full search (`orders/capz_su{3,4,5}_d3_scaling.out`): `d3_s9z.json`, $(\bar b,\bar a)=(2^{-7},2^{-6})$,
$\varpi=b^{5/2}$, fronts to 100 edges, $F$ raise $10^{-9}$. The Arb evaluation sums fronts to 300 edges and clusters/polymers to 110 plaquettes, with rigorous
remainders, and checks $\beta_*\ge\beta_0(3)$ in balls (`outputs/arb_eval_decision_d3_g{3,4,5}.out`):

| case | $D(\beta_*)/\beta_*^4$ (= authors' to 12 digits) | $(1+\kappa)^3D(\beta_*)/\beta_*^4\le$ (Arb) | requested $c'_-$ | $\delta'_*$ |
|---|---|---|---|---|
| $SU(3)$, $d=3$ | $0.184249751630$ | $0.1844723899$ | $0.185$ | $4.1\cdot10^{-15}$ |
| $SU(4)$, $d=3$ | $0.0930582734620$ | $0.09317072049$ | $0.094$ | $1.7\cdot10^{-15}$ |
| $SU(5)$, $d=3$ | $0.0661361900644$ | $0.06621610577$ | $0.067$ | $1.0\cdot10^{-15}$ |

* **Every hypothesis holds in all three cases.**
  * decoration conditions (room $\ge9.7\cdot10^{-16}$) and barred ones ($\ge7.6\cdot10^{-16}$);
  * the six inequalities, actual and barred, with the referee's own marked constants;
  * the conditions on $\delta'$, $F$, $F_6$ (room $\ge10^{-9}$), and the third and fourth hypotheses (room $\approx1$);
  * $\varpi_*=\hat b_*^{5/2}$, checked in Arb.
* **The stronger form of Lemma 16.15 holds:** the conditions on $F$ and $F_6$ with $e^{-2\beta_*}-\hat a$ (room $1.000\cdot10^{-9}$, as $\delta'_*\approx10^{-15}$), the
  third hypothesis $\le1/(1+\kappa)$, the fourth $\le\varpi_*/(1+\kappa)$, and $\delta'_*\le\kappa$.
* **So the requested constants $0.185$, $0.094$, $0.067$ are certified.** By Lemma 16.15 they bound $D/\beta^4$ for all $\beta\le\beta_*$.
* **Reproductions.** The coordinator's `orders/capz_su{3,4,5}_d3_scaling.out` (all 17 scaled points pass) and `orders/error_budget_sectors_d3.out` (tails
  $\le8.53\cdot10^{-23}$, room $\ge2.717\cdot10^{-16}$) reproduce identically (`outputs/repro_nvdec_*.out`, added to `run_reproductions.sh`).
* **One rounding to fix.** The follow-up message quotes "$SU(5)$: $\le0.066216$", which is below the certified $0.06621610577$ by $1.1\cdot10^{-7}$. The output file
  prints $0.0662162$, correctly rounded up. If the text quotes the bound on $D$ to six digits, write $0.0662162$ (or $0.06622$); the table constant $0.067$ is
  unaffected. $SU(3)$'s $0.184473$ and $SU(4)$'s $0.0931708$ are correctly rounded up.
* **What is left in $D$.** With the decision bound the corner majorant becomes negligible ($0.0017$, $0.00017$, $0.00003\,\beta^4$). $D$ is now mostly three terms:
  the shifts and jumps ($0.101$, $0.057$, $0.037\,\beta^4$), the closed-polymer loss $8(K+1)\bar\theta=0.0274\,\beta^4$ (the same for all three), and the exits
  and one-plaquette junctions ($0.051$, $0.0075$, $0.0018\,\beta^4$). Further gains in $d=3$ would have to come from these.

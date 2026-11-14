# Referee report: Lemma 16.11 (`mg:lem:u2sun`) and Remark 16.12 (`mg:rem:su5`) of Part IV

State refereed: the final version, commit `79f8183` (working tree). In it the count of the admissible sets with 11 plaquettes in $\mathbb Z^4$ is
complete, and Remark 16.12 states the $d=4$ result from both sides. A first pass was made on the checkpoint `7bef055`; the issues it raised that the
update settled are marked *resolved*. File: `theory/mass-gap-strong-coupling.tex` (5199 lines). The PDF is current: it contains the final text, and
`pdftotext` shows no undefined references. Line numbers are those of this `.tex`.

Items refereed:
* Lemma 16.11 (lines 4906–4928) and Remark 16.12 (4930–4973);
* what depends on them: the abstract (136–139), the last sentence of the introduction of Section 16 (4217–4218), and the two rows of the
  summary table (5131–5132);
* the new rows of `../README.md` (`sun_group.py`, `--group N`, `orders/d4_s10x_92.json`, `admissible_split.c`) and `orders/admissible_d4_n11.out`.

Nothing here bears on the open problem. These are strong-coupling statements about lattice $SU(5)$ (and $SU(3)$) for $\beta<\beta_0(d)$, and the
numerical results are correctly labelled Numerical.

## Recommendation

**Minor revisions.** No error affects the results of Remark 16.12:
* $SU(5)$, $d=3$: $-0.346\beta^4\le m-4\log\frac1u+\log\Upsilon\le0.318\beta^4$;
* $SU(5)$, $d=4$: $-0.519\beta^4\le\dots\le0.252\beta^4$;
* $SU(3)$, $d=4$: $-3.94\beta^4\le\dots\le0.033\beta^4$.

Every number of the remark that I checked reproduces, most of them also with code of my own. In particular:
* I recounted the 11-plaquette sets for $(S_A,\{p_0\})$ by a different algorithm: $760$.
* My own evaluator of the upper bound gives $E(\beta_0(4))=0.2511\beta^4$ ($SU(5)$) and $0.0327\beta^4$ ($SU(3)$).

There is one genuine error, in the *general* statement of Lemma 16.11: it is false for $N=4$, the one case the proof overlooks, and it is never
used there. There is one small gap: the upper bound uses $\Upsilon\ge1$ for $SU(5)$, but the text never states it. The rest are imprecisions and
suggestions. One observation (issue 7) makes the $d=4$ upper half independent of the enumeration: the rooted binomial already proved in Lemma 5.5
bounds the 11-plaquette count, and with it $E(\beta_0(4))=0.360\beta^4<\beta^3/6000=1.084\beta^4$.

## Issues

1. **Error (the general statement of Lemma 16.11; no consequence for $SU(5)$).** Lines 4910–4916 and the proof at 4923–4926. For $N=4$,
   $\Lambda^2\rho$ (dimension 6) is self-conjugate ($\Lambda^2\otimes\Lambda^2\to\Lambda^4=\det$), so it occurs both in $\rho\otimes\rho$ and in
   $\bar\rho\otimes\bar\rho$: $m_{\Lambda^2}(2,0)=m_{\Lambda^2}(0,2)=1$. Its order-2 coefficient is therefore $t^2/2+t^2/2=t^2=\beta^2/(4N^2)$, not
   $\beta^2/(8N^2)$. The sentence "$\bar\rho\otimes\bar\rho$ is the conjugate of $\rho\otimes\rho$; for $N=3$, $\Lambda^2\rho\cong\bar\rho$" misses this
   coincidence. Numerically, by Weyl integration (`weyl_sun.py 4 40`, `weyl_sun_N4.out`), the multiplicities $[m(2,0),m(1,1),m(0,2)]$ of
   $\Lambda^2$ are $[1,0,1]$, and the bound fails:

   | $\beta$ | $u_{\Lambda^2}$ (Weyl) | Lemma 16.11 bound | ratio |
   |---|---|---|---|
   | $10^{-3}$ | $2.604\cdot10^{-9}$ | $1.307\cdot10^{-9}$ | 0.502 |
   | $10^{-2}$ | $2.604\cdot10^{-7}$ | $1.348\cdot10^{-7}$ | 0.518 |
   | $0.1$ | $2.604\cdot10^{-5}$ | $1.777\cdot10^{-5}$ | 0.682 |
   | $0.3$ | $2.344\cdot10^{-4}$ | $3.037\cdot10^{-4}$ | 1.30 |

   For $N=3$ and $N\ge5$ there is no other coincidence among the irreps of order 2: $\mathrm{Sym}^2$ has highest weight $2\omega_1$ against
   $2\omega_{N-1}$ for its conjugate, $\Lambda^2$ has $\omega_2$ against $\omega_{N-2}$, and $\mathrm{ad}$ occurs only in $\rho\otimes\bar\rho$. The lemma is
   correct there (item A below).

   Fix: for $N=4$ put $\beta^2/(4N^2)$ in the $\Lambda^2$ entry, or state the lemma for $N=3$ and $N\ge5$. $SU(4)$ has no third-order term, so
   nothing else changes. The same applies to `sun_group.u2_bound` (code and docstring, the `N >= 4` branch). Table row 5131 is right as worded
   ("exact order-2 coefficient"), but it inherits the lemma's $N\ge3$.

2. **Gap (minor): $\Upsilon\ge1$ for $SU(5)$ is used but not stated.** Lines 4938 and 4953–4958, through Remark 16.9(ii) (lines 4852–4858).
   * Where it is used. The majorant $E(\beta)$ of the upper bound, and the argument that $E/\beta^5$ increases (a power series with non-negative
     coefficients), replace $\Upsilon$ by $1$. This needs $\Upsilon\ge1$.
   * Why it is not automatic. For $SU(3)$ the text proves it ($\kappa-18u^2>0$, line 4854). For $SU(5)$, $\Upsilon_+-1=\kappa-50u^2$ cancels to
     second order, and the bounds of Lemmas 9.1(ii) and 5.1(iii) do not decide its sign.
   * How it is checked now. The code asserts `Ups >= 1` (`upper_exact.E2`, at 401 points), and `monotone_check.py --group 5` checks it at 400
     points, so it holds numerically; but the remark never says so.
   * The scaling of 16.7(vii) also uses $\bar\Upsilon\ge1$ (line 4704). There only $\bar\Upsilon(\beta_*)/\bar\Upsilon(\beta)\le1+\kappa$ is needed,
     and that holds trivially.

   A sentence suffices, for example from the exact series that (i) already uses (truncation error below $10^{-180}$). Why the first term is
   $\beta^3$: with $X=\mathrm{Re}\,\chi_\rho$ and $\Lambda(s)=\log\int e^{sX}dU$, $U_+=2\Lambda''(\beta/5)$. The cumulants are $\kappa_2=\frac12$,
   $\kappa_3=\kappa_4=\kappa_6=0$ and $\kappa_5=r_5=\frac1{16}$, so $U_+=1+\beta^3/6000+O(\beta^5)$.

   My exact series (`moments_check.py`, item B):
   * $U_+-1=\frac{\beta^3}{6000}-\frac{7\beta^5}{12\cdot10^6}+\frac{\beta^7}{1.4\cdot10^9}-\frac{\beta^8}{1.92\cdot10^{10}}-\dots$
   * $U_+-U_-=\frac{\beta^3}{3000}-\frac{\beta^5}{1.2\cdot10^6}+\dots$

   Both are positive on $(0,\beta_0(3)]$, with a relative correction of $3.5\cdot10^{-3}\beta^2$, and Weyl integration gives $U_+>1$ up to $\beta=5$.

3. **Resolved by the update ("in $d=4$ from below", first pass).** On `7bef055` the $d=4$ statement was one-sided. The $SU(5)$ term is
   *negative*, so the lower half alone ($m-4\log\frac1u\ge-\frac{\beta^3}{6000}-0.519\beta^4$) did not exhibit it. The update supplies the upper half,
   $m-4\log\frac1u\le-\frac{\beta^3}{6000}+0.252\beta^4<0$, which does. Optional: a clause in (v) noting that for $SU(5)$ it is the upper half that
   shows the sign of the term (for $SU(3)$ it is the lower half).

4. **Imprecision: "with $E(\beta)/\beta^5$ increasing (checked on 401 points)" (line 4953; also "and $E/\beta^5$ increases", line 4957).** Once
   $\Upsilon\ge1$ (issue 2), the argument of Remark 16.9(ii) applies verbatim. $E$ depends on the group only through $N$, $\iota$ and $\bar u$, and
   $E/\beta^5$ is a power series with non-negative coefficients; the counts are constants. Cite that argument, and keep the 401-point check as a
   confirmation, as for $SU(3)$.

5. **Imprecision (line 4907): "this bound loses a factor of about $N^2d_r$".** The true coefficients are $u_r\approx\beta^2/(8N^2d_r)$ for
   $\mathrm{Sym}^2,\Lambda^2$ and $\beta^2/(4N^2d_A)$ for $\mathrm{ad}$, against $e^\beta-1-\beta\approx\beta^2/2$. So the factor is $4N^2d_r$ ($2N^2d_A$ for
   the adjoint). For $\Lambda^2$ of $SU(5)$ that is $1000$, which matches the "$10^{-3}$" of the lemma.

6. **Suggestion (line 4926): "Every other $r\ne0$ has $d_r\ge N$".** This is the standard fact that the non-trivial irreps of $SU(N)$ have
   dimension $\ge N$; give a reference or one line. It is not needed for the application: with $d_r\ge1$ alone, $u_r\le R\le\beta^2/2000+R/100$
   for $\beta\le3\cdot10^{-3}$, which covers $\beta<\beta_0(3)$.

7. **Suggestion (Remark 16.12(iv)): the $d=4$ upper half does not depend on the enumeration.**
   * The proof of Lemma 5.5 (`mg:lem:twopoint`) bounds the number of admissible sets with $n$ plaquettes by the rooted binomial $\binom{J'+Kn}n$
     (Part III, Lemma `al:lem:rooted`; $K=20$, and $J'=48$ for $S_Q$, $20$ for $S_A$, in $d=4$). The crude tail of Lemma 16.8 only replaces this
     binomial by its exponential majorant $(1+\frac1K)^{J'}(eK)^n$.
   * So $N_{11}:=\binom{J'+11K}{11}$, i.e. $1.042\cdot10^{19}$ ($S_Q$) and $3.021\cdot10^{18}$ ($S_A$), is admissible in Lemma 16.8 with $n_0=11$ as it
     stands. My evaluator (`upper_su5.py`) then gives $E(\beta_0(4))=0.3603\beta^4<1.084\beta^4$, with $E/\beta^5$ increasing.
   * Any $N_{11}\le7.9\cdot10^{19}$ suffices. "Up to $10^9$" (line 4957) is very conservative.
   * The enumeration ($760$, $2616$) therefore only sharpens the constant from $0.361$ to $0.252$. It would strengthen (iv) to say so: the
     conclusion then rests on a proved bound, not only on the enumeration.

8. **Suggestion (line 63, the `\date` block).** It lists the additions up to 13 November, but not the $SU(5)$ addition of 14 November (Lemma 16.11,
   Remark 16.12). On `7bef055` the introduction's "extend both to $SU(5)$" (line 4218) was premature in $d=4$; it is true now (*resolved*).

9. **Suggestion (Remark 16.12(vi)).** (vi) sharpens the refereed $SU(3)$, $d=4$ statement from $[-11.4,0.21]\beta^4$ to $[-3.94,0.033]\beta^4$. The
   update put this into table row 5132. The abstract (line 135), table rows 5127 and 5130, and Remarks 16.7(x) and 16.9(iv) still state
   $[-11.4,0.21]$, which is true but superseded. Add a forward reference, or update the abstract.

10. **Imprecision (line 4936):** "$u$, $U_\pm$ … and $\Upsilon=\max(U_+,U_-)$ are power series in $\beta$ with exact rational coefficients". In fact $u$
    and $U_\pm$ are quotients of such series (analytic, with rational Taylor coefficients), and $\Upsilon$ is their maximum. The code evaluates the
    truncated numerator and denominator series. The wording should say so; the numbers are unaffected.

11. **Suggestion (lines 4970–4972).** The $SU(3)$ value $E(\beta_0(4))=0.033\beta^4$ (`upper_exact.py 4`) has no output file. Its run is
    `reruns/su3_upper_d4_final.out` here ($0.03273$). The parenthetical "(Outputs in …)" could list an `orders/` file for it, and
    `orders/admissible_d4_n11.out`.

## What I checked, and how

While the authors' count was running, my jobs used two cores (`nice`; a third job overlapped for a few minutes). `PYTHONPYCACHEPREFIX` pointed to a
scratch directory, so that no file outside this folder changed.
Reruns of the authors' scripts are in `reruns/`. $\beta_0(3)=4.0262143352\cdot10^{-4}$, $\beta_0(4)=1.5376644804\cdot10^{-4}$.

**A. Lemma 16.11, line by line.**
* Expansion: $w_\beta=e^{(\beta/N)\mathrm{Re}\chi}=e^{t(\chi+\bar\chi)}$ with $t=\beta/(2N)$, so $c_r=\sum t^{a+b}m_r(a,b)/(a!b!)$. This is the paper's
  normalisation (Lemma 5.1(iii),(iv): $y=\frac\beta N\mathrm{Re}\,\chi_\rho$; $u_r=c_r/(d_rc_0)$, $c_r=\int\bar\chi_rw_\beta$).
* $c_0\ge1$ by Jensen. $d_rm_r(a,b)\le N^{a+b}$ by counting dimensions. $\sum_{a+b=k}1/(a!b!)=2^k/k!$, so the tail is $R/d_r$. Correct.
* Order 2, by Weyl integration (`weyl_sun.py N M`, characters as Schur polynomials by Jacobi–Trudi, nothing divided by $\Delta$; the Weyl weight
  has mean $1$, and the characters are orthonormal on the grid to $10^{-15}$). Multiplicities $[m(2,0),m(1,1),m(0,2)]$:
  * $SU(5)$: $\Lambda^2=(1,1)$: $[1,0,0]$; $\Lambda^3=\overline{\Lambda^2}$: $[0,0,1]$; $\mathrm{Sym}^2$: $[1,0,0]$; $\overline{\mathrm{Sym}^2}=(2,2,2,2)$:
    $[0,0,1]$; $\mathrm{ad}=(2,1,1,1)$: $[0,1,0]$. As claimed.
  * $SU(3)$: as claimed ($\Lambda^2=\bar\rho$).
  * $SU(4)$: $\Lambda^2$ is $[1,0,1]$ (issue 1).
* $SU(5)$ specialisation: entries $\beta^2/3000+R/225$, $\beta^2/2400+R/576$, $\beta^2/2000+R/100$ and $R/25$. The third dominates when
  $R\le\beta^2/60$, and $R\le\beta^3e^\beta/6\le\beta^3/5\le\beta^2/100$ for $\beta\le1/20$. Correct. "About $10^{-3}$ times $e^\beta-1-\beta$": correct.
* Numerically ($SU(5)$, `weyl_sun_N5.out`, 70 irreps with parts $\le4$): $u_2$ is attained at $\Lambda^2$ (and its conjugate) for every $\beta$ tried.
  Bound/true ratio: $1.0033$ ($\beta=10^{-3}$), $1.033$ ($10^{-2}$), $1.17$ ($0.05$), $1.36$ ($0.1$), $7.5$ ($0.5$), $356$ ($5$). So the bound holds and
  is tight at small $\beta$; at $\beta_0(3)$ the excess is about $1+\frac{10}3\beta\approx1.0013$. `sun_group.u2_bound` equals my transcription of the
  lemma to 12 digits, and equals $\beta^2/2000+R/100$ for $\beta\le1/20$ (`compare_group_data.out`, part (c)). For $SU(3)$ the ratio is
  $1.0006$–$125$ for $\beta\in[10^{-3},5]$ (`weyl_sun_N3.out`).
* Where $u_2$ enters (all the uses I found: Lemma 10.6 $n_B$, $(u_2/u)^4$ and $\sigma_6$; $\theta_7$ and $b_2$ through $\sigma_6$; $\sigma_5$ of Lemma 16.8):
  each needs only an upper bound for $\sup_{r\notin\{0,\rho,\bar\rho\}}u_r$, so Lemma 16.11 may replace $e^\beta-1-\beta$ everywhere. `dressed` uses
  `g['u2']` for $n_B$ and $\sigma_6$ (hence for $\theta$ and $b_2$ too). `upper_exact` keeps the crude $u_2$ in $\bar\sigma_5$, which is valid.
* Is the lemma needed? `ablation_su5.out`: with $u_2\le e^\beta-1-\beta$, $SU(5)$, $d=3$: $\log(\bar\Upsilon/\Upsilon)=0.2536\beta^4$ and
  $D=0.5986\beta^4>0.414\beta^4$. The remark's $0.254$ and $0.599$ are confirmed, so in $d=3$ the lemma is needed. In $d=4$ it is not:
  $\log(\bar\Upsilon/\Upsilon)=0.0370\beta^4$ and $D=0.5557\beta^4<1.084\beta^4$. The remark claims the need only for $d=3$, correctly.

**B. The $SU(5)$ group data, by independent methods.**
* Exact Haar moments by alias-free Weyl integration (`moments_check.py 5 19`: grid $28^4$, exact for $a+b\le19$; the largest distance to an
  integer is $5.9\cdot10^{-10}$). All 210 moments $\int\chi^a\bar\chi^b$, $a+b\le19$, equal `sun_group.moment(5,a,b)`. For example
  $E|\chi|^{2k}=1,1,2,6,24,120,719,5003,39429,344837$ and $E\chi^5=1$.
* From these moments, my own exact series (Fractions, written from the definitions of Section 10):
  * $u=\frac\beta{50}+\frac{\beta^4}{1.2\cdot10^6}-\frac{7\beta^6}{3.6\cdot10^9}+\dots$
  * $-\log U_+=-\frac{\beta^3}{6000}+\frac{7\beta^5}{12\cdot10^6}+\frac{\beta^6}{72\cdot10^6}-\dots$
  * $-\log U_-=+\frac{\beta^3}{6000}-\frac{\beta^5}{4\cdot10^6}+\dots$

  So $\Upsilon=U_+$, and $-\log\Upsilon=-\beta^3/6000+O(\beta^5)$ with no $\beta^4$ term, as stated in (i). These coefficients agree with
  `orders/sun_group.out`.
* Values (`compare_group_data.out`):
  * (a) Direct double-precision Weyl integration of the weighted integrals at $\beta=10^{-3},\dots,5$ against `sun_group.sun(5)`: $u$ agrees to
    $10^{-13}$–$10^{-16}$ relative; $U_\pm$ agree to $\le7\cdot10^{-16}$ absolute; $\Upsilon=U_+$ at every point.
  * (b) My degree-17 exact series against `sun(5)` (degree 40) at 60 digits, for $\beta=10^{-12},\dots,10^{-2}$, including $\beta_0(4)$ and
    $\beta_0(3)$: $u$ agrees to $10^{-57}$ or better, and $U_\pm-1$ agree to all printed digits, e.g. $U_+-1=1.08777587274\cdot10^{-14}$ at
    $\beta_0(3)$.
* The truncation claim (line 4937) holds: the tail beyond degree 40 is at most $N^2\sum_{k>40}\beta^k/k!\approx5\cdot10^{-188}$ at $\beta_0(3)$.

**C. The lower bound (authors' evaluators, rerun; `reruns/`).**
* $SU(5)$, $d=3$, `practical_beta.py --group 5 --terms orders/d3_s9.json --closed 10:10,11:11,12:24,13:0,14:70 --closed-from 15 3 β0(3)`:
  * $D/\beta^4=0.3449$, made of the $\theta$-loss $0.0274$, $\log(\bar\Upsilon/\Upsilon)$ $2.6\cdot10^{-13}$, and $\log(1+\delta')$ $0.3175$.
  * Parts of $\delta'$: majorant $0.2793$, shifts and jumps $0.0364$, exits $0.0015$, one-plaquette $0.0004$.
  * The remark's $0.345$, $0.279$, $0.036$, $0.027$ and $3\cdot10^{-13}$ are confirmed.
* $SU(5)$, $d=4$, with `orders/d4_s10x_92.json` and the closed polymers $10..16$, crude from 17: $D/\beta^4=0.5187$. Parts of $\delta'$: exits $0.1801$,
  two four-edge front-pairs $0.1448$, shifts and jumps $0.0898$, majorant $0.0862$, one-plaquette $0.0178$. The remark's $0.18$, $0.145$, $0.09$ and $0.086$
  are confirmed. With `orders/d4_s10x.json` the result is $D/\beta^4=1.481$ (the remark's $1.48$), and the majorant ($1.05$) dominates.
* Scaling: `scaling_check.py --group 5 β0(3) 1` and `--group 5 --d4 orders/d4_s10x_92.json ... β0(4) 1`. All hypotheses hold at all 17
  values of $t$, with bounds $0.345352$ ($d=3$) and $0.51896$ ($d=4$). The output is identical to `orders/su5_d3_scaling.out` and
  `orders/su5_d4_scaling.out`.
* Monotonicity facts: `monotone_check.py --group 5` (and `--d4`) give "all facts hold: True". The output is identical to the committed files. In
  $d=3$, $\bar\Upsilon-a$ falls from $1$ to $0.99999195$ ($=1-\beta_0/50$), and $\bar\Upsilon(\beta_0)-1=1.09\cdot10^{-14}$.
* The modified scaling argument (16.12(ii) against 16.7(vii)) is correct. At $\beta=t\beta_*$ the right sides of the conditions on $F$, $F_6$ are
  $(1+\kappa)tF_*(\bar\Upsilon(1+\delta')-a)\ge(1+\kappa)tF_*(\bar\Upsilon_*-a_*)$, since $\bar\Upsilon-a$ decreases. The left sides are at most
  $t\cdot$(their value at $\beta_*$)$\cdot(1+O(\kappa\beta_*^3))$, and at $\beta_*$ they were $<F_*(\bar\Upsilon_*-a_*+\bar\Upsilon_*\delta'_*)$. So it
  suffices that $\bar\Upsilon_*\delta'_*+O(\kappa\beta_*^3)\le\kappa(\bar\Upsilon_*-a_*)$, i.e. about $10^{-14}\ll4\cdot10^{-4}$.
  * For $\delta'=(1+\kappa)^2t^4\delta'_*$ one needs only $\bar\Upsilon_*/\bar\Upsilon(\beta)\le1+\kappa$. Here $\bar\Upsilon\in[1,1+1.1\cdot10^{-14}]$, given
    $\Upsilon\ge1$ (issue 2).
  * The remaining steps of (vii) ($\omega=O(t^{4-p})$, $c_J$ non-increasing, $\theta/\beta^6$, $\log(\bar\Upsilon/\Upsilon)/\beta^4$ increasing) do not
    involve the group beyond the checked monotonicity facts.
* Regression ($SU(3)$, no `--group`):
  * $D/\beta^4=11.38$ at $\beta_0(4)$ with `d4_s10x.json`, and $2.554$ at $\beta_0(3)$ with `d3_s9.json` (`practical_beta.py`).
  * `scaling_check.py β0(3) 1`: output identical to `orders/d3_scaling.out` ($D_*=2.55433\beta^4$, bound $2.55742$).
  * `scaling_check.py --d4 orders/d4_s10x.json ... β0(4) 1`: output identical to `orders/d4_x_scaling.out` ($D_*=11.3827\beta^4$, the refereed
    $11.383$; bound $11.3879$).
  * `monotone_check.py` and `--d4`: all facts true, numbers identical to `orders/monotone_check.out` and `_d4.out` (only the line labels gained
    "d=3, SU(3):").
  * `upper_exact.py 3`: $E(\beta_0(3))=0.04146\beta^4$; `upper_exact.py 4`: $0.2066\beta^4$. These are the refereed $0.0415$ and $0.21$.
  * `sun_group.py`: output identical to `orders/sun_group.out`.

**D. The set $S^\times_{10}\cup\{(9,2)\}$ and the count $148894$.**
* $E=2$ forces $K=(4,4)$. Fronts have $4$ or $\ge6$ edges, and $|K|\ge2$. If $k\ge6$, then $E\ge(k-4)+2\ge4$. If $k=4$, then $\frac12(\sum K-4)=2$
  forces one further front of 4 edges. `check_set.py` confirms this over all multisets with $|K|\le6$ and sizes $\le24$.
* `check_set.py` rebuilds the set from its definition, without `order_counts.py`/`assemble_set.py`: its own class enumeration, the four cached
  class files (checked consistent where they overlap), the terms and the corners.
  * The 22 classes with a pair in the set are all counted far enough.
  * The 74 terms equal `orders/d4_s10x_92.json` exactly, and the corners are $(2,7),(3,\frac{13}2),(4,6),(5,\frac{11}2),(6,5),(7,4),(8,3),(10,2)$.
  * The same holds for `orders/d4_s10x.json` (73 terms, corner $(9,2)$).
  * The new term is $(297788,9,2,4,1)$: two pairs per junction, $\kappa=e^{c_J}$.
* The README's `assemble_set.py` command reproduces the JSON exactly (equal as parsed JSON). `order_counts.py --check-classes 11`: the pruned class
  generator equals the unpruned one for $m=4..11$ (92 classes at $m=11$, the order `assemble_set.py` uses for this set; `reruns/check_classes_11.out`).
* $148894$, reproduced by `junctions.c 3 9 4,4` (25 s) and by `referee-section16g/junc_partition_big.c` on the fronts of `fronts_enum 3 4` (8 s,
  0 duplicates). Both give $876,42,912,1866,21552,22038,148894$ for $n=3..9$.
* Ablation (`ablation_su5.out`): multiplying the $(9,2)$ count by $0$, $1$ or $10^9$ leaves $D(\beta_0(4))/\beta^4=0.518721$ ($SU(5)$) and $3.93203$
  ($SU(3)$) unchanged to six digits. Even $\times10^{20}$ gives only $0.5525$ and $4.193$. What matters is that $(9,2)$ is no longer a corner, i.e. the
  bookkeeping confirmed above, and not the value of the count.

**E. The upper bound.**
* `upper_exact.py --group 5 3`: $E(\beta_0(3))=0.3178\beta^4$ ($n_0=12$), with $E/\beta^5$ increasing on 401 points. Identical to
  `orders/su5_upper.out`.
* `upper_exact.py --group 5 4`:
  * On the checkpoint (`reruns/su5_upper_d4.out`): $n_0=8,9,10$ give $E(\beta_0(4))=4860$, $81.5$ and $1.587\beta^4$, with
    $T'_Q=6.33\cdot10^{-4}\beta^4$ at $n_0=10$.
  * On the final version, with `COUNTS[4]` extended by $760$, $2616$ (`reruns/su5_upper_d4_final.out`): $n_0=11$ gives $E(\beta_0(4))=0.2511\beta^4$,
    with $E/\beta^5$ increasing on 401 points. Identical to the $d=4$ part of the new `orders/su5_upper.out`.
* `upper_exact.py 4` ($SU(3)$, final version; `reruns/su3_upper_d4_final.out`, 5 min): $E(\beta_0(4))=0.2066\beta^4$ ($n_0=10$) and $0.03273\beta^4$
  ($n_0=11$). This is the $0.033$ of (vi), and the README's "$0.0327$ instead of $0.207$".
* My own evaluator (`upper_su5.py`, written from Lemma 16.8 and Remark 16.9(ii); no code shared):
  * $d=3$: $E(\beta_0(3))=0.317825\beta^4$, equal to the actual $D'_2+\log\Upsilon$ computed with my $u$, $\Upsilon$.
  * $d=4$, $n_0=10$: $E=1.58722\beta^4$. The crude tail of $T'_{S_Q}$ is $5.419\cdot10^{-4}\beta^4$, as stated (line 4954).
  * With $n=11$ entries $N_{11}\in\{0,10^6,10^9,10^{12}\}$ for both $S_Q$ and $S_A$: $E(\beta_0(4))=0.251095\beta^4$, and $E/\beta^5$ increases (401
    points). The threshold is $N_{11}\le7.9\cdot10^{19}$. "$0.2511$ for any counts up to $10^9$" is confirmed.
  * With the proved rooted binomials $N_{11}=\binom{J'+11K}{11}$ (issue 7): $E(\beta_0(4))=0.360293\beta^4$.
  * With the final counts $760$, $2616$ (`upper_final.py`, `upper_final.out`):
    * $SU(5)$: $E(\beta_0(4))=0.251095\beta^4$ against $1.0839\beta^4$.
    * $SU(3)$ ($N^4/\tau_0^2=324$, $\Upsilon\ge1$ proved in Remark 16.9(ii)): $0.206589\beta^4$ with the counts to 10, and $0.0327273\beta^4$ with 11,
      against $15.05\beta^4$.
    * In both cases $E/\beta^5$ increases on 401 points.
* The derivation of $E$ (Remark 16.9(ii)) checked by hand: $V_2\le\Upsilon(1+T'_A)+T'_A(1+3a_1+3a_1^2)$ and $L_2\ge\Upsilon^2(1-B)$, given $\Upsilon\ge1$.
  It depends on the group only through $N^4/\tau_0^2=2500$, $\iota=0$ and $\bar u$.
* The split enumeration (`admissible_split.c` against unsplit; `split/split_check.out`):
  * The split logic is right. Every first-level Redelmeier branch sees the same state as in the unsplit run: all roots are pre-marked, the
    branch of `roots[i]` gets `roots[0..i)` as untried, and marks and counts are restored after each branch.
  * The pieces add up to the unsplit counts in $\mathbb Z^4$ for $S_A$ ($J=2,3,4$) and $S_Q$ ($J=2,3,16$), at $n_{\max}=8$ and $9$. They also
    reproduce the paper's counts ($S_A$: 4 caps, $32,0,0,60$; $S_Q$: 18 caps, $60,0,24,262$).
  * At $n_{\max}=10$ (`split/split_check.out`):
    * $S_A$: unsplit $4,32,0,0,60,506$, equal to the sum of 4 pieces.
    * $S_Q$: the unsplit rerun gives $18,60,0,24,262,1143$, with 27 sets of two components at 10 (1117 s), equal to the sum of its 16 pieces.
    * Every run and piece has `boundary_hit=0`.
  * `orders/admissible_d4_n11.out` sums to $4,32,0,0,60,506,760$ ($S_A$) and $18,60,0,24,262,1143,2616$ ($S_Q$), with two components $27$ at 10 and
    $624$ at 11. Every piece has `boundary_hit=0`.
* **A second count of the 11-plaquette sets for $(S_A,\{p_0\})$, by a different method** (`adm_forced.c`; `forced/`). The program is my own and is
  not Redelmeier's algorithm.
  * Method. A *forced-link search*: from each root, branch at the smallest link outside $S$ covered exactly once over the non-empty sets of the
    plaquettes through it, excluding the rest. Closed states (connected admissible sets) are recorded in a hash table and extended by single
    plaquettes. The sets with two components are formed afterwards as unions of link-disjoint connected ones. The header proves completeness.
  * Validation against every published count:
    * $\mathbb Z^3$ to 12: $S_A$ $2,8,0,0,10,47,60,96$ (16 with two components at 12); $S_Q$ $6,12,0,0,30,93,144,162$ (two components $3,48,42$).
    * $\mathbb Z^4$ to 10: $S_A$ $4,32,0,0,60,506$; $S_Q$ $18,60,0,24,262,1143$ (27 with two components).
  * Result: **$N_{11}=760$ for $(S_A,\{p_0\})$, none with two components** (3.9 s, against about 22 CPU-minutes for the Redelmeier pieces). It
    agrees with `orders/admissible_d4_n11.out`.
  * As instructed, I did not run $(S_Q,\mathcal Q)$ at 11. At $n_{\max}=10$ the program took 3.2 s for $S_Q$, so it could recount $2616$ in seconds if
    the authors want a second count.

**F. Group-specific hypotheses (item 8 of the brief).**
* $u_*=u$: Corollary 6.5 needs $N\le49$ and $\beta e^{2\beta}\le(1+\iota)/N^2=1/25$; true on the range. The remark uses only this part of the
  corollary. The failing condition $N^4\le1200\tau_0^2$ (625 > 300) belongs to the *upper* inequality of Corollaries 6.5 and 5.8 (crude tails, Lemma
  5.5), which the remark does not use. Its upper bound is Theorem 10.2 with the tails of Lemma 16.8, whose only hypothesis is $L_2>0$. That holds,
  since $B<1$.
* $\gamma_g\le\gamma_2$ ($\rho$ without trivial summand), $\epsilon$, $b$, the counts, $\theta$: these do not depend on the group except through $\sigma_6$.
* Lemma 10.6 (`mg:lem:dressed`): its two denominators are positive (asserted in `dressed`).
* Theorem 14.2 (`mg:thm:lower5`) and Lemma 14.1: these need only $u_*=u$.
* Section 15: not used. $\Upsilon=\max(U_+,U_-)=U_+$ directly, and the upper bound takes the eigenvector $\sqrt2\,\mathrm{Re}\,\chi$, with
  $a_1=\sqrt2Nu\le a=\sqrt2N\bar u$.
* The only group-specific fact used and not stated is $\Upsilon\ge1$ (issue 2).
* The abstract (136–139), table rows 5131–5132 and the README rows agree with the evaluations, after the rounding up: $0.345352\to0.346$,
  $0.317825\to0.318$, $0.51896\to0.519$, $0.251095\to0.252$, $3.93384\to3.94$, $0.03273\to0.033$.

## What I could not check

* A second count of the 11-plaquette sets for $(S_Q,\mathcal Q)$ ($2616$). As instructed, I did not start it; my program could do it in seconds
  (see E). By issue 7 its value does not matter for the conclusions: any value up to $7.9\cdot10^{19}$ gives $E<\beta^3/6000$ for $SU(5)$.
* A third, independent method for $148894$. I reran both existing programs instead, and showed by ablation that the value is immaterial (D above).
* The proofs of the theorems that Remark 16.12 applies to $SU(5)$ (Theorems 10.2, 14.2, Lemmas 16.6, 16.8, the scaling of 16.7(vii)) were
  refereed earlier. I checked only that nothing in them is specific to $SU(3)$ beyond $u$, $\Upsilon$, $u_2$, $N$ and $\iota$.

## Files

| file | what |
|---|---|
| `weyl_sun.py`, `weyl_sun_N{3,4,5}.out` | Weyl integration on the maximal torus: $u$, $U_\pm$, the $u_r$ of all irreps with parts $\le4$, the order-2 multiplicities, Lemma 16.11 against the true $u_2$ ($N=3,4,5$) |
| `moments_check.py`, `moments_check_N5.out` | exact moments $a+b\le19$ by alias-free Weyl integration against `sun_group.moment`; exact series of $u$, $U_\pm$, $-\log U_\pm$ |
| `compare_group_data.py`, `.out` | `sun_group.sun(5)` against double-precision Weyl integration and against the exact series at 60 digits; `u2_bound` against the lemma |
| `upper_su5.py`, `.out` | the upper bound's $E(\beta)$ from the formulas of the paper, $d=3,4$; $n=11$ entries, the threshold $7.9\cdot10^{19}$, and the rooted binomial of issue 7 |
| `upper_final.py`, `.out` | $E(\beta_0(4))$ for $SU(5)$ and $SU(3)$ with the final counts $760$, $2616$ |
| `adm_forced.c`, `forced/` | the admissible sets by a forced-link search (second method); outputs for $\mathbb Z^3$ to 12, $\mathbb Z^4$ to 10, and $(S_A,\{p_0\})$ to 11 |
| `check_set.py`, `.out` | $S^\times_{10}\cup\{(9,2)\}$ rebuilt from its definition: classes, terms, corners |
| `ablation_su5.py`, `.out` | the $(9,2)$ count $\times0,1,10^9,10^{20}$; $SU(5)$ with the crude $u_2$; all parts of $\delta'$ |
| `reruns/` | outputs of the authors' scripts as rerun here (the commands are the first line of each file) |
| `split/` | `split_check.sh`, `split10.sh`, `sum_split.py`, `split_check.out` (pieces against unsplit, $n_{\max}=8,9,10$), the unsplit outputs `un_*.txt`, `split_times.txt` |
| `junc/` | `run_junc.sh`, `junctions_3_9_44.out`, `partition_3_9_44.out` (the class $4,4,4$ to nine plaquettes, two existing programs) |

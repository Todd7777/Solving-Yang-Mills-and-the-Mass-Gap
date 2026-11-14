# Referee report (tenth pass on Section 16): Corollary 16.20, the two plaquette glueballs on the whole range

Independent referee checks, made on 14 November 2026, of the items added at the end of Section 16 of `theory/mass-gap-strong-coupling.tex`. I worked on the
working tree as found (`.tex` and `.pdf` both of 04:57; no git metadata in the container). I ran no git command and wrote no file outside this folder.
Bytecode caches were suppressed with `PYTHONDONTWRITEBYTECODE`, and the authors' `run_cap.sh` commands were rerun with their outputs redirected here.

The items refereed:
* Corollary `mg:cor:sectorscap` (Corollary 16.20, lines 5239–5267) and its proof (lines 5269–5305), with the paragraph after it (lines 5307–5308);
* the subsection "The two plaquette glueballs" (lines 5231–5237);
* Lemma `mg:lem:upsbounds`(i) (Lemma 16.13(i), lines 4991–4992, and the added sentence of its proof, lines 5000–5001);
* the new sentence of the abstract (lines 145–148) and the summary-table rows (lines 5445, 5449).

Theorem 16.18 and Lemmas 16.14–16.17 I checked only where the corollary uses or extends them: the sectors, $SU(4)$, and Lemma 16.16(ii) for all $n\le40$.

Nothing here bears on the Yang–Mills problem itself, which is open. The corollary concerns lattice Yang–Mills theory at strong coupling,
$\beta<\beta_0(d)\approx4\cdot10^{-4}$ ($d=3$) and $1.5\cdot10^{-4}$ ($d=4$). Nothing in the refereed text claims or implies more.

## Summary

The corollary is correct, and its constants are right. Each claim was checked by an independent argument or an independent program:

* **The lower bound per sector.** Its mathematics holds: Lemma 15.3 and Lemma 16.13(i) give $e^{-2\beta}\le\Upsilon_C\le\bar\Upsilon_C\le\Upsilon_C+\Delta$ and
  $\Upsilon_C\le e^\beta$ with the same $\Delta$. So the hypotheses in bound mode do not depend on the sector, and Lemmas 16.14, 16.15 carry over verbatim.
* **The $d=3$ set.** The empty exponents follow from Lemma 16.16(ii). I recomputed the corners, $(2,7),(3,6),(5,5),(6,4),(41,2)$, and confirmed the face
  argument by brute force. On the tori with sides $\ge8$ of the finite-volume step, the extra exponents are empty as well.
* **The upper bound per sector.** I derived it myself from Theorem 10.2. It agrees with `sector_upper.py` to 40 digits, and $E_C/\beta^5$ has non-negative
  coefficients.
* **The series and the leading term $2(\beta/2N)^{N-2}/(N-2)!$.** I recomputed both from Haar moments obtained by a different method (Toeplitz
  determinants), exactly, for $N=3,\dots,8$. The bound on $Q$ I recomputed with my own tail estimate, in ball arithmetic.
* **Ball arithmetic.** I adapted the ninth referee's Arb evaluators. They certify:
  * every hypothesis at $\beta_*$ in the four new cases ($d=3$ with `d3_s9z.json` for $SU(3)$, $SU(4)$, $SU(5)$; $SU(4)$ in $d=4$), including the stronger form of
    Lemma 16.15, which I check explicitly;
  * $(1+\kappa)^3D(\beta_*)/\beta_*^4\le0.5993438$, $0.1332402$, $0.0728881$, $1.0821542$, below $c'_-=0.600$, $0.134$, $0.073$, $1.083$;
  * $E_C/\beta^4$ at an exact point above $\beta_0(d)$, in all twelve cases (three groups, two dimensions, two sectors), below the table's $c'_+$;
  * $Q\ge q_Nt^{N-2}(1-10^{-8})$, and the positivity of the lower bound for $m_--m_+$ on the whole range.

  Every constant of the table is thereby certified, and the theorem's four cases are re-certified with the stronger-form checks added.

There is no error that affects a stated result. The issues are minor:
* three intermediate bounds in the proof are rounded down instead of up (Issue 1);
* one general sentence about $SU(N)$ is wrong for $N\ge7$ in $d=4$ (Issue 2);
* the proof leans on the proofs, not the statements, of two earlier theorems (Issue 3);
* the new evaluations were made only in mpmath, without full tail coverage, until the runs here (Issue 4);
* the rest is notation and wording.

## Recommendation

**Minor revisions.** The corollary stands with the stated constants; the "Computer-assisted" label is justified once the ball-arithmetic evaluations of
this folder (or equivalent) are cited for the new cases, as was done for Theorem 16.18.

## Issues

Line numbers refer to `theory/mass-gap-strong-coupling.tex`.

1. **Minor (an intermediate bound stated too low).** Line 5279: "Lemma 16.15 from $\beta_*$ then gives $D\le0.59934\,\beta^4$, $0.13324\,\beta^4$ and
   $0.072888\,\beta^4$". The certified values of $(1+\kappa)^3D(\beta_*)/\beta_*^4$ are higher: $0.59934371$, $0.13324018$, $0.07288809$
   (`outputs/arb_eval_cor_su{3,4,5}_d3_s9z.out`; the authors' own outputs print $0.599344$, $0.13324$, $0.0728881$). So the three inequalities as written are
   false, by $4\cdot10^{-6}$, $2\cdot10^{-7}$ and $9\cdot10^{-8}$.
   * The table constants $0.600$, $0.134$, $0.073$ are unaffected.
   * Write $0.59935$, $0.13325$, $0.072889$, or round up to six digits.
   * `scaling_check.py` prints the bound rounded to nearest ("so D <= 0.13324 beta^4"); it should round up. The SU(4), $d=4$ value $1.0822$ (line 5274) is
     correctly rounded up.

2. **Minor (a general claim that is wrong as stated).** Lines 5236–5237: "The two glueballs of $SU(N)$ split at order $\beta^{N-2}$, through the invariant in
   $\rho^{\otimes N}$". This holds for the plaquette variances, $\log(\Upsilon_+/\Upsilon_-)$, but not for the masses once $N-2\ge4$. In $d=4$ the hopping
   matrices of the two sectors differ at order $\beta^4$: $h_+=34$, $h_-=18$. Corollary 15.7(i) and Theorem 15.5 give, for every $N$,
   $o(\beta^4)\le m_--m_+-\log(\Upsilon_+/\Upsilon_-)\le(h_+-h_-)u^4+o(\beta^4)$.
   * Under Conjecture 15.11 the splitting for $N\ge7$, $d=4$ is $16u^4+o(\beta^4)$, of order $\beta^4$, not $\beta^{N-2}$.
   * For $N=6$, $d=4$ both mechanisms enter at order $\beta^4$.
   * In $d=3$ ($h_+=h_-=8$) the order for $N\ge7$ is open.

   The sentence is also unlabelled (the remark after Lemma 15.6 that it cites is Numerical for $N\le9$). The corollary itself is restricted to $N\le5$, where
   $N-2<4$, and is correct. Restrict the sentence to $N=3,4,5$, and mention the $\beta^4$ hopping split for $N\ge6$, $d\ge4$.

3. **Minor (the proof rests on proofs, not statements).**
   * Lines 5270–5273: "By the proof of Theorem 15.5(i), the lower bound for $m_C$ is that of Theorem 14.2 with $\bar\Upsilon_C$ in place of $\bar\Upsilon$".
     Theorem 15.5(i) is an asymptotic statement ($\liminf$). The fixed-$\beta$ bound that the corollary needs exists only inside its proof, as
     $\sup\spec(P|_{\mathcal H_C\ominus\C\mathbf1})\le u^4\bar\Upsilon_C(1+\delta')e^{8(K+1)\theta}$ whenever the hypotheses of Theorem 14.2 hold with $\bar\Upsilon_C$.
     I checked that the argument is right (§A). But for a computer-assisted Theorem the fixed-$\beta$ statement should be stated: one sentence in Theorem 15.5,
     or a sentence in Lemma 16.14 "and the same for $m_C$ with $\Upsilon_C$". This is the ninth referee's Issue 2 in a milder form.
   * Lines 5284–5286: Theorem 10.2 is stated for "a real eigenvector of $\mathsf V$ for the eigenvalue $\Upsilon$", the largest. The corollary uses it with the
     eigenvector of $\Upsilon_C$, which for $C=-$ is the smaller eigenvalue. The proof of Theorem 10.2 uses only $\mathsf V\psi=\Upsilon_C\psi$ and
     $\langle(1+f_\beta)\psi^2\rangle=\Upsilon_C+a_1^2$ (I checked; Theorem 12.2 already allows either eigenvector). Say so.

4. **Minor (what the evaluation of the new cases covers).** Lines 5279–5282 cite only the authors' mpmath evaluation and `error_budget.py`.
   * As the ninth referee noted for the theorem (their Issue 3), `error_budget.py` bounds only the families with the least room, not the full sum $\Sigma_4$,
     the four-edge part of $\Sigma^{(6)}$, the one-plaquette sums, or $e^{\mathsf A}-1-\mathsf A$.
   * `sector_upper.py` evaluates $E_C$ at $\beta_0(d)$ in mpmath, not at a point above it in balls.
   * For Theorem 16.18 these gaps were closed by the ninth referee's Arb evaluation. For the corollary's new cases they are closed here, as described in the
     second bullet of the Summary.
   * Cite these (or your own ball-arithmetic run), and extend Remark 16.19 ("What is computer-assisted") to the corollary. The new computer steps are: the
     $d=3$ evaluation with the new set; the $SU(4)$ evaluations; the sector upper bounds; and the exact series and the bound on $Q$ (exact rationals with
     majorant tails).

5. **Minor (wording).** Line 5234: "in bound mode only through $e^{-2\beta}\le\Upsilon_C\le\bar\Upsilon_C\le\Upsilon_C+\Delta$". Bound mode also uses
   $\Upsilon_C\le e^\beta$: $\bar\Upsilon'\le e^\beta+\Delta$ enters $\eta$ and the third and fourth hypotheses. The proof (line 5271) has it; add it here.

6. **Minor (notation).**
   * Line 5254: "$c_\pm\beta^3$ with $c_+=\frac1{432}$ and $c_-=\frac1{1296}$" clashes with the constants $c_\pm$ of Theorem 16.18, used two lines later
     in (iii) (line 5264: "$c_-=0.600$") and in the abstract. Rename them, or write the two series out.
   * Lines 5286–5291: $\mathcal T'_p$, $\mathcal T'_Q$ are not defined; Lemma 16.9 and Lemma 16.17 use $\mathcal T'_{S_A}$, $\mathcal T'_{S_Q}$.
   * Line 5292: "$V_2-\Upsilon_C\le\mathcal T'_p(\Upsilon_C+1+4a_1+3a_1^2)$" is correct but weaker than the case analysis of Lemma 16.17's proof, which gives
     $2a_1$ (and $3a_1$ as written there). Harmless: the effect is below $10^{-9}$ relative (`outputs/upper_arb_sectors.out`, "mine, 2a"). Say where it comes
     from, for consistency with Lemma 16.17.

7. **Minor (a sentence that does not apply to $SU(4)$).** Line 5256: "So the third-order term of the mass of each glueball is determined from both sides".
   For $SU(4)$, $-\log\Upsilon_\pm$ is even in $\beta$ and there is no third-order term, so "residual smaller than the term" cannot hold. Restrict the
   sentence to $SU(3)$ and $SU(5)$. For $SU(4)$ what holds is that the second-order terms $\mp\beta^2/128$ are determined with an $O(\beta^4)$ residual.

8. **Minor (what "determined" covers).** Lines 5261–5262 and 5303–5304: "for $SU(3)$ also the term $-\beta^3/648$ is determined". The proof compares
   $|r|\le(c'_-+c'_+)\beta^4$ with $\beta^3/648$ ($0.642$ against $3.83$; $3.984$ against $10.04$); these numbers are right. The remainder of the series,
   $\log(\Upsilon_+/\Upsilon_-)-\beta/3+\beta^3/648=-\beta^4/15552+\dots$, is not bounded. This is the convention of Theorem 16.18 ("the residual is smaller
   than the term"), but here the term is a coefficient of a series. Either say that the convention is meant, or add the one-line bound on the remainder
   (about $6.4\cdot10^{-5}\beta^4$; the exact coefficients are in `outputs/series_check.out`).

9. **Minor (trivial roundings, both in the safe direction).**
   * Line 5303: $(c'_-+c'_+)\beta_0(4)=1.3547\cdot10^{-4}$ for $SU(5)$, quoted as $1.36\cdot10^{-4}$.
   * Line 5308: the factor is $2.459$, quoted as $2.4$.

   `outputs/numbers_check.out` lists all 49 numbers checked. The only other flags are those of Issue 1.

10. **Minor (abstract).** Line 145: "The same bounds hold for each sector of charge conjugation". The constants differ ($c'_+=0.319$ against $0.318$ for
    $SU(5)$, $d=3$; the $d=3$ lower constants are smaller; $SU(4)$ is new). "Bounds of the same form" would be exact.

11. **Minor (summary table).**
    * Line 5445: "$e^{-2\beta}\le\Upsilon,\Upsilon_\pm\le e^\beta$ for every $G,\rho$". $\Upsilon_\pm$ are defined only for $\iota=0$; say so.
    * Line 5437: the planar argument is still attributed to "Rem. 16.7(iv)". It is now Lemma 16.16, whose part (ii) the corollary uses for all $n\le40$.
    * Optionally, row 5449 could list the twelve constants, as the row of Theorem 16.18 does.

12. **Minor (wording in Lemma 16.16(ii), which the corollary now uses for all $n\le40$).** Line 5122: "Each front square, and the domino, is a bounded face of
    $W$" is not literally true. If the middle link of the domino lies in $W$ (always for a seven-edge front), the domino is not a face. If it does not, the
    two squares are not faces. The argument needs only this: a bounded face of $W$ that contains a front square is that square or the domino; any other face
    has four front squares around it, distinct and pairwise without a common link. The conclusion is right: brute force over all placements of three fronts
    of $4,4,4$, $4,4,6$, $4,4,7$ edges in a window, 20480 in all, finds no other bounded face, and at most 3, 4, 4 plaquettes in the union of the bounded
    faces (`corners_check.py` (c)).

13. **Suggestion (finite volume).** Lemma 16.16 is stated in infinite volume, while Theorems 13.4/14.2 run on tori with sides $\ge8$. There the remark after
    the lemma applies: complements of the planar sets also qualify. Those complements have $n\ge L_1L_2-4\ge60>40$ plaquettes, so the exponents added in the
    corollary ($n\le40$) are empty on these tori as well, and the complements lie above the corner $(41,2)$, inside the majorant. The choice $40$ is therefore
    safe in finite volume too; one sentence in the proof would say so. (Outside the scope of this report: the same remark is needed for
    "$Z_{n,4}=0$ for all $n\ge3$" in $d=3$, which `third_order_effective.py` and the ninth referee's `arb_eval.py` use for all $n$. On a torus the complements
    with $n=L_1L_2-2$ add terms of relative size below $10^{-70}$ for $L\ge8$, far inside the room $\ge10^{-9}$, so nothing changes. But the passage to the
    limit should mention them.)

## What I checked, and how

### A. The lower bound per sector (item 1 of the brief)

* **Lemma 15.3 and the formula of Lemma 10.7.** With $\Upsilon_C$ in place of $\Upsilon$:
  * $\upsilon_{1,C}=\Upsilon_C+n_B^2/(\Upsilon_C-(u_2/u)^4e^\beta)$ and $\bar\Upsilon_C=\upsilon_{1,C}+2e^{2\beta}\sigma_6/(\upsilon_{1,C}-u^4(1+\gamma_g))$.
  * The two-by-two block bound of Lemma 10.7's proof applies on the range of $Q_C$. The $\{\rho,\bar\rho\}$ block becomes the line of
    $e^-_\rho+Ce^-_{\bar\rho}$, with exact compression $u^4\Upsilon_C$. The other blocks are restrictions of those bounded there, and $Q_C$ commutes with
    $P^qD_qP^q$ and $S$.
  * Both denominators are positive once $\Upsilon_C\ge e^{-2\beta}>\bar c_2$ and $\upsilon_{1,C}\ge\Upsilon_C>\bar u^4(1+\gamma_2)$, which bound mode asserts.
  * Then $\bar\Upsilon_C-\Upsilon_C\le\bar n_B^2/(e^{-2\beta}-\bar c_2)+2e^{2\beta}\bar\sigma_6/(e^{-2\beta}-\bar u^4(1+\gamma_2))=\Delta$, the same $\Delta$ for both
    sectors.
* **Lemma 16.13(i) for $\Upsilon_\pm$.** Correct: $\sqrt2\,\mathrm{Re}\,\chi_\rho$ and $\sqrt2\,\mathrm{Im}\,\chi_\rho$ have $\int\psi^2\dU=1$ when $\iota=0$
  (both $\int(\mathrm{Re}\,\chi)^2$ and $\int(\mathrm{Im}\,\chi)^2$ are $\frac12$), mean zero, and variances $\Upsilon_\pm$ (Lemma 10.1(i)).
  * Numerically, by Weyl integration (`weyl_check.py`), $e^{-2\beta}\le\Upsilon_\pm\le e^\beta$ for $SU(3)$, $SU(4)$, $SU(5)$ at
    $\beta=0.01,0.3,1,3,10$.
  * For example, $SU(3)$ at $\beta=10$: $\Upsilon_+=0.626$, $\Upsilon_-=0.062$. So $\Upsilon_\pm<1$ does happen, which is why only $e^{-2\beta}$ is used.
* **Where $\bar\Upsilon_C$ enters Theorems 13.4/14.2.** It enters through:
  * the rate $a^4\bar\Upsilon(1+\delta')$ and $w_\varpi$;
  * $C_\varpi$, via $(1+\gamma_g)/\bar\Upsilon$;
  * the right sides of the conditions on $\delta'$, $F$, $F_6$;
  * $\eta$ and the third and fourth hypotheses, via $2\bar\Upsilon(1+\delta')$.

  In bound mode each occurrence is replaced by $e^{-2\beta}$ or $e^\beta+\Delta$, in the safe direction (the ninth referee's table, §C there). Lemma 15.2(ii),
  (iii) supply $P_-Q_C$ and $\nu^{-C}$. In Lemma 10.9 with $Q_C$, a pattern in which the dressed tube is in the sector $(-,-C)$ has another piece and is
  charged $\nu(\mu)\le u^4(1+\gamma_g)$, which $C_\varpi$ covers. So the hypotheses checked in bound mode are literally the same for both sectors, and
  $D=8(K+1)\bar\theta+e^{2\beta}\Delta+\log(1+\delta')$ is the same, since $\log\bar\Upsilon'_C-\log\Upsilon_C\le\Delta/\Upsilon_C\le e^{2\beta}\Delta$.
* **Lemma 16.15 unchanged.** Its hypotheses and monotonicity facts involve only bound-mode quantities ($e^{-2\beta}$, $e^\beta+\Delta$), not $\Upsilon_C$.
* **The passage to $\mu_\beta$.** The proof of Theorem 15.5(i) bounds $\langle\tilde\psi,P^t\tilde\psi\rangle_\nu$ for $C$-eigenfunctions $\psi$, which are dense
  in $\mathcal H_C$ (by $\mathsf Q_C$). Correct; see Issue 3 for the presentation.

### B. The $d=3$ empty exponents (item 2)

* **Lemma 16.16(ii) gives exactly what is used.** With $E=(|p'|-4)+\sum_\psi|\psi|/2$ (Lemma 13.2), and no five-edge fronts:
  * $E=2$ forces three front-pairs of $4,4,4$ edges;
  * $E\in\{3,\frac72\}$ forces $4,4,6$ or $4,4,7$ with $|p'|=4$;
  * four or more front-pairs give $E\ge4$.

  So $(n,2)$, $n\ge4$, and $(n,3)$, $(n,\frac72)$, $n\ge5$, carry no pairs in the sense of Lemma 16.6, and $W_{n,E}=0$ there is legitimate. $S$ stays finite
  ($n\le40$), as Lemma 16.6 requires.
* **The face argument.** I rechecked it by hand. In $\Z^2$ the fronts of 4, 6, 7 edges are the boundaries of a square, of a domino, and of a domino with its
  middle link (minimum degree 2, bipartite). A bounded face other than a front square or the domino has four distinct front squares across its extreme rows and
  columns, pairwise without a common link. The wording point is Issue 12. Brute force confirms it (`corners_check.py` (c)).
* **The corners** (`corners_check.py` (a), my own code, brute force on a window $n\le60$, $E\le30$, with every exponent outside $S$ shown to dominate a corner):
  * $S''_9$: $(2,7),(3,6),(5,5),(6,4),(7,3),(8,2)$;
  * $S''_9\cup Z$: $(2,7),(3,6),(5,5),(6,4),(41,2)$.

  Both agree with `orders/d3_s9.json`, `orders/d3_s9z.json`. The terms of `d3_s9z.json` equal those of `d3_s9.json`, none lies in $Z$, and the only exponents
  with $E\le\frac72$ that carry pairs are $(3,2)$, $(n,3)$ and $(n,\frac72)$ with $n\le4$, as the lemma says.
* **The corner $(8,2)$ set the old majorant.** For $SU(3)$, $(b/\bar b)^n(a/\bar a)^E$ is about $2.7\cdot10^{-11}$ at $(8,2)$ against $1.2\cdot10^{-11}$
  at $(7,3)$. Now $(6,4)$ sets it.
* **The change of the majorant** $M_S$ (exact pairs plus the corner term) at $\beta_*$, in Arb: $SU(3)$ $2.35378\to0.399421$, $SU(5)$
  $0.279548\to0.00700058$, as stated (2.354, 0.399, 0.280, 0.007). Each is at the authors' grid point, $\bar a=2^{-11}$, $2^{-12}$ before and $2^{-10}$ after.
  The exact pairs contribute $0.00262\beta^4$ and $0.00034\beta^4$.
* **Finite volume.** See Issue 13.

### C. The upper bound per sector (item 3)

My derivation from Theorem 10.2 and its proof, with $\psi$ the eigenvector of $\Upsilon_C$, $a_1=|\langle f_\beta\psi\rangle|$ ($=\sqrt2Nu$ for $C=+$, $0$ for
$C=-$, as $P_\rho f_\beta=Nu(\chi_\rho+\chi_{\bar\rho})$) and only $\Upsilon_C\ge e^{-2\beta}$:

* **The variance.** $V_2=(1+\mathcal T)(\Upsilon_C+a_1^2)+\mathcal T-((1-\mathcal T)a_1-\mathcal T)_+^2\le(1+\mathcal T)\Upsilon_C+\mathcal T(1+2a_1+3a_1^2)$.
  * If $(1-\mathcal T)a_1\ge\mathcal T$, use $((1-\mathcal T)a_1-\mathcal T)^2\ge(1-2\mathcal T)a_1^2-2\mathcal Ta_1$.
  * Otherwise $a_1\le2\mathcal T$, so $a_1^2\le2\mathcal Ta_1$.

  Hence $\log(V_2/\Upsilon_C)\le\mathcal T(1+e^{2\beta}(1+2a+3a^2))$. The corollary's $4a$ is weaker and also valid.
* **The covariance.** $\|P_\rho\tilde a\|\ge\Upsilon_C-x$ with $x=\mathcal T(1+a_1)\sqrt2Nu\le\mathcal T(1+a)\sqrt2N\bar u$. Then $(\Upsilon_C-x)_+^2\ge\Upsilon_C^2(1-2e^{2\beta}x)$ and
  $(1-\mathcal T_Q)/(1+\epsilon^6)\ge1-\mathcal T_Q-\epsilon^6$. The last term of $L_2$, divided by $\Upsilon_C^2$, is at most $(1+a'^2)e^{20\beta}N^4\tau_0^{-2}\mathcal T_Q$.
  With $(1-z_1)(1-z_2)\ge1-z_1-z_2$, $L_2/\Upsilon_C^2\ge1-B_C$, the corollary's $B_C$ exactly. So $-\log(L_2/\Upsilon_C^2)\le B_C/(1-B_C)$, and $L_2>0$ as
  $B_C<1$ ($B_C\le8.4\cdot10^{-15}$).
* **Nothing uses maximality.** Theorem 10.2's proof needs only $\mathsf V\psi=\Upsilon_C\psi$ (Issue 3). Theorem 15.5(ii) turns the trial function into a bound
  on $m_C$.
* **`sector_upper.py`** implements the corollary's $E_C$ exactly. My Arb implementation (`upper_arb_sectors.py`, from the ninth referee's `upper_arb.py`
  tails) agrees with it to $3\cdot10^{-40}$ relative in all twelve cases (`sector_upper_compare.py`).
* **$E_C/\beta^5$ is a power series with non-negative coefficients.** Structurally: sums, products and $1/(1-B_C)$ of such series, with $\mathcal T'_S$, $B_C$ of
  order $\beta^5$. Also checked by exact expansion to $\beta^{14}$: coefficients of $\beta^0..\beta^4$ zero, all others $\ge0$, in all twelve cases.
* **Certified in Arb** at the exact point $\beta_0(d)(1+10^{-30})$ above $\beta_0(d)$:

| $E_C/\beta^4$ | $d=3$, $C=+$ / $C=-$ | $d=4$, $C=+$ / $C=-$ | $c'_+$ ($d=3$, $d=4$) |
|---|---|---|---|
| $SU(3)$ | $0.0415283815$ / $0.0415283494$ | $0.0469679411$ / $0.0469679319$ | $0.042$, $0.047$ |
| $SU(4)$ | $0.1305720694$ / $0.1305720449$ | $0.1478297027$ / $0.1478296957$ | $0.131$, $0.148$ |
| $SU(5)$ | $0.3183366296$ / $0.3183366095$ | $0.3605147849$ / $0.3605147792$ | $0.319$, $0.361$ |

  All are below $c'_+$. The closest margin is $SU(3)$, $d=4$: $3.2\cdot10^{-5}$. "The same for both sectors to these digits" is right.

### D. The series and the splitting (item 4)

* **Haar moments by a different method** (`series_check.py`). I used
  $\int_{SU(N)}e^{x\chi+y\bar\chi}\dU=\sum_k\det[c^{(k)}_{i-j}(x,y)]$ (Heine's identity on $U(N)$, summed over powers of $\det$), in exact rationals.
  * All moments with $a+b\le10$ agree with `sun_group.moment` (Schur–Weyl), for $N=3,4,5$.
  * Spot values: $m(4,4)=23,24,24$ and $m(N,0)=1$.
* **The series**, exact to $\beta^{10}$:
  * $SU(3)$: $-\log\Upsilon_\pm=\mp\frac\beta6+\frac{\beta^2}{72}+\{\frac1{432},\frac1{1296}\}\beta^3\pm\frac{\beta^4}{31104}+\dots$, and
    $\log(\Upsilon_+/\Upsilon_-)=\frac\beta3-\frac{\beta^3}{648}-\frac{\beta^4}{15552}+\dots$ (the $\beta^2$ terms cancel).
  * $SU(4)$: $\mp\frac{\beta^2}{128}+\{\frac3{32768},\frac1{98304}\}\beta^4$, and $\frac{\beta^2}{64}-\frac{\beta^4}{12288}+\frac{\beta^6}{2359296}$; even in $\beta$.
  * $SU(5)$: $\mp\frac{\beta^3}{6000}+\{\frac7{12000000},-\frac1{4000000}\}\beta^5$, and $\frac{\beta^3}{3000}-\frac{\beta^5}{1200000}$.

  Identical to `rigorous_bounds.sector_series` to $\beta^6$, and to the Weyl integration at $\beta=0.05,0.2$ (`weyl_check.py` (b)).
* **The leading term.** $\log(\Upsilon_+/\Upsilon_-)=2(\beta/2N)^{N-2}/(N-2)!+\dots$ holds exactly for $N=3,\dots,8$, and $q_N=2^{-(N-2)}/(N-2)!$:
  $\frac12,\frac18,\frac1{48},\frac1{384},\frac1{3840},\frac1{46080}$.
* **$Q$.** The coefficients of $Q=Z_2Z-Z_1^2-Z_{\rm im}Z$ equal the authors' to $t^{10}$.
* **The identities.**
  * $\Upsilon_+=2(Z_2Z-Z_1^2)/Z^2$ and $\Upsilon_-=2Z_{\rm im}/Z$, since $\int\mathrm{Im}\,\chi\,e^{t\mathrm{Re}\chi}=0$ by $U\mapsto\bar U$; so $\Upsilon_+-\Upsilon_-=2Q/Z^2$.
  * $Z_{\rm im}=-(G(2,0)-2G(1,1)+G(0,2))/4$.
* **The majorants for $Z_{\rm im}$.** Its coefficients do change sign: $SU(3)$: $\frac12,-\frac14,\frac18,-\frac1{32},\dots$; $SU(5)$ has negative ones too.
  * The majorant $N^{k+2}/k!$ is valid, as $|\mathrm{Im}\,\chi|,|\mathrm{Re}\,\chi|\le N$.
  * `split_lower` bounds the degree $>K$ part of each product with absolute values, so there is no cancellation.
  * The sum truncated at $K+79$ times $(1+10^{-20})$ covers the exponential tail with an enormous margin.
* **My own rigorous bound** (`split_check_arb.py`, Arb). Exact $q_k$ to $t^{14}$ from my moments; beyond that,
  $|q_k|\le3N^2(2N)^k/k!$ (each of the three products has $|[\cdot]_k|\le N^{k+2}2^k/k!$). On $(0,\beta_0(3)/N]$, $Q\ge q_Nt^{N-2}(1-\delta_N)$ with
  $\delta_N=6.75\cdot10^{-9}$, $4.22\cdot10^{-9}$, $2.84\cdot10^{-9}$. These are the authors' numbers ($3.38\cdot10^{-9}/0.5$, $5.28\cdot10^{-10}/0.125$,
  $5.91\cdot10^{-11}/(1/48)$), and all are $<10^{-8}$.
* **The inequality $\log(\Upsilon_+/\Upsilon_-)\ge2e^{-3\beta}Q$.** It holds: $\log x\ge1-1/x$, $Z\le e^\beta$ ($|\mathrm{Re}\chi|\le N$, $t=\beta/N$),
  $\Upsilon_+\le e^\beta$ (Lemma 16.13(i)), and $Q\ge0$.
* **The positivity of (ii)** on the whole range: after division by $\beta^{N-2}$ the bracket decreases, so it suffices at $\beta_0(d)$. Certified in Arb for all
  six cases. For $SU(5)$: $3.3293\cdot10^{-4}$ against $1.5783\cdot10^{-4}$ ($d=3$, factor $2.109$), and $3.3318\cdot10^{-4}$ against $1.3547\cdot10^{-4}$
  ($d=4$, factor $2.459$; Issue 9).

### E. $SU(4)$ (item 5)

* **Lemma 16.11 for $N=4$.**
  * $\Lambda^2\rho$ is self-conjugate ($\omega_2=\omega_{N-2}$) and occurs in $\rho\otimes\rho$ and in $\bar\rho\otimes\bar\rho$, so its order-2 coefficient is
    $t^2=\beta^2/(4N^2)$, i.e. $c_\Lambda=2$.
  * $\mathrm{Sym}^2$ (dimension 10) occurs once ($\beta^2/(8N^2)$), the adjoint (15) once in $\rho\otimes\bar\rho$ ($\beta^2/(4N^2)$).
  * All other $r$ have $u_r\le R/d_r^2$.
  * Both `sun_group.u2_bound` and the Arb code implement this.
* **The rest of Lemma 16.14.** It holds for $N=4$: $\iota=0$, $u_*=u$ (Corollary 6.4, $N\le49$), $\bar u$, $u_{\rm low}$ (Lemma 5.1(iii)), $\bar\sigma_6$
  (group-free beyond $N$).
* **Lemma 16.15.** $\bar u_2/\beta^2$ is a maximum of $\beta^2$-normalised increasing functions.
* **The evaluations** in $d=3$ (new set) and $d=4$ are certified in Arb (§F).
* **$\Upsilon=\Upsilon_+$ on the whole range.** It follows from (ii): $Q\ge t^2\cdot0.125(1-4.2\cdot10^{-9})>0$ for $t\le\beta_0(3)/4$, which covers $d=4$ too.
  So $m\le m_+\le4\log(1/u)-\log\Upsilon+c'_+\beta^4$. The lower bound $m=\min(m_+,m_-)\ge4\log(1/u)-\log\max\Upsilon_C-c'_-\beta^4$ does not need it.
  Lemma 16.17's hypothesis $\Upsilon\ge1$ is not needed for $SU(4)$, since the sector bound with $\Upsilon_+\ge e^{-2\beta}$ is used instead.

### F. The lower bounds in ball arithmetic (items 6, 7)

`arb_eval_cor.py` is the ninth referee's `arb_eval.py` with the changes marked `COR`:
* the set and its corners read from the JSON of the case;
* the authors' $(\bar b,\bar a)$ grid point;
* $GROUP=4$ allowed;
* the stronger form of Lemma 16.15 checked explicitly: the conditions on $F$ and $F_6$ with $e^{-2\beta_*}-\hat a$, the third hypothesis $\le1/(1+\kappa)$, the
  fourth $\le\varpi_*/(1+\kappa)$, and $\delta'_*\le\kappa$;
* the majorant printed.

Everything else as there: 400-bit balls, fronts to 200 edges and clusters/polymers to 110 plaquettes with rigorous remainders, the referee's marked constants
(six inequalities verified in Arb), $\beta_*\ge\beta_0(d)$ checked in balls. `my_eval.py` is a verbatim copy, used for the parameters.

| case | every hypothesis | stronger form of Lemma 16.15 | $D(\beta_*)/\beta_*^4$ (= authors' to 12 digits) | $(1+\kappa)^3D/\beta_*^4\le$ | $c'_-$ |
|---|---|---|---|---|---|
| $SU(3)$, $d=3$, `d3_s9z` | holds | holds | $0.598620366262$ | $0.5993437096$ | $0.600$ |
| $SU(4)$, $d=3$, `d3_s9z` | holds | holds | $0.133079370401$ | $0.1332401769$ | $0.134$ |
| $SU(5)$, $d=3$, `d3_s9z` | holds | holds | $0.072800124612$ | $0.07288809269$ | $0.073$ |
| $SU(4)$, $d=4$ | holds | holds | $1.081655087292$ | $1.082154131$ | $1.083$ |
| $SU(3)$, $d=4$ (theorem) | holds | holds | $3.935071589926$ | $3.936887115$ | $3.937$ |
| $SU(5)$, $d=4$ (theorem) | holds | holds | $0.519171567638$ | $0.5194110980$ | $0.520$ |
| $SU(3)$, $d=3$, `d3_s9` (theorem) | holds | holds | $2.554552301740$ | $2.557639097$ | $2.558$ |
| $SU(5)$, $d=3$, `d3_s9` (theorem) | holds | holds | $0.345567127983$ | $0.3459846943$ | $0.346$ |

Rooms in the new cases:
* decoration $\ge9.7\cdot10^{-16}$, barred $\ge6.7\cdot10^{-16}$;
* the six inequalities hold (the referee's construction);
* the conditions on $\delta'$, $F$, $F_6$ $\ge1.0\cdot10^{-9}$, in the stronger form too;
* the third and fourth hypotheses $\approx1$;
* $\delta'_*\approx10^{-14}\le\kappa$.

The authors' error budget of the new cases (`orders/error_budget_sectors.out`, reproduced) is right as quoted: tails $\le2.76\cdot10^{-27}$, room
$\ge2.84\cdot10^{-16}$ ($SU(4)$, $d=4$, `jc_bar5`).

### G. Reproductions (item 6)

`run_reproductions.sh` reruns the corollary part of `run_cap.sh` with the outputs in `outputs/repro_*.out`:
* the four scalings, `error_budget.py --d3 orders/d3_s9z.json --groups 3,4,5`, `sector_upper.py` and `rigorous_bounds.py`;
* all are identical to the files in `../orders/` (apart from timing lines);
* about 5 minutes with 2 jobs.

`numbers_check.py` checks all 49 numbers of the corollary, its proof, the paragraph after it and the subsection against these outputs and the Arb outputs. 44
agree; the 5 flags are Issue 1 (three bounds rounded down) and Issue 9 (two safe roundings).

### H. Labels and claims (item 8)

* **"Computer-assisted" and "Theorem (computer-assisted)" are justified.** The corollary rests on:
  * the enumerations of Theorem 16.18, each by two programs;
  * exact rational computations with rigorous majorant tails (the series, $Q$), now confirmed independently;
  * one evaluation at $\beta_*$ per new case and the twelve $E_C$, now certified in ball arithmetic (Issue 4: cite it).
* **No grid or numerical monotonicity is used.**
* **No mathematical gap was found.** Issue 3 asks for two fixed-$\beta$ facts that the proofs contain to be stated.
* **No overclaim.** Nothing claims or implies that the Yang–Mills mass gap problem is solved. The corollary, the subsection, the abstract sentence and the table
  rows all concern lattice theory at $\beta<\beta_0(d)$, and "plaquette glueball" denotes the sector gaps $m_\pm$ of Definition 15.4, as in Section 15.
* **The README status board** ("Corollary 16.20, computer-assisted, not yet refereed") is accurate.
* **One unlabelled claim:** Issue 2.

## What I could not check

* **The enumerations.** The counts inherited from Theorem 16.18 ($S''_9$ terms, exits, $N^v_k$, closed polymers, admissible sets) I did not recount. They were
  checked by earlier referees, and the corollary adds no new count. For the extra exponents it uses a theorem (Lemma 16.16(ii)), which I checked.
* **The earlier sections.** Sections 10–15 (Theorems 10.2, 13.4, 14.2, 15.5, Lemmas 10.7, 15.2, 15.3) I read only as far as the corollary uses them, and I
  did not re-derive them.

## Files

Run from this folder: `bash run_checks.sh` (about 20 s on 4 cores). From `theory/checks/third-order`: `bash referee-section16j/run_reproductions.sh 2`
(about 5 min).

| file | what |
|---|---|
| `run_reproductions.sh` | the corollary part of `run_cap.sh`, outputs to `outputs/repro_*.out` (identical to `../orders/`) |
| `run_checks.sh` | all my checks below |
| `corners_check.py` | §B: corners of $S''_9$ and $S''_9\cup Z$ (own code), the JSON files, brute force of the planar face argument |
| `series_check.py` | §D: Haar moments by Toeplitz determinants (exact), the series, $Q$, the leading term for $N=3..8$ |
| `split_check_arb.py` | §D: the bound $Q\ge q_Nt^{N-2}(1-10^{-8})$ with own tail bound, positivity of (ii), the quoted comparisons, in Arb |
| `weyl_check.py` | §A, §D: Lemma 16.13(i) for $\Upsilon_\pm$ at $\beta$ up to 10, and the series against the Weyl integration (Numerical, a sanity check) |
| `upper_arb_sectors.py` | §C: $E_C$ in Arb above $\beta_0(d)$, all twelve cases; my $2a$ variant; Lemma 16.17's $E$; the power-series structure to $\beta^{14}$ |
| `sector_upper_compare.py` | §C: `sector_upper.E_sector` against `upper_arb_sectors.E_sector`, to 40 digits |
| `arb_eval_cor.py` | §F: the lower bound in Arb, adapted from `referee-section16i/arb_eval.py` (changes marked `COR`) |
| `my_eval.py` | verbatim copy of `referee-section16i/my_eval.py` (parameters, settings, group bounds) |
| `numbers_check.py` | §G: the 49 quoted numbers against the outputs |
| `outputs/` | all outputs, each with its command line |

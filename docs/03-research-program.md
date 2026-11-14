# 3. Research program

This is the plan of attack. It is deliberately conservative about what is proved and deliberately
ambitious about what is attempted. It will be revised at the end of every session
(`journal/`), and milestones move from "open" to "done" only when a complete proof is in `theory/`.

## 3.1 Principles

1. **Everything is labelled.** Every mathematical statement in this repository is one of
   *Theorem* (complete proof in `theory/`), *Cited theorem* (published proof, reference given and verified),
   *Conjecture* (precise statement, no proof), *Heuristic* (physical argument), or *Numerical* (data in
   `results/` produced by code in `lattice/`). An unlabelled claim is a bug.
2. **The $U(1)$ test.** Any proposed proof of the mass gap must visibly fail for $G=U(1)$ in four
   dimensions and for the $O(2)$ model in two dimensions, where massless phases are theorems. If an argument
   does not use a non-abelian feature at an essential step, it is wrong. We apply this test to every idea
   before investing in it.
3. **The strong-coupling test.** Any proposed argument, restricted to $\beta<\beta_0$, must reduce to
   something equivalent to the cluster expansion, with the correct asymptotics $m(\beta)\sim 4\log(1/u(\beta))$.
4. **Numerics as a compass, never as evidence.** Every conjecture is checked against the Monte Carlo code
   before it is pursued. Numerical agreement is necessary, never sufficient.
5. **Reduce before attacking.** The problem is split into independent lemmas with precise hypotheses,
   so that partial progress is cumulative and so that we can measure how far any idea actually goes.

## 3.2 The target in lattice form

From `docs/01-problem-statement.md`, §1.4–1.6: fix $G=SU(N)$ and the Wilson action. Let
$m_L(\beta)$ be the gap of the transfer-matrix Hamiltonian on the spatial lattice $(\mathbb{Z}/L)^3$ and
$m(\beta)=\liminf_{L\to\infty}m_L(\beta)$.

* **(G1)** $m(\beta)>0$ for all $\beta$.
* **(G2)** $\liminf_{\beta\to\infty}m(\beta)/a\Lambda_L(\beta)>0$ with $a\Lambda_L$ the two-loop asymptotic-scaling function.
* **(UV)** Hypothesis (L1): convergence of renormalised gauge-invariant Schwinger functions as $a\to0$.

(G2)+(UV) give the Clay statement via the reduction theorem (monograph §5). (G1) is the "confinement at all
couplings" statement and is logically independent, but any method that proves (G2) will almost certainly
prove (G1).

**Infinite-volume form.** Part II (`theory/infinite-volume-strong-coupling.tex`) shows that for every $\beta$
every infinite-volume limit point $\mu$ of the periodic Wilson measures (with even time extents) is a
reflection-positive, translation-invariant Gibbs measure whose Osterwalder–Schrader reconstruction has a
transfer operator $\hat T_\mu$, so the mass gap $m(\mu)=-\log\lVert\hat T_\mu|_{\Omega^\perp}\rVert$ is defined at
every coupling. The statement that matters for a quantum field theory is

* **(G1$_\infty$)** for $G=SU(N)$, $d=4$ and every $\beta>0$, every such limit point has $m(\mu)>0$.

Its relation to (G1) (spatial volume to infinity at infinite time extent) is not settled; at strong coupling
both hold.

## 3.3 Milestones

### Tier 0 — Foundations (this session; monograph `theory/lattice-foundations.tex`)

| ID | Statement | Status |
|---|---|---|
| T0.1 | Wilson measure, gauge invariance, Elitzur-type remark | Theorem (§2) |
| T0.2 | Reflection positivity for site reflections | Theorem (§3) |
| T0.3 | Transfer matrix: self-adjoint, trace class, strictly positive on the gauge-invariant space, positivity-improving; unique vacuum; $Z=\operatorname{Tr}\hat T^{\,T}$ | Theorem (§4) |
| T0.4 | $m_L(\beta)>0$ for every finite $L$; gap $\Leftrightarrow$ exponential clustering in time; the ground-state Markov chain reformulation $m_L=-\log(1-\mathrm{gap}(P_L))$ | Theorem (§4) |
| T0.5 | Reduction theorem: (L1)+(L2)+(L3) $\Rightarrow$ OS axioms + mass gap | Theorem (§5) |
| T0.6 | Strong coupling: for $\beta<\beta_D(d)$ exponential clustering uniform in volume, hence $\inf_L m_L(\beta)>0$ (Dobrushin uniqueness, fully self-contained proof) | Theorem (§6); threshold improved to $\beta_*(d)=\operatorname{artanh}\frac1{6(d-1)}$ ($\beta_*(4)\approx0.0556$, about $4\times\beta_D$) in Part II, Cor. `iv:cor:finitevol` |
| T0.6′ | Strong-coupling gap with the correct leading asymptotics, $d\ge3$: $4\log(1/q(\beta))\le m\le4\log(1/u(\beta))+D(\beta)$ on the whole range $\beta<\beta_0(d)$, hence $m=4\log(1/\beta)+O(1)$ (principle 3, "the strong-coupling test") | Theorem (Part IV, `theory/mass-gap-strong-coupling.tex`). Lower bound (refereed 2026-11-09): duplicate variables and four plaquettes per time layer; finite volume needs spatial sides $\ge4$. Upper bound (2026-11-10; it replaces the 2026-11-09 version, which needed $\beta<10^{-16}$): a positivity property of the transfer matrix isolates the minimal tube (the four sides of a unit cube), which is then evaluated exactly; $D=O(\beta)$; on the whole range when $\dim\rho\le4$ ($\le5$ if self-conjugate), with $D<0.32$ when $\dim\rho\le3$ (numerically $D<0.05$ for $U(1)$, $SU(2)$, $SU(3)$ in $d=3,4$); finite volume: sides $\ge3$. Refereed 2026-11-10. The asymptotics, indeed the full expansion $-4\log\beta+r(\beta)$ with $r$ analytic, is classical (Schor 1983–84; O'Carroll et al. 1985; as reported by Faria da Veiga–O'Carroll, arXiv:2509.03513; originals not yet read); the contribution here is a self-contained proof with explicit constants and range for the infinite-volume theory. Lower bound sharpened on 2026-11-10: temporal plaquettes at the operator norm of their transfer kernel (window $4\log(eK)+O(\beta)$, refereed), and counting by fronts, giving $m(\mu_\beta)=4\log(1/u(\beta))+O(\beta)$ (Section 7, refereed; revised the same day to be effective on the whole range: for $SU(2)$, $-1.38\le m-4\log(1/u)\le0.0014$ for all $\beta<\beta_0(4)$ in $d=4$ and $-0.63\le\dots\le0.005$ in $d=3$, Numerical evaluation of the explicit constants; the revision refereed the same day, accept with minor revisions, applied). First-order constant by junctions (Section 8, refereed, minor revisions applied): $m\ge4\log(1/u)-2\beta+O(\beta^2)$. First-order term by characters (Section 9, refereed, minor revisions applied): $m=4\log(1/u)-\kappa(\beta)+O(\beta^2)$, $\kappa=\tau_0^{-1}\int(\mathrm{Re}\,\chi_\rho)^2f_\beta$, which is $O(\beta^2)$ for $SU(2)$, $U(1)$, $SU(N)$ with $N\ne3$, and $\beta/6+O(\beta^2)$ for $SU(3)$; for $SU(2)$, $-1.3\cdot10^{-4}\le m-4\log(1/u)\le4.8\cdot10^{-4}$ on the whole range in $d=4$ (Numerical evaluation of the explicit constants). Second-order term (Section 10, 2026-11-11, refereed: minor revisions, applied): $m=4\log(1/u)-\log\Upsilon(\beta)+O(\beta^3)$, $\Upsilon$ the largest variance of a real unit plaquette function in $\mathrm{span}(\chi_\rho,\bar\chi_\rho)$ under the one-plaquette weight; $+\beta^2/8$ for $SU(2)$ and $U(1)$, $-\beta/6+\beta^2/72$ for $SU(3)$; the shifts and jumps of the tube vanish on the swap-odd sector of the duplicated system. Third order (Section 11, 2026-11-11, refereed: minor revisions, applied): nothing new enters, $m=4\log(1/u)-\log\Upsilon(\beta)+O(\beta^4)$; the $\beta^3$ term is that of $-\log\Upsilon$ ($\beta^3/432$ for $SU(3)$, $-\beta^3/6000$ for $SU(5)$, $0$ for $U(1)$, $SU(2)$, $SU(4)$, $SU(N\ge6)$); asymptotic only. Fourth order (Section 12, 2026-11-12; refereed: minor revisions, applied): the plaquette state moves (bulges and hops through six-edge fronts, turns and jumps through unit cubes of a slice); the hopping matrix has top $h_+=26d-70$ ($h_-=10d-22$ for $C$-odd). Theorem: $\limsup\beta^{-4}(m-4\log(1/u)+\log\Upsilon)\le-h_Cu_1^4$ (e.g. $-\frac{17}{128}$ for $SU(2)$, $d=4$). Conjecture: $m=4\log(1/u)-\log\Upsilon-h_Cu^4+O(\beta^5)$, supported by exact $\mathbb Z_2$, $\mathbb Z_3$ spectra and by the classical series of Münster and Seo (as quoted by Langelage–Münster–Philipsen), which agree at order $u^4$ in four of five channels. Junctions of order $\beta^4$ (Section 13, 2026-11-12; two referees: minor revisions, applied): the shifts and jumps, the turns, jumps and U-turns through unit cubes
($Z_{4,4}\le12(d-3)$, proved) and the widenings to six-edge fronts (two narrowings each) enter at order $\beta^4$; with $\varpi=b^{5/2}$ the shifts and jumps drop out, and the
U-turns too. Theorem: the third-order lower bound holds with $\delta'=(\delta'_4+o(1))\beta^4$, $\delta'_4=160(d-3)+8K_Su_1^2$ ($2$ in $d=3$, $166$ in $d=4$ for $SU(2)$,
instead of $3.4\cdot10^8$ and $8.4\cdot10^{10}$): the moves of the hopping matrix at the sup norm. So $\liminf$ and $\limsup$ of $\beta^{-4}(m-4\log(1/u)+\log\Upsilon)$ lie
in $[-166,-\frac{17}{128}]$ for $SU(2)$, $d=4$. Lower bound to fourth order (Section 14, 2026-11-12; refereed: minor revisions, applied): the junctions of the moves are closed surfaces, on which the transfer
operator is diagonal with norm $u$ per plaquette; so the lower bound has $-h_+u^4$ and, with the upper bound, $m=4\log(1/u)-\log\Upsilon-(26d-70)u^4+o(\beta^4)$ for $SU(N)$,
$2\le N\le5$, and for every self-conjugate $\rho$ (Conjecture 12.4 there, with $o(\beta^4)$), and $-8u^4$ for every $\rho$ in $d=3$; Section 15 (12 Nov): $SU(N)$, $N\ge6$ (so every $SU(N)$), and the two sectors of charge conjugation, each to third order and the even one to fourth ($SU(3)$: $m_--m_+=\beta/3+O(\beta^2)$); Section 16 (12 Nov): closed polymers have 6 or $\ge8$ plaquettes ($\theta_7=O(\beta^6)$), the junctions with three front-pairs of lowest order classified and counted exactly; the $SU(3)$ third-order term is resolved for $\beta\le2.6\cdot10^{-6}$ ($d=4$) and $\beta\le2.3\cdot10^{-5}$ ($d=3$) (numerical; two referees), and with the next order (Lemma 16.6, Remark 16.7; counts by computer, $\theta_{10}$) for $\beta\le1.2\cdot10^{-5}$ and $\beta\le9\cdot10^{-5}$, in $d=3$ two orders more for $\beta\le3.2\cdot10^{-4}$, and with order 9 in part and a scaling of the parameters in $\beta$ (Remark 16.7(vii), 13 Nov) on the whole range $\beta<\beta_0(3)$ (numerical), and with the upper bound's tails counted exactly (Lemma 16.8, Remark 16.9) from both sides: $-2.56\beta^4\le m-4\log(1/u)+\log\Upsilon\le0.042\beta^4$ for $SU(3)$, $d=3$ (refereed); in $d=4$ also both sides on the whole range, $-11.4\beta^4\le m-4\log(1/u)+\log\Upsilon\le0.21\beta^4$ (Remark 16.7(ix)–(x), Lemma 16.10: a bound for fronts by decisions, $N^v_k\le0.734\cdot7.355^k$ in $\mathbb Z^3$; searches in C; refereed 14 Nov: minor revisions, no error, applied). The practical-$\beta$ programme for the $SU(3)$ third-order term is complete in $d=3$ and $d=4$ (numerical evaluations); so is that for $SU(5)$ in $d=3$, the only other $SU(N)$ with a third-order term ($-\beta^3/6000$; Lemma 16.11, a sharper $u_2$; Remark 16.12: $-0.346\beta^4\le m-4\log(1/u)+\log\Upsilon\le0.318\beta^4$), and in $d=4$ the lower half ($-0.519\beta^4$), the upper half pending a count (not yet refereed). Open: the $C$-odd fourth-order term ($U(1)$, $d\ge4$; Conjecture 15.11: the hop/turn cancellation is an operator identity, Lemma 15.8, but the grouping of sets along the path is not done); the $u^6$ term; analyticity of $m+4\log u$ (classical) |
| T0.7 | Polymer (closed-surface) representation of the Wilson measure; convergence of the cluster expansion for $\beta<\beta_0$ via Kotecký–Preiss; analyticity; area law | Theorem modulo cited KP theorem (§7); area law for every representation of non-zero $N$-ality: Theorem (Part III, `theory/area-law-strong-coupling.tex`, checked by two referees 2026-11-09); analyticity not done |

### Tier 1 — Consolidation (next 2–4 sessions)

| ID | Statement | Why |
|---|---|---|
| T1.1 | Infinite-volume limit at strong coupling: existence, uniqueness, translation invariance, OS axioms of the lattice theory for $\beta<\beta_0$ | completes the strong-coupling picture; sets the template for (L1). **Done** for $\beta<\beta_*(d)\supset(0,\beta_0)$ in Part II (`theory/infinite-volume-strong-coupling.tex`, 2026-11-09): uniqueness, convergence with rate, all symmetries, clustering, global Markov property, OS reconstruction, unique vacuum, infinite-volume gap $\ge\log(1/(6(d-1)\tanh\beta))$; checked by an independent referee and revised (2026-11-09) |
| T1.2 | Hamiltonian lattice gauge theory: anisotropic transfer matrix, $a_t\to0$ limit, Kogut–Susskind Hamiltonian as a self-adjoint operator; uniform gap at strong coupling in the Hamiltonian picture | makes the Dirichlet-form reformulation (§3.5) rigorous in continuous time |
| T1.3 | Numerical: $a\,m_{0^{++}}(\beta)$ and $a\sqrt\sigma(\beta)$ for $SU(2)$, $\beta\in[2.2,2.6]$, several volumes; test of (G2) in the form $m/(a\Lambda_L)\to$ const; finite-volume behaviour of $m_L$ | gives the quantitative shape of what must be proved and checks the toolkit against published data |
| T1.4 | Numerical: the ground-state Markov chain spectral gap on small spatial lattices ($2^3$, $3^3$) by direct diagonalisation in a truncated character basis, as a function of $\beta$ | first look at the reformulation of §3.5 in a computable setting |
| T1.5 | Re-read and summarise Balaban's program: the exact theorems, the norm on the space of effective actions, and the precise point at which the argument stops | required before any UV work |

### Tier 2 — Model problems (where the crossover must be crossed in a simpler setting)

| ID | Problem | Relation to Yang–Mills |
|---|---|---|
| T2.1 | $O(3)$ non-linear sigma model in $d=2$: mass gap for all $\beta$ (Polyakov's conjecture) | same structure: asymptotic freedom, dimensional transmutation, abelian counterpart ($O(2)$) massless by a theorem. Any method that works here is a candidate for $SU(N)_4$ |
| T2.2 | $SU(2)$ Yang–Mills in $d=3$: existence of the continuum measure in finite volume | super-renormalisable UV; stochastic quantisation and Balaban's 3-d results give partial control |
| T2.3 | $SU(2)$ Yang–Mills in $d=3$: mass gap | Polyakov's monopole mechanism is rigorous for $U(1)_3$ (Göpfert–Mack); the non-abelian version is open |
| T2.4 | 4-d $SU(N)$ lattice theory: existence of the thermodynamic limit (full convergence as $L\to\infty$) for all $\beta$ | strictly weaker than (G1); would already be new |

### Tier 3 — The two hard problems

| ID | Problem |
|---|---|
| T3.1 | (UV) Convergence of renormalised gauge-invariant correlation functions as $a\to0$ in finite volume (upgrade of Balaban/MRS from bounds to limits), with restoration of Euclidean invariance |
| T3.2 | (G2) A non-perturbative lower bound $m(\beta)\ge c\,a\Lambda_L(\beta)$ for large $\beta$ |

## 3.4 The crossover, stated as a mathematical problem

Write the Wilson measure at coupling $\beta$ on lattice spacing $a$ as $\mu_{\beta}$. A block-spin
renormalisation map $\mathcal{R}$ (Balaban's, or any other) sends $\mu_\beta$ to an effective measure
$\mathcal{R}\mu_\beta$ on the lattice of spacing $2a$ for gauge-invariant block observables. Along the
asymptotic-freedom trajectory the effective coupling grows under iteration:
$g^2_{k}\approx g^2/(1-2b_0g^2k\log2)$. After $k_\star\approx1/(2b_0g^2\log2)$ steps the effective coupling is
$O(1)$, and the block lattice has spacing $2^{k_\star}a\approx1/\Lambda$.

* For $k\ll k_\star$ the effective action is controlled by the weak-coupling (small-field) expansion plus a
  large-field estimate (Balaban).
* If, for $k\ge k_\star+k_1$, the effective action were in the convergence region of a strong-coupling
  expansion with uniform constants, (G1) and (G2) would follow (the strong-coupling gap in block units is
  $O(1)$, i.e. $O(\Lambda)$ in physical units).

**The crossover problem.** Find a space of effective actions $\mathcal{A}$ with a norm $\|\cdot\|$ and a
renormalisation map $\mathcal{R}$ such that (i) the weak-coupling expansion controls $\mathcal{R}^k\mu_\beta$ in
$\mathcal{A}$ for $k\le k_\star-k_0$, (ii) a strong-coupling expansion converges uniformly for all actions in a
ball $B\subset\mathcal{A}$, and (iii) finitely many steps of $\mathcal{R}$ map the weak-coupling region into $B$.
Step (iii) is the whole difficulty: at the crossover there is no small parameter, and (iii) is a statement
about finitely many iterations of a non-linear map on an infinite-dimensional space, with no expansion
available. It is exactly the step that the $U(1)$ test says must fail for abelian groups: for $U(1)_4$ the
flow at large $\beta$ goes to the Gaussian (Coulomb) fixed point instead.

Known partial substitutes for (iii):

* **Migdal–Kadanoff bounds** give rigorous one-sided inequalities for free energies under approximate
  decimations, and the approximate flow does go to strong coupling for non-abelian groups in $d=4$ (and for
  $U(1)$ it does not — the approximation respects the $U(1)$ test). The inequalities do not control correlation
  functions; Tomboulis' 2007 attempt to upgrade them is incomplete.
* **Tomboulis–Yaffe inequalities** bound Wilson loops by electric-flux free energies, moving the problem to a
  statement about twisted partition functions (center vortex free energies). The required 't Hooft area law for
  the vortex free energy is open.
* **Large $N$**: at $N=\infty$ the crossover of the Wilson action becomes a bulk first-order transition and both
  sides confine; the loop equations are exact but their solution at weak coupling is not controlled.

## 3.5 An exact reformulation: the gap as a Markov-chain spectral gap

Let $\hat T$ be the transfer matrix on the spatial lattice, $\lambda_0>\lambda_1\ge\dots$ its eigenvalues on the
gauge-invariant space, $\Omega>0$ the vacuum. Define the Markov kernel
$$P(U,U')=\frac{\hat T(U,U')\,\Omega(U')}{\lambda_0\,\Omega(U)} .$$
Then $P$ is a reversible Markov chain on $G^{\text{spatial links}}$ with stationary measure
$\nu=\Omega^2\,dU$ (the vacuum measure), and (monograph, Theorem 4.4)
$$m_L(\beta)=-\log\bigl(1-\operatorname{gap}(P_L)\bigr),\qquad
\operatorname{gap}(P)=\inf_\varphi\frac{\langle\varphi,(1-P)\varphi\rangle_\nu}{\operatorname{Var}_\nu(\varphi)} .$$
The mass gap problem is therefore *exactly* the question whether a reversible, local, gauge-averaged
"parallel update" chain whose stationary measure is the Yang–Mills vacuum measure on a 3-dimensional lattice
has a spectral gap uniform in the volume. In continuous time (Kogut–Susskind Hamiltonian, milestone T1.2) the
same statement reads: the gap equals the Poincaré constant of the Dirichlet form
$\mathcal{E}(\varphi)=\tfrac{g^2}{2}\sum_l\int|\nabla_l\varphi|^2\,d\nu$.

Why this is worth recording:

* For statistical-mechanics spin systems, uniform spectral gaps of Glauber-type dynamics are *equivalent* to
  strong spatial mixing of the stationary measure (Stroock–Zegarlinski, Martinelli–Olivieri). Transplanted
  here, the 4-d mass gap becomes a mixing property of a 3-d measure (the vacuum measure).
* Physically this is the "dimensional reduction" picture of the Yang–Mills vacuum (Greensite; Karabali–Kim–Nair
  in $2+1$): at long distances $\nu\approx$ a 3-d Yang–Mills Gibbs measure. If that picture could be made rigorous
  with uniform constants, it would reduce the 4-d gap to a 3-d statement, i.e. to an induction on dimension.
* It gives a well-posed finite problem (T1.4): compute $\operatorname{gap}(P_L)$ numerically on tiny lattices.

This reformulation passes the $U(1)$ test only trivially (it is an identity); it is a tool, not a mechanism.

**In infinite volume (Part II).** For $\beta<\beta_*(d)$ the unique infinite-volume Gibbs measure has the global
Markov property with respect to time slices, the spatial configurations of successive slices form a stationary
reversible Markov chain with stationary law $\nu$ (the infinite-volume vacuum measure), its transition operator
$P$ is unitarily equivalent to the OS transfer operator, $0\le P\le 1$, and $m(\mu_\beta)=-\log(1-\operatorname{gap}P)$
exactly. So the reformulation is not an artefact of finite volume. Whether the Markov property (and hence the
Markov-chain picture) holds at every $\beta$ is open.

## 3.6 Candidate mechanisms, triaged

| Mechanism | Rigorous status | $U(1)$ test | Assessment |
|---|---|---|---|
| Strong-coupling expansion continued into the weak-coupling region | converges only for $\beta<\beta_0$ | fails (good) | cannot be extended as is; serves as the endpoint of an RG flow |
| Monopole plasma / Debye screening (Polyakov; Göpfert–Mack) | theorem in $U(1)_3$ | passes trivially in $d=3$ (abelian theory is gapped there) | the right template for $d=3$; no non-abelian duality in $d=4$ |
| Center vortices (Mack–Petkova; Tomboulis–Yaffe) | rigorous inequalities only | passes (needs centre) | the best-developed rigorous *structure*; the dynamical statement (vortex condensation) is open |
| Dimensional reduction of the vacuum measure (§3.5) | exact reformulation | identity | promising as an organising principle; requires a 3-d gap result as input |
| Large-$N$ loop equations | strong coupling only | $U(N)$ vs $SU(N)$ differ by a $U(1)$ factor only at $N=\infty$ | not a route to finite $N$ |
| Gribov–Zwanziger horizon, Coulomb-gauge confinement, AdS/CFT | heuristic | — | not rigorous; not pursued |

## 3.7 Anti-patterns (arguments that cannot work)

* "The transfer matrix is positivity improving, hence the vacuum is unique, hence there is a gap." True on every
  finite lattice; says nothing uniform in $L$. ($U(1)$ test fails.)
* "Every lattice theory is gapped; take the limit." The gap may close as $a\to0$ or $L\to\infty$; it does for $U(1)_4$.
* "Use a convergent expansion in $\beta$." No expansion in $\beta$ alone can reach $\beta\to\infty$.
* "Prove the gap in a small box and extrapolate." In a small box the gap is $O(1/L)$ and perturbative (femto-universe); the non-perturbative gap is a large-volume phenomenon.
* "Variational upper bounds." They bound energies from above; the gap needs a lower bound on $E_1-E_0$.
* "Numerical evidence plus a finite-size criterion." No finite-size criterion is known for non-frustration-free transfer matrices; proving one would itself be a major result.

## 3.8 Next-session plan (revised 2026-11-09)

1. Finish the referee check of Part II and revise; then decide whether to merge Parts II and III into the
   monograph or keep them as sequels.
2. T1.3: production runs for $SU(2)$ glueball masses at $\beta=2.3,2.4,2.5$ on $L^3\times2L$, $L=8,10,12$ (the
   container is too slow for the larger lattices; use `lattice/slurm/` on a cluster); compare $a\,m_{0^{++}}$ with
   the two-loop scaling function.
3. T1.2: the Hamiltonian (Kogut–Susskind) version, where the gap is the Poincaré constant of a Dirichlet form on
   the vacuum measure; prove the volume-uniform gap at strong coupling in that setting.
4. T1.4: truncated-character-basis transfer matrix on a $2^3$ spatial lattice; compute $\operatorname{gap}(P)$.
5. T1.5: Balaban summary. T2.1: the $O(3)$ model foundations.
6. Open question raised by Part II: does the global Markov property for time slices hold for every
   infinite-volume limit point at every $\beta$? (It holds at strong coupling, where it follows from uniqueness
   with a rate.) A positive answer would make the Markov-chain reformulation of §3.5 available at all couplings.

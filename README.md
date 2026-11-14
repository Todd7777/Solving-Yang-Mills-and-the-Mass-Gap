# Solving Yang–Mills and the Mass Gap

A long-horizon research repository on the Clay Millennium Prize problem:

> For any compact simple gauge group $G$, prove that a non-trivial quantum Yang–Mills theory exists on
> $\mathbb{R}^4$ and has a mass gap $\Delta>0$.

The problem is open. This repository does not claim to solve it. It is built so that progress
towards a solution is **cumulative, verifiable, and honest**: every mathematical statement is labelled
(Theorem with a complete proof in `theory/`, cited theorem, conjecture, heuristic, or numerical), every
numerical claim is reproducible from `lattice/`, and every session ends with a journal entry recording what
was established, what failed, and what is next.

## Map

| Path | Contents |
|---|---|
| `docs/01-problem-statement.md` | The problem stated precisely: Clay statement, Wightman/Osterwalder–Schrader axioms, definition of the mass gap, reduction to lattice hypotheses (L1)–(L3), the quantitative lattice conjectures (G1), (G2) |
| `docs/02-state-of-the-art.md` | What is rigorously known, numerically known, and believed; known no-go facts; where each technique stops; bibliography |
| `docs/03-research-program.md` | Principles, milestones (Tiers 0–3), the crossover problem in mathematical form, an exact reformulation of the gap as a Markov-chain spectral gap, triage of mechanisms, anti-patterns, next steps |
| `theory/lattice-foundations.tex` | Rigorous monograph: Wilson measure, gauge invariance, reflection positivity, transfer matrix and Hamiltonian, gap ⇔ clustering, ground-state chain, reduction theorem, strong-coupling mass gap (self-contained Dobrushin proof), polymer representation and cluster expansion |
| `theory/infinite-volume-strong-coupling.tex` | Part II: the infinite-volume lattice theory. For every $\beta$, limit points are reflection-positive Gibbs measures with an OS transfer operator, so the infinite-volume gap $m(\mu)$ is defined. For $\beta<\beta_*(d)=\operatorname{artanh}\frac1{6(d-1)}$: unique Gibbs measure, convergence with rate, invariance under all gauge transformations, clustering, global Markov property, transfer operator $=$ a reversible Markov chain of time slices, unique vacuum, mass gap $\ge\log(1/(6(d-1)\tanh\beta))$ |
| `theory/mass-gap-strong-coupling.tex` | Part IV: for $d\ge3$ and all $\beta<\beta_0(d)$ the infinite-volume mass gap satisfies $4\log(1/q(\beta))\le m\le4\log(1/u(\beta))+D(\beta)$, $D=O(\beta)$ ($D<0.32$ when $\dim\rho\le3$; the upper bound needs $\dim\rho\le4$, or $\le5$ if $\rho$ is self-conjugate, on the whole range), so $m=4\log(1/\beta)+O(1)$, the leading behaviour of the strong-coupling expansion. Lower bound: duplicate variables plus a "four plaquettes per time layer" lemma (also finite volume, spatial sides $\ge4$), sharpened by weighing temporal plaquettes at the operator norm of their transfer kernel and by counting surfaces through their fronts, so that $m(\mu_\beta)=4\log(1/u(\beta))+O(\beta)$, with explicit constants on the whole range, grouping by junctions gives $m\ge4\log(1/u)-2\beta+O(\beta^2)$, and evaluating the cap by characters gives the first-order term, $m=4\log(1/u)-\kappa(\beta)+O(\beta^2)$, with $\kappa=O(\beta^2)$ for $SU(2)$, $U(1)$ and $SU(N)$, $N\ne3$, and $\kappa=\beta/6+O(\beta^2)$ for $SU(3)$, and the second-order term, $m=4\log(1/u)-\log\Upsilon(\beta)+O(\beta^3)$ with $\Upsilon$ the largest variance of a
plaquette function under the one-plaquette weight ($+\beta^2/8$ for $SU(2)$) (Sections 6–10; for $SU(2)$, $-1.3\cdot10^{-4}\le m-4\log(1/u)\le4.8\cdot10^{-4}$ for all $\beta<\beta_0(4)$ in $d=4$, a numerical evaluation of the constants). Upper bound: a positivity property of the transfer matrix isolates the minimal tube, which is evaluated exactly by character orthogonality. Checks in `theory/checks/` |
| `theory/area-law-strong-coupling.tex` | Part III: Wilson's area law at strong coupling for every representation of non-zero $N$-ality (and the same for $U(1)$: strong-coupling confinement is not non-abelian) |
| `lattice/` | `ymlattice`: tested SU(2)/SU(3) Wilson-action Monte Carlo (heat bath, overrelaxation, Cabibbo–Marinari), observables (Wilson loops, Polyakov loops, APE smearing, $0^{++}$ correlators), exact formulas (character coefficients, exact 2-d solution, strong/weak-coupling series, asymptotic scaling) |
| `results/` | Validation data and analyses produced by the scripts in `lattice/scripts` |
| `journal/` | One entry per session |

## Standards

1. **Labels.** Theorem / Cited / Conjecture / Heuristic / Numerical. Unlabelled claims are bugs.
2. **The $U(1)$ test.** Any argument for the mass gap must fail for $U(1)$ in four dimensions, where a
   massless phase is a theorem (Guth; Fröhlich–Spencer). Arguments that do not use the non-abelian structure
   at an essential step are discarded.
3. **The strong-coupling test.** Restricted to small $\beta$, any argument must reproduce the cluster-expansion
   results.
4. **Numerics are a compass, not evidence.** Conjectures are checked against the Monte Carlo code before they are
   pursued; agreement is necessary, never sufficient.
5. **No citation without verification.** References not yet re-read in this project are marked (verify).

## Status board

| Item | Status |
|---|---|
| Precise problem statement and lattice reduction | done (`docs/01`) |
| Survey of rigorous and numerical knowledge | done (`docs/02`) |
| Research program with milestones | done (`docs/03`) |
| Foundations monograph (Tier 0 theorems) | see `theory/` and `journal/` for the current state of each section |
| Part II: infinite-volume theory at strong coupling (T1.1) | done for $\beta<\beta_*(d)$ ($\beta_*(4)\approx0.056$); independent referee report addressed (2026-11-09) |
| Part IV: strong-coupling gap with the correct leading asymptotics, $m=4\log(1/\beta)+O(1)$ | done for $d\ge3$ on the whole range $\beta<\beta_0(d)$ (upper bound for $\dim\rho\le4$, or $\le5$ if $\rho$ is self-conjugate; for every $\rho$ when $\beta$ is small); lower bound refereed 2026-11-09, upper bound (new proof) and the operator-norm sharpening (Section 6) refereed 2026-11-10; the sharp constant $m=4\log(1/u)+O(\beta)$ by counting fronts (Section 7) refereed 2026-11-10 (accept with minor revisions, applied) and then made effective on the whole range (closed polymers, exact small fronts, vertex-rooted counting; revision refereed the same day, accept with minor revisions, applied); the asymptotics itself is classical (Schor 1983–84, via Faria da Veiga–O'Carroll 2025), what is new is a self-contained proof with explicit constants and range; first-order constant $\delta=b(1+O(\beta))\approx2\beta$ by junctions (Section 8, 10 November, refereed: minor revisions, applied); the first-order term exactly, $m=4\log(1/u)-\kappa(\beta)+O(\beta^2)$, by evaluating the cap by characters (Section 9, 10 November; refereed: minor revisions, applied): $O(\beta^2)$ for $SU(2)$, $-\beta/6$ for $SU(3)$; the second-order term, $m=4\log(1/u)-\log\Upsilon(\beta)+O(\beta^3)$, $+\beta^2/8$ for $SU(2)$ (Section 10, 11 November; refereed: minor revisions, applied); third order: nothing new, $m=4\log(1/u)-\log\Upsilon(\beta)+O(\beta^4)$ (Section 11, 11 November; refereed: minor revisions, applied); fourth order: the plaquette state moves (bulges, hops, turns and jumps); proved $\limsup\beta^{-4}(m-4\log(1/u)+\log\Upsilon)\le-hu_1^4$, $h=26d-70$, and conjectured equality, consistent with the classical series ($-\frac{98}3u^4$ for $SU(2)$, $d=4$) (Section 12, 12 November; refereed: minor revisions, applied); the junctions of order $\beta^4$ identified ($Z_{4,4}\le12(d-3)$: turns, jumps and U-turns through unit cubes) and the third-order lower bound made effective, $\delta'=(\delta'_4+o(1))\beta^4$ with $\delta'_4=160(d-3)+8K_Su_1^2$ (the moves of the hopping matrix at the sup norm; $166$ for $SU(2)$, $d=4$, instead of $8.4\cdot10^{10}$) (Section 13, 12 November; two referees: minor revisions, applied); the lower bound to fourth order by characters (the moves are closed-surface junctions, diagonal with norm $u$ per plaquette) and **the fourth-order term** $m=4\log(1/u)-\log\Upsilon-(26d-70)u^4+o(\beta^4)$ for $SU(N)$, $2\le N\le5$ ($-\frac{98}3u^4$ for $SU(2)$, $d=4$), and for every self-conjugate representation, and every representation in $d=3$ (Section 14, 12 November; refereed: minor revisions, applied); $SU(N)$, $N\ge6$, via $\Upsilon_+\ge\Upsilon_-$ up to $O(\beta^5)$, so every $SU(N)$, and the two sectors of charge conjugation separately, each to third order and the even one to fourth, with the $SU(3)$ splitting $\beta/3$ (Section 15, 12 November; a first draft claiming the odd sector at fourth order had major gaps found by the referee and was withdrawn; second referee: minor revisions, applied); closed polymers have 6 or $\ge8$ plaquettes, so every lower bound loses only $O(\beta^6)$, and the junctions with three front-pairs of lowest order are classified and counted exactly, which resolves the $SU(3)$ third-order term for $\beta\le2.6\cdot10^{-6}$ ($d=4$) and $\beta\le2.3\cdot10^{-5}$ ($d=3$) (Section 16, 12 November, numerical evaluation; two referees: minor revisions, applied); with the next order counted by computer and $\theta_{10}$, for $\beta\le1.2\cdot10^{-5}$ and $\beta\le9\cdot10^{-5}$, and with two orders more in $d=3$ for $\beta\le3.2\cdot10^{-4}$ (13 November; third referee: minor revisions, applied); with order 9 in part and a scaling of the parameters in $\beta$, in $d=3$ on the whole range $\beta<\beta_0(3)$; with the tails of the upper bound counted exactly as well (Lemma 16.8), for $SU(3)$ in $d=3$ $-2.56\beta^4\le m-4\log(1/u)+\log\Upsilon\le0.042\beta^4$ for all $\beta<\beta_0(3)$, the third-order term from both sides; in $d=4$, $-11.4\beta^4\le m-4\log(1/u)+\log\Upsilon\le0.21\beta^4$ for all $\beta<\beta_0(4)$, the third-order term from both sides on the whole range too (13 November; numerical; the $d=4$ lower bound to $5\cdot10^{-5}$ refereed, minor revisions, applied; beyond, with a bound for fronts by decisions and searches in C, refereed 14 November: minor revisions, no error, applied); open: the fourth-order term in the $C$-odd sector ($U(1)$, $d\ge4$; Conjecture 15.11, with the hop/turn cancellation proved as an operator identity), the $u^6$ term |
| Part III: area law at strong coupling | done; two independent referee reports addressed (2026-11-09) |
| Lattice toolkit, 41 fast tests + 4 slow tests | passing (`cd lattice && pytest` / `pytest --runslow`) |
| Numerical validation against exact 2-d solution, strong/weak coupling, reference values | done (`results/`) |
| $0^{++}$ mass measurements (SU(2)) | in progress (`results/`) |
| (G1) mass gap at all $\beta$ | **open** |
| (G2) gap survives the continuum limit | **open** |
| (UV) convergence of the continuum limit | **open** |

## Running the code

```bash
cd lattice
pip install -e .            # numpy, scipy
pytest                      # fast tests (~1 min): algebra, exact 2-d torus, heat bath, gauge invariance
pytest --runslow            # 4-d checks against strong/weak coupling and the SU(3) beta=6.0 plaquette
python scripts/plaquette_scan.py --group su2 --L 8 --betas 0.5 1 2 2.3 3 6 12 --out ../results/su2.csv
python scripts/glueball.py --group su2 --beta 2.3 --shape 16 8 8 8 --nmeas 2000 --out ../results/gb.npz
python scripts/analyze_glueball.py ../results/gb.npz --md ../results/gb.md
```

Building the monograph requires a TeX distribution: `cd theory && latexmk -pdf lattice-foundations.tex`.

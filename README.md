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
| `theory/mass-gap-strong-coupling.tex` | Part IV: for $d\ge3$ and all $\beta<\beta_0(d)$ the infinite-volume mass gap satisfies $4\log(1/q(\beta))\le m\le4\log(1/u(\beta))+D(\beta)$, $D=O(\beta)$ ($D<0.32$ when $\dim\rho\le3$; the upper bound needs $\dim\rho\le4$, or $\le5$ if $\rho$ is self-conjugate, on the whole range), so $m=4\log(1/\beta)+O(1)$, the leading behaviour of the strong-coupling expansion. Lower bound: duplicate variables plus a "four plaquettes per time layer" lemma (also finite volume, spatial sides $\ge4$), sharpened by weighing temporal plaquettes at the operator norm of their transfer kernel and by counting surfaces through their fronts, so that $m(\mu_\beta)=4\log(1/u(\beta))+O(\beta)$, with explicit constants on the whole range, grouping by junctions gives $m\ge4\log(1/u)-2\beta+O(\beta^2)$, and evaluating the cap by characters gives the first-order term, $m=4\log(1/u)-\kappa(\beta)+O(\beta^2)$, with $\kappa=O(\beta^2)$ for $SU(2)$, $U(1)$ and $SU(N)$, $N\ne3$, and $\kappa=\beta/6+O(\beta^2)$ for $SU(3)$ (Sections 6–9; for $SU(2)$, $-1.3\cdot10^{-4}\le m-4\log(1/u)\le4.8\cdot10^{-4}$ for all $\beta<\beta_0(4)$ in $d=4$, a numerical evaluation of the constants). Upper bound: a positivity property of the transfer matrix isolates the minimal tube, which is evaluated exactly by character orthogonality. Checks in `theory/checks/` |
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
| Part IV: strong-coupling gap with the correct leading asymptotics, $m=4\log(1/\beta)+O(1)$ | done for $d\ge3$ on the whole range $\beta<\beta_0(d)$ (upper bound for $\dim\rho\le4$, or $\le5$ if $\rho$ is self-conjugate; for every $\rho$ when $\beta$ is small); lower bound refereed 2026-11-09, upper bound (new proof) and the operator-norm sharpening (Section 6) refereed 2026-11-10; the sharp constant $m=4\log(1/u)+O(\beta)$ by counting fronts (Section 7) refereed 2026-11-10 (accept with minor revisions, applied) and then made effective on the whole range (closed polymers, exact small fronts, vertex-rooted counting; revision refereed the same day, accept with minor revisions, applied); the asymptotics itself is classical (Schor 1983–84, via Faria da Veiga–O'Carroll 2025), what is new is a self-contained proof with explicit constants and range; first-order constant $\delta=b(1+O(\beta))\approx2\beta$ by junctions (Section 8, 10 November, refereed: minor revisions, applied); the first-order term exactly, $m=4\log(1/u)-\kappa(\beta)+O(\beta^2)$, by evaluating the cap by characters (Section 9, 10 November; refereed: minor revisions, applied): $O(\beta^2)$ for $SU(2)$, $-\beta/6$ for $SU(3)$; open: the second-order coefficient |
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

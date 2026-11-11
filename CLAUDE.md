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

## State on 2026-11-11

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
  minor revisions, applied. Next: the fourth order (hopping of the plaquette state, junctions with four plaquettes), the $C$-odd sector, near $\beta_0$ (journal
  2026-11-11, Part 7).
- Glueball runs: raw `.npz` files and checkpoints in `results/` are gitignored and do not survive the container; only
  the analyses (`results/*.md`) are committed. `glueball.py` resumes from its checkpoint when rerun with the same
  arguments.
- The newest journal entry is `journal/2026-11-11.md` (session 3, Parts 6–7: the second- and third-order terms); `journal/2026-11-10.md` has Parts 1–5,
  `journal/2026-11-09b.md` is the first part of session 3 and `journal/2026-11-09.md` is session 2.

## Where the next steps are

`docs/03-research-program.md` and the newest entry in `journal/`.

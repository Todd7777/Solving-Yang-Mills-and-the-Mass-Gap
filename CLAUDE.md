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

## State on 2026-11-10

- `theory/` holds four monographs: the foundations [F] (`lattice-foundations.tex`, errata 2026-11-09); Part II
  (`infinite-volume-strong-coupling.tex`, the infinite-volume theory, milestone T1.1; refereed and revised 2026-11-09);
  Part III (`area-law-strong-coupling.tex`, the area law; two referees, revised 2026-11-09); Part IV
  (`mass-gap-strong-coupling.tex`, the strong-coupling gap is $4\log(1/\beta)+O(1)$: $4\log(1/q)\le m\le4\log(1/u)+O(\beta)$ for all
  $\beta<\beta_0(d)$ when $\dim\rho\le4$ ($\le5$ if self-conjugate); new on 2026-11-09; the upper bound was reproved by transfer-matrix positivity on 2026-11-10, refereed the same day, with checks in
  `theory/checks/upper-bound/`; the lower bound sharpened the same day to $m=4\log(1/u)+O(\beta)$ by transfer-kernel norms and counting by fronts, Sections 6–7, checks in
  `theory/checks/lower-bound/`; the asymptotics is classical (Schor), see Part IV's Literature paragraph).
- Glueball runs: raw `.npz` files and checkpoints in `results/` are gitignored and do not survive the container; only
  the analyses (`results/*.md`) are committed. `glueball.py` resumes from its checkpoint when rerun with the same
  arguments.
- The newest journal entry is `journal/2026-11-10.md` (session 3, continued); `journal/2026-11-09b.md` is its first part and
  `journal/2026-11-09.md` is session 2.

## Where the next steps are

`docs/03-research-program.md` and the newest entry in `journal/`.

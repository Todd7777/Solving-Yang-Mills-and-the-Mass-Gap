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
cd lattice && pytest -q      # 2026-11-08: 37 pass, 4 slow tests skipped (add --runslow to run them)
```

To build the monographs, first install the TeX packages listed in `.github/workflows/ci.yml` if `latexmk`
is missing, then:

```bash
cd theory && latexmk -pdf -interaction=nonstopmode lattice-foundations.tex
```

The PDFs in `theory/` are committed. Rebuild and commit them whenever their `.tex` source changes.

## State on 2026-11-08

- `theory/area-law-strong-coupling.tex` (Part III of the strong-coupling series) is a draft under review by
  two independent referees. Revisions are pending.
- The last commit (2026-11-07) was made while agents were still running. Nothing later was pushed.
- The newest journal entry is `journal/2026-11-05.md`.

## Where the next steps are

`docs/03-research-program.md` and the newest entry in `journal/`.

#!/usr/bin/env python3
"""Plot a plaquette scan against the strong-coupling (leading and next-to-leading) and one-loop
weak-coupling series.  Usage: python scripts/plot_plaquette.py results/plaquette_su2_8x4.csv out.png"""
import csv
import sys

import numpy as np

from ymlattice import exact


def main():
    src, out = sys.argv[1], sys.argv[2]
    rows = list(csv.DictReader(open(src)))
    N = 3 if "su3" in src else 2
    beta = np.array([float(r["beta"]) for r in rows])
    P = np.array([float(r["plaquette"]) for r in rows])
    err = np.array([float(r["error"]) for r in rows])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bb = np.linspace(0.05, beta.max(), 400)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(bb, [exact.plaquette_strong_leading(b, N) for b in bb], "--", color="#888", label="strong coupling, leading $u(\\beta)$")
    ax.plot(bb, [exact.plaquette_strong_nlo(b, N) for b in bb], "-", color="#c44", label="strong coupling, with cube polymer")
    weak = exact.plaquette_weak(bb, N)
    ax.plot(bb[weak > 0], weak[weak > 0], "-", color="#27c", label="weak coupling, one loop $1-(N^2-1)/(4\\beta)$")
    ax.errorbar(beta, P, yerr=err, fmt="o", color="k", ms=4, label="Monte Carlo, $8^4$ lattice")
    ax.set_xscale("log")
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\beta = 2N/g^2$")
    ax.set_ylabel(r"$\langle P\rangle = \langle \frac{1}{N}\,\mathrm{Re\,Tr}\,U_p\rangle$")
    ax.set_title(f"SU({N}) Wilson action, d=4: plaquette vs coupling")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("wrote", out)


if __name__ == "__main__":
    main()

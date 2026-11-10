# Checks for the sharper lower bound of Part IV (Section 6)

Numerical checks of `theory/mass-gap-strong-coupling.tex`, Section "A sharper lower bound: temporal plaquettes at their character weight"
(Lemmas `mg:lem:layerop`, `mg:lem:transfer`, Theorem `mg:thm:sharp`, Corollary `mg:cor:window`, Remark `mg:rem:sharp`). A compass, not
part of the proofs. Run with `python3 <script>` from this directory (needs `numpy`, `mpmath`).

| Script | What it checks | Runtime |
|---|---|---|
| `s3_transfer_bound.py` | $G=S_3$ (two-dimensional representation) on $(\mathbb Z/2)\times$(one square). The single-copy transfer representation used in the proof of Lemma `mg:lem:transfer` agrees with brute-force summation over $6^9$ gauge-fixed configurations (40 random cases, to $10^{-17}$); $\norm{F_f}=u_*$; the duplicated terms $\Phi_F(X)$ for six sets $X$ satisfy the bound of Lemma `mg:lem:transfer`. | ~5 min |
| `lb_constants.py` | The new lower bound $4\log(1/q_u)$ against the old one $4\log(1/q)$, Part II's $\log(1/(6(d-1)\tanh\beta))$ and the upper bound $4\log(1/u)+D$, for $SU(2)$, $SU(3)$ and $U(1)$ in $d=3,4$, with $u$ computed exactly (Bessel ratios; Weyl integration for $SU(3)$). | ~10 min |

Results on 2026-11-10 are recorded in `journal/2026-11-10.md`.

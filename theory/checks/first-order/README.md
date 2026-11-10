# Checks for the first-order term of Part IV (Section 9, "The cap by characters")

Numerical checks of `theory/mass-gap-strong-coupling.tex`, Section 9 (Lemmas `mg:lem:capelement`, `mg:lem:pieces`, `mg:lem:capchar`,
Theorems `mg:thm:caplower`, `mg:thm:upperchar`, Corollary `mg:cor:firstorder2`, Remark `mg:rem:cap`). A compass, not part of the proofs.
Run with `python3 <script>` from this directory (needs `numpy`, `scipy`, `mpmath`; the constants script imports
`../upper-bound/constants_upper.py` and `../lower-bound/junction_constants.py`).

| Script | What it checks | Runtime |
|---|---|---|
| `cap_operator.py` | Lemma `mg:lem:capchar` by brute force: the norm of the gauge-projected cap operator $P^JO_JP^J$ on the full two-copy space of the four links of a plaquette (matrix-free Lanczos; $\mathbb Z_3$: dimension $3^8$, $S_3$ with its two-dimensional representation: $6^8$) equals the class-function formula $\Lambda M_g\Lambda$ of the proof, lies below $u^4\bar\kappa_{\rm cap}$, and is $u^4(\kappa+(2-\iota)N^2u^2)$ to six digits. For $S_3$, $\int\chi_\rho^3=1$ and $\kappa\approx\beta/2$. | ~15 min |
| `weyl_cubic.py` | Lemma `mg:lem:capelement`(iii) by Weyl integration: $\int(\mathrm{Re}\,\chi)^3=0$ for $SU(2)$, $SU(4)$, $SU(5)$ and $\frac14$ for $SU(3)$ ($\int\chi^3=1$); $\kappa_{SU(3)}(\beta)/(\beta/6)=1+O(\beta)$. | seconds |
| `zn_first_order.py` | Corollary `mg:cor:firstorder2` against exact finite-volume spectra, $d=3$: $\mathbb Z_3$ on the $3\times3$ torus, winding-zero sector ($3^8$ states, $\kappa=u$ exactly): the lowest $C$-even and $C$-odd states have $(m-4\log(1/u))/u\to\mp1$; $\mathbb Z_2$ on the $4\times4$ torus, winding zero ($2^{15}$ states, $\kappa=0$): $(m-4\log(1/u))/u^2\to1$, no first-order term. | ~1 min |
| `z3_4x4.py` | The independent referee's check of the $\mathbb Z_3$ spectrum on the $4\times4$ torus, which satisfies the hypotheses (sides $\ge4$): winding-zero sector ($3^{15}$ states), power iteration with FFTs over $\mathbb Z_3^{16}$ in the sectors of zero momentum, rotation invariance and $C=\pm1$. At $\beta=0.05$ and $0.02$, $(m-4\log(1/u))/u=-0.9382993$, $-0.9751159$ ($C$-even) and $+1.0127402$, $+1.0050505$ ($C$-odd), the $3\times3$ values to six or seven digits. About 2 GB of memory. | ~90 s per $\beta$ |
| `first_order_constants.py` | The constants of Section 9 for $SU(2)$ and $SU(3)$ (fundamental), $d=3,4$, infinite and finite volume: $\kappa$, $\bar\kappa_{\rm cap}$, the smallest admissible $\delta'$ of Theorem `mg:thm:caplower` (with its parts), the loss, and $D^\chi$, $D_0$ of Theorem `mg:thm:upperchar`. $SU(2)$: $\delta_\chi/\beta^2\to128.4$ ($d=4$), $32.4$ ($d=3$), loss $1.28\cdot10^{-4}$ at $\beta_0(4)$, $D^\chi/\beta^2\to0.375$; $SU(3)$: both bounds $\to-\kappa\approx-\beta/6$. With the argument `monotone`: $b_{\rm cap}$ increases on a grid of $(0,\beta_0(d)]$. | ~10 min |

Results on 2026-11-10 are recorded in `journal/2026-11-10.md` (Part 5).

# Checks for the sharper lower bounds of Part IV (Sections 6 and 7)

Numerical checks of `theory/mass-gap-strong-coupling.tex`, Section "A sharper lower bound: temporal plaquettes at their character weight"
(Lemmas `mg:lem:layerop`, `mg:lem:transfer`, Theorem `mg:thm:sharp`, Corollary `mg:cor:window`, Remark `mg:rem:sharp`). A compass, not
part of the proofs. Run with `python3 <script>` from this directory (needs `numpy`, `mpmath`).

| Script | What it checks | Runtime |
|---|---|---|
| `s3_transfer_bound.py` | $G=S_3$ (two-dimensional representation) on $(\mathbb Z/2)\times$(one square). The single-copy transfer representation used in the proof of Lemma `mg:lem:transfer` agrees with brute-force summation over $6^9$ gauge-fixed configurations (40 random cases, to $10^{-17}$); $\norm{F_f}=u_*$; the duplicated terms $\Phi_F(X)$ for six sets $X$ satisfy the bound of Lemma `mg:lem:transfer`. | ~5 min |
| `fronts_check.py` | Lemmas `mg:lem:fronts` and `mg:lem:incidence` (Section 7) on 400 random sets $X$ without free links outside $S$ (unions and mod-2 sums of tubes and cube boundaries; branched configurations included) on the periodic lattice $8\times5\times5$: vertex degrees in fronts, the column lemma, every segment touched at both ends, adjacency in $X$ reflected in the incidence graph, a path from $S_A$ to $S_B$ covering the inner layers, and $|X|=\sum|\sigma|\ell(\sigma)+\sum|C|$. | ~1 min |
| `local_counts.py` | The constants of Lemma `mg:lem:localcounts` by enumeration on spatial tori of dimension 2 and 3 (sides $\ge4$): $D_e=4d-5$, $N_4\le2d-3$, $N_k$ against $e(e(4d-6))^{k-1}$, $2(d-2)$ plaquettes per link, $K_S=4(2d-5)$. | ~2 min |
| `fronts_constants.py` | The constants of Theorem `mg:thm:fronts` ($\lambda_T,\lambda_C$ by iteration, $\delta$) and the loss $8(K+1)\theta+\log(1+\delta)$ for $SU(2)$, $d=3,4$, against Theorem `mg:thm:sharp` (Remark `mg:rem:fronts`). | seconds |
| `z2_4x4.py` | Exact finite-volume gap of $\mathbb Z_2$ gauge theory, $d=3$, on the $4\times4$ spatial torus (physical space $2^{17}$, Lanczos): $m-4\log(1/u)$ is $-10^{-7}$ at $\beta=0.01$ and smaller below, i.e.\ the true gap is $4\log(1/u)$ up to tiny corrections, as the strong-coupling expansion predicts. | ~5 min |
| `lb_constants.py` | The new lower bound $4\log(1/q_u)$ against the old one $4\log(1/q)$, Part II's $\log(1/(6(d-1)\tanh\beta))$ and the upper bound $4\log(1/u)+D$, for $SU(2)$, $SU(3)$ and $U(1)$ in $d=3,4$, with $u$ computed exactly (Bessel ratios; Weyl integration for $SU(3)$). | ~10 min |

Results on 2026-11-10 are recorded in `journal/2026-11-10.md`.

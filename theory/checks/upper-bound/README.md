# Checks for the upper bound of Part IV (Section 5)

Numerical checks of `theory/mass-gap-strong-coupling.tex`, Section "The matching upper bound"
(Lemmas `mg:lem:char`, `mg:lem:pos`, `mg:lem:five`, `mg:lem:twopoint`, `mg:lem:range`, Remark `mg:rem:threshold`).
They are a compass, not part of the proofs. Run each with `python3 <script>` from this directory
(needs `numpy`, `scipy`, `mpmath`).

| Script | What it checks | Runtime |
|---|---|---|
| `constants_upper.py` | The constants of Lemma `mg:lem:twopoint` ($\mathcal T_Q$, $\mathcal T_p$, $\eta$, $L$, $V$, $D=\log(V/L)$) for $SU(2)$, $SU(3)$, $U(1)$ in $d=3,4,5$ across $(0,\beta_0(d))$, and the largest $N^4/\tau_0^2$ for which $L>0$ at $\beta_0(d)$. | seconds |
| `z2_torus.py` | Exact transfer matrix of $\mathbb Z_2$ gauge theory, $d=3$, spatial torus $3\times3$ (physical space of dimension $2^{10}$), infinite time extent: the sixteen terms of the decomposition in the proof of Lemma `mg:lem:pos` are non-negative and sum to the covariance; $\langle H_Q\rangle/(\tau_0u^4)\to1$; for $\beta<\beta_0(3)$ the bounds $\langle H_Q\rangle\ge\tau_0u^4L$, $\mathrm{Var}\le\tau_0V$ and $\log(\mathrm{Var}/\mathrm{Cov})\le4\log(1/u)+D$; the true gap is $\approx3\log(1/u)$ (torelon). | ~1 min |
| `z3_torus.py` | The same for $\mathbb Z_3$ (complex representation, $\iota=0$, $\tau_0=1/2$; physical space $3^{10}$, sparse Lanczos). | ~1 min |
| `s3_ladder.py` | Non-abelian check: $S_3$ with its two-dimensional representation on $(\mathbb Z/L_0)\times$(one square). Positivity of every $\hat T_F$, non-negativity of all terms for $L_0=2,\dots,10$, and agreement of the transfer-matrix terms with brute-force summation over gauge-fixed configurations ($6^9$) at $L_0=2$. | ~2 min |
| `five_plaq.py` | Lemma `mg:lem:five` by exhaustive enumeration of connected plaquette sets: tori $(3,3,3)$, $(4,4,4)$ up to five plaquettes (only the caps occur at five) and $(3,3,3,3)$, $(4,3,3,3)$ up to four, for $B=E(Q)$ and $B=E(p_0)$. | ~1 min |

Results on 2026-11-10 are recorded in `journal/2026-11-10.md`.

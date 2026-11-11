# Checks for the second-order term of Part IV (Section 10)

Numerical checks of `theory/mass-gap-strong-coupling.tex`, Section 10 ("The second-order term"). A compass, not part of the proofs.
Run with `python3 <script>` from this directory (needs `numpy`, `scipy`, `mpmath`).

The quantity: $\Lambda(\beta)$, the largest eigenvalue of $\mathcal T=P_\rho(M_{1+f_\beta}-|f_\beta\rangle\langle f_\beta|)P_\rho$ on the span of $\chi_\rho,\chi_{\bar\rho}$,
i.e. the largest variance, under the single-plaquette weight $w_\beta/c_0$, of a real unit plaquette function $\varphi$ in that span (the $C$-even
$\sqrt2\,\mathrm{Re}\,\chi_\rho$ or the $C$-odd $\sqrt2\,\mathrm{Im}\,\chi_\rho$; $\chi_\rho$ itself if $\rho\cong\bar\rho$). Claim (Section 10):
$m=4\log(1/u)-\log\Lambda(\beta)+O(\beta^3)$.

| Script | What it checks | Runtime |
|---|---|---|
| `lambda_values.py` | $-\log\Lambda_C$ for $U(1)$, $SU(2)$, ..., $SU(5)$ (fundamental) by Weyl integration, and its expansion: $SU(2)$: $+\beta^2/8$; $U(1)$: $+3\beta^2/8$ ($C$-even), $+\beta^2/8$ ($C$-odd, the lighter); $SU(3)$: $\mp\beta/6+\beta^2/72$; $SU(4)$: $\mp\beta^2/128$; $SU(5)$: $O(\beta^3)$. | ~1 min |
| `exact_spectra.py` | Exact transfer-matrix spectra against $4\log(1/u)-\log\Lambda_C$: (1) $\mathbb Z_n$ gauge theory, $d=3$, winding-zero sector, both $C$-sectors ($\mathbb Z_2$ on $4\times4$; $\mathbb Z_3$, $\mathbb Z_4$, $\mathbb Z_5$ on $3\times3$): agreement to $O(u^4)$; (2) a single square, class functions ($SU(2)$, $S_3$, $Q_8$): agreement to all printed digits; (3) two squares sharing an edge, full link space with gauge projection ($S_3$: $6^7$, $Q_8$: $8^7$ states): the symmetric and antisymmetric plaquette states agree with each other and with the formula to $O(\beta^4)$, so the shift of the tube between neighbouring plaquettes does not contribute at second order. | ~30 min |

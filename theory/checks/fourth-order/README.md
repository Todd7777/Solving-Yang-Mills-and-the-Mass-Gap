# Checks for the fourth-order term of Part IV (Section 12)

Numerical checks of `theory/mass-gap-strong-coupling.tex`, Section 12 ("The fourth-order term: the plaquette state moves"). They are a compass, not part
of the proofs. Run each script with `python3 <script>` from this directory; they need `numpy` and `scipy`. `fourth_order_values.py` imports
`../third-order/upsilon_series.py`, which was written by the referee of Section 11.

The claims under test:
- Conjecture 12.4: $m=4\log(1/u)-\log\Upsilon_C-h_Cu^4+O(\beta^5)$ with $h_+=26d-70$ and $h_-=10d-22$, the top of the spectrum of the hopping matrix $H_C$
  (Definition 12.1, Lemma 12.2). For $SU(2)$ in $d=4$ this is $m=-4\log u+2u^2-\frac{98}3u^4+\dots$.
- Theorem 12.8 (proved): the upper bound $m\le4\log(1/u)-\log\Upsilon_C-\log(1+u^4a^\top H_{C,W}a/a^\top a)+O(\beta^5)$ for every finite set $W$ of plaquettes.

The effective model: per unit of time the plaquette state has $K_S$ *bulges*, the cube $c\times[s,s+1]$ attached to a side of the tube. It *hops* to each
plaquette $c$ sharing a link, through one layer of the six-edge front $\partial q\triangle\partial c$. In $d\ge4$ it *turns* or *jumps* through the four
other faces of a unit cube of a time slice. Each process has relative weight $u^4$, and in the $C$-odd sector a hop or turn carries the orientation sign.

| Script | What it checks | Runtime |
|---|---|---|
| `block_spectra.py` | Exact transfer matrix of $\mathbb Z_2$ gauge theory in dual plaquette variables. The state space is the configurations reachable by link flips, i.e. the Bianchi constraints and, on tori, the global ones. Geometries: spatial tori $3\times3$ and $4\times4$ ($d=3$), and blocks of $1\times1\times1$, $1\times1\times2$, $1\times2\times2$, $1\times2\times3$ unit cubes with open boundary ($d=4$). Every level $m_i$ of the plaquette band gives $(m_i-4\log(1/u)+\log(1-u^2))/u^4$, which tends to $-\log(1+h_iu^4)/u^4$ with $h_i$ the eigenvalues of $H$ for the same geometry; the corrections are $O(u^2)$. On $4\times4$ there are 16 levels $8,6^4,4^6,2^4,0$. For one cube the top value is $13$ ($8$ without turns and jumps); for $1\times2\times2$ all 20 levels agree; for $1\times2\times3$ ($2^{23}$ states) the top 4. | 1 s to 12 min |
| `zn_block_spectra.py` | The same for $\mathbb Z_n$, $n\ge3$ (run `python3 zn_block_spectra.py 3 1x1x1 1x1x2 3x3t`), in both $C$-sectors, with the exact $\Upsilon_\pm$ and the signed matrix $H_-$. On one cube the $C$-odd levels follow $5,5,5,3,3,3$ and the $C$-even ones $13,3,3,3,1,1$; on $1\times1\times2$ all 11 levels of each sector agree. | ~1 min |
| `hopping_matrix.py` | Lemma 12.2, building $H_\pm$ from the orientations on tori ($d=3,4,5$): coplanar hops and jumps have sign $+1$, and turns cancel the perpendicular hops for $C=-$; all row sums are $h_+=26d-70$ and $h_-=10d-22$, which are the top eigenvalues. The $k=0$ blocks are $34,10,10$ ($C+$) and $18,18,18$ ($C-$) for $d=4$. Uniform trial vectors on growing boxes approach $h$ (e.g. $31.5$ and top eigenvalue $33.2$ of $34$ at side $8$). | ~2 min |
| `slice_lemmas.py` | Exhaustive enumeration in spatial dimensions 2, 3, 4. (A) Lemma 12.5: a set of plaquettes with no link in exactly one of them has at least 6 plaquettes, and exactly 6 only for a unit cube; none of 7 in dimension 3. (B) Lemma 12.6(i): for a four-edge front, the configurations with $\le4$ further plaquettes are only $q'=q$ with nothing, and the $5\cdot2(n-2)$ cube partners with their four faces. (C) Lemma 12.6(ii): for a six-edge front (coplanar or perpendicular), the configurations with $\le2$ plaquettes are only the bulge and hop completions. The search is complete; see the docstring. | ~1 min |
| `fourth_order_values.py` | The predicted series of the gap, in exact rational arithmetic, for $U(1)$, $SU(2)$–$SU(5)$, $d=3,4$, in $\beta$ and in $u$. For $SU(2)$, $d=4$: $m+4\log u=2u^2-\frac{98}3u^4$, i.e. $m-4\log(4/\beta)=\frac7{24}\beta^2-\frac{167}{1152}\beta^4$. For $SU(3)$: $-3u+9u^2-\frac{27}2u^3-7u^4$ ($C+$) and $3u+\frac92u^3-\frac{99}4u^4$ ($C-$, the $T_1^{+-}$). For $U(1)$ ($C-$): $\frac12u^2-\frac{425}{24}u^4$. | seconds |

Written by the independent referee of Section 12 (added with credit; separate code from the scripts above):

| Script | What it checks | Runtime |
|---|---|---|
| `zn_links_tm.py`, `zn_links_mp.py` | Exact $\mathbb Z_n$ transfer matrices in **link** variables (gauge orbits by label propagation; centre shifts on tori). For each level, the plaquette band against $-\log(1+h_iu^4)$. The equal-time and unit-time covariance matrices of the plaquette functions are checked entry by entry against Steps 2 and 5 of the proof of Theorem 12.8: $G=\Upsilon\delta+u^4\Upsilon^2C$ and $T^{(1)}=u^4\Upsilon^2(\delta+u^4(K_S\delta+A+2C))$, all signs included. `zn_links_mp.py` runs at 50 digits for $\mathbb Z_3$ on one cube. | minutes |
| `geometry_H.py` | Orientations and the matrices $H_\pm$ from scratch. $\epsilon_\Gamma$ is solved from $\sum_p\epsilon(p)\partial p=0$, without outward normals. Covers all sign identities, row sums and $k=0$ blocks. | seconds |
| `enum_lemmas.py` | Lemmas 12.5 and 12.6 by a different exhaustive algorithm: growth of rooted connected sets with a pruning bound. | ~1 min |
| `series_check.py` | Exact series of $u$, $\Upsilon_\pm$ and the predicted masses from Weyl-integrated moments, for $U(1)$, $SU(2)$, $SU(3)$. Also the identification of the tube factors $a,b,c$ of LMP07 eq. (2.9) with $\Upsilon$, $\Upsilon_-$, $\Upsilon_+$. | seconds |
| `lmp_fractions.py` | The printed coefficients of LMP07 eqs. (2.14)–(2.18): which are not in lowest terms, the integer residuals $-34,-10,-34,-10$ against $-71/4$ for $T_1^{+-}$, and the sensitivity to errors in $u_6$, $u_8$. | seconds |
| `moves_by_characters.py` | Lemma 14.1 (Section 14), by exact group sums on $\mathbb Z_3$ and on $S_3$ (non-abelian, two-dimensional representation), $\beta=0.3$ and $1$, for a widening, a turn and a jump. The single-copy integrals $I(a,b,S)$ are $1$, $u_a^{|P|}$ or $0$ exactly as stated, for all irreducible $a,b$ and all subsets $S$; and the two-copy operator $P^JO_JP^J$ on the class functions of the front holonomies has norm exactly $u_*^6$, resp. $u_*^8$ (to machine precision). The domino holonomy is built by hand here; the referee's checks in `referee-section14/` build it from the lattice geometry and compute the norm from the definitions without the reduction. | ~2 min |

## Comparison with the classical series (Remark 12.10(iii))

Langelage, Münster and Philipsen, PoS(LATTICE 2007) 201 (arXiv:0710.0512), eqs. (2.14)–(2.18), quote for $d=4$ the strong-coupling series of Münster
(Nucl. Phys. B190 (1981) 439) and Seo (Nucl. Phys. B209 (1982) 200). The original papers were not consulted. The model gives, at momentum $0$, $h=34$ ($A_1^{++}$),
$10$ ($E^{++}$) and $18$ ($T_1^{+-}$):

| Channel | quoted $u^4$ coefficient | predicted |
|---|---|---|
| $SU(2)$ $A_1^{++}$ | $-98/3$ | $-98/3$ |
| $SU(2)$ $E^{++}$ | $-26/3$ | $-26/3$ |
| $SU(3)$ $A_1^{++}$ | $-7$ | $-7$ |
| $SU(3)$ $E^{++}$ | $+17$ | $+17$ |
| $SU(3)$ $T_1^{+-}$ | $-98/4$ (both arXiv versions) | $-99/4$ |

The terms of orders $u$, $u^2$, $u^3$ agree in all five channels. The $SU(3)$ series of $\Upsilon_\pm$ were checked by Weyl integration to about $10^{-16}$.
The $T_1^{+-}$ difference of $u^4/4$ is very probably a misprint in the printed series (Remark 12.10(iii)): it is the only printed fraction not in lowest
terms, and in the band picture the difference from $-\log\Upsilon_C$ must be an integer on $T_1$.

The independent referee of Section 14 checked Lemma 14.1 and the values of Corollary 14.3 with its own code; see `referee-section14/README.md`.

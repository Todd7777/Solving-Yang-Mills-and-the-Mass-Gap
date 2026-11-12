# Independent checks of Section 14 of Part IV (the referee's code)

Written by the independent referee of Section 14 of `theory/mass-gap-strong-coupling.tex` ("The fourth-order term: the lower bound"),
12 November 2026, and added with credit. A compass, not part of the proofs.

| Script | What it checks |
|---|---|
| `groups.py` | Character tables of $Q_8$, $A_4$, $F_{21}=\mathbb Z_7\rtimes\mathbb Z_3$ (three-dimensional irreducible representations that are not self-conjugate), $S_3$, $\mathbb Z_n$, computed from the class algebra; orthonormality checked. |
| `check_moves.py` | Lemma 14.1 on these groups, built from the lattice geometry: holonomies along oriented cycles, so the domino holonomy is not assumed. Configurations: flat and bent dominos (widening and narrowing), cube turns on two different adjacent faces, and the jump; both front orientations; $\beta=0.4$, $1.3$. For all irreducible $a,b$ and all subsets $S$, every $I(a,b,S)$ is as stated, and $\sigma$ is a bijection with $\sigma(a)\in\{a,\bar a\}$. The two-copy matrix is built from all $S_1\cup S_2=P$; its norm is $u_*^6$, resp. $u_*^8$. For the dominos a direct sum over both copies agrees. |
| `full_operator_zn2.py` | $\lVert P^JO_JP^J\rVert$ computed from the definitions of Lemma 9.2, without the reduction to class functions: $G_g^{1/2}$ as a numerical matrix square root, $M_J=\prod g_p$, and $P^J$ as the average over all gauge transformations in each copy. It also checks that $P^J$ commutes with $O_J$. $\mathbb Z_2$ and $\mathbb Z_3$ dominos, and $\mathbb Z_2$ cube turn and jump (a $2^{24}$-dimensional space): $u_*^6$ and $u_*^8$ to $10^{-12}$ (`outputs/full_operator_zn2.out`). |
| `series_values.py` | The values of Corollary 14.3 and Conjecture 12.4: $-\log\Upsilon_C$ as a series in $u$, from moments by Weyl integration. $SU(2)$: $-\frac{20}3u^4$ ($d=3$), $-\frac{98}3u^4$ ($d=4$), $\frac7{24}\beta^2-\frac{167}{1152}\beta^4$. $SU(3)$: $+19u^4$ ($d=3$), $-7u^4$ ($d=4$). $U(1)$: $-\frac{185}{24}$, $-\frac{425}{24}$. Needs `sympy`. |

Run from this directory (`check_moves.py` imports `groups.py`).

The referee noted that our own check, `../moves_by_characters.py`, takes the hexagon holonomy of the domino as $W_qW_c^{-1}$ by hand, which is the claim
under test, and takes the norm from the class-function formula of the proof. The scripts here avoid both.

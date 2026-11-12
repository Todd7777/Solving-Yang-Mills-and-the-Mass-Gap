# Referee checks for Part IV, Section 16 ("Closed polymers and the third-order term at practical β")

Independent checks by the referee of Section 16 of `theory/mass-gap-strong-coupling.tex` (Lemmas `mg:lem:seven`, `mg:lem:theta7`, Remark
`mg:rem:practical`). The enumerations were written from the definitions in the paper without the authors' scripts (`closed_sets.py`,
`fronts_count.py`, `junction_count.py`), and were compared with them only afterwards. A compass, not part of the proofs.

Run from this directory with `python3` (needs `mpmath`; `ref_numbers.py` imports the authors' constant modules from `../../first-order`,
`../../second-order`, `../../lower-bound`, `../../upper-bound` and `..`). Outputs of the runs are in `outputs/`.

| Script | What it checks | Result | Runtime |
|---|---|---|---|
| `ref_closed.py D NMAX [L]` | Lemma `mg:lem:seven` (seven plaquettes). All link-connected sets of at most `NMAX` plaquettes through a fixed plaquette of $\mathbb Z^D$, or of the torus $(\mathbb Z/L)^D$, in which every link lies in at least two of them. Decision search that also extends closed sets, so closed sets containing a closed proper subset are found too (complete without the caveat of `closed_sets.py`). | $\mathbb Z^3$, $\mathbb Z^4$ up to 10, $\mathbb Z^5$ up to 9: only the $2(D-2)$ unit cubes up to 9 plaquettes; size 10: $10$ ($\mathbb Z^3$) and $60$ ($\mathbb Z^4$) sets. $\mathbb Z^3$ up to 12: $\{6{:}\,2,\ 10{:}\,10,\ 11{:}\,11,\ 12{:}\,24\}$, as counted by hand (box boundaries, two cubes sharing a face, two cubes sharing an edge). Tori $(\mathbb Z/L)^3$, $L=3,\dots,8$ and $(\mathbb Z/L)^4$, $L=3,4,5,8$, up to 9: only cubes, except the flat $3\times3$ torus (9 plaquettes) when $L=3$. | seconds |
| `ref_closed_extra.py` | (1) a plaquette that is not a face of a unit cube contains at most one of its links ($\mathbb Z^3$–$\mathbb Z^5$), the step "two links of $E(\Gamma)$ in a common plaquette span a face" of the case $t=1$; (2) the 60 ten-plaquette sets of $\mathbb Z^4$ are the 20 straight and 40 bent two-cube boxes. | both confirmed | 1 s |
| `ref_fronts.py n KMAX` | Remark `mg:rem:practical`(ii): $N^v_k$, fronts (connected edge sets, all degrees $\ge2$) with $k$ edges through a vertex of $\mathbb Z^n$; decision search (not Redelmeier). | $\mathbb Z^2$: $4,12,12,70,160,634$ ($k=4,6,\dots,10$); $\mathbb Z^3$: $12,132,108,2034,4344$ ($k=4,6,7,8,9$). Cycles: $4,12,56,280$ and $12,132,1656$, i.e. $p_{2n}\cdot2n$ with the known polygon counts $1,2,7,28$ (square) and $3,22,207$ (cubic). | 1 s; 27 s ($\mathbb Z^3$, $k\le9$) |
| `ref_junctions.py m NMAX [list] [nosym] [nofrontadj]` | Remark `mg:rem:practical`(ii): $Z_{n,4}$ from Definitions `mg:def:fronts`, `mg:def:junction` and Lemma `mg:lem:junction` (local admissibility $n_J(e)\ge2$, same-layer fronts vertex-disjoint, both fronts touched, $H_s$ connected), both layers of $p'$. Decision search over $P$ for every $q'$ up to the 16 (resp. 8) lattice symmetries fixing $q_0$; `nosym` searches every $q'$. | $d=3$: $1,8,0,0,0,0,0$; $d=4$: $1,32,0,12,26,302,288$ ($n=1,\dots,7$); the same without symmetry reduction for $n\le5$. The 26 junctions with five plaquettes: $12$ with $q$ added to $\Gamma\setminus\{q,q'\}$, $12$ with $q'$ added, $2$ tubes through the five other faces of a cube. | 1.5 min ($d=4$, $n\le7$) |
| `ref_numbers.py A` | Remark `mg:rem:practical`(i): the losses $8(K+1)\theta_5$ and $8(K+1)\theta_7$; $\theta_7$ increasing in $\beta$ ($SU(2)$). | see below | 2 s |
| `ref_numbers.py B` | Remark `mg:rem:practical`(i): $SU(2)$ lower bound of Theorem `mg:thm:caplower` at $\beta_0(d)$ with $\theta_7$ (copy of `lower_chi` with $\theta$ as a parameter); monotonicity on a grid. | see below | 1 min |
| `ref_numbers.py C`, `E`, `D` | Remark `mg:rem:practical`(iii): $\delta'$ of Theorem `mg:thm:lower5` (`chars`) for $SU(3)$ with $\theta_7$ and the exact counts (copy of `lower3e` with $\theta$, $N^v_k$ and $Z_{n,4}$ as parameters); $d=4$ (C), $d=3$ (E); threshold against $\beta^3/432$ (D). | see below | ~10 min per $\beta$ |

## Results of `ref_numbers.py` (outputs in `outputs/numbers_*.txt`)

**(A) losses.** $8(K+1)\theta_5$ / $8(K+1)\theta_7$: $3.415\cdot10^{-5}$ / $3.207\cdot10^{-9}$ at $\beta_0(4)$; $2.082\cdot10^{-4}$ / $7.589\cdot10^{-8}$ at $\beta_0(3)$;
$3.80\cdot10^{-11}$ / $9.8\cdot10^{-19}$ at $\beta=10^{-5}$, $d=4$, all as in Remark `mg:rem:practical`(i). At these $\beta$ the term $x^8/(1-x)$ dominates $\theta_7$ (the cube
term is $10^{-12}$ smaller at $\beta_0(4)$; for $SU(2)$, $d=4$ the cube term wins only below $\beta\approx1.6\cdot10^{-10}$). $8(K+1)\theta_5=3.79\cdot10^7\beta^4$ at $\beta=10^{-7}$
($d=4$), against $\beta^3/432=2.31\cdot10^4\beta^4$ (ratio $1.6\cdot10^3$). $\theta_7$ ($SU(2)$, $\sigma_6$ exact or by Lemma `mg:lem:capelement`(iv)) increases on a 60-point grid of
$(0,\beta_0(d)]$.

**(B) $SU(2)$, Theorem `mg:thm:caplower`, infinite volume.** Loss at $\beta_0(d)$: $9.39\cdot10^{-5}$ ($d=4$) and $3.56\cdot10^{-5}$ ($d=3$) with $\theta_7$, against $1.280\cdot10^{-4}$ and
$2.438\cdot10^{-4}$ with $\theta_5$ (Remark `mg:rem:practical`(i): $9.4\cdot10^{-5}$, $3.6\cdot10^{-5}$, $1.28\cdot10^{-4}$, $2.44\cdot10^{-4}$). With $\theta_7$ the loss and $b_{\rm cap}$ increase on a
20-point grid of $(0,\beta_0(d)]$, so the whole-range statement follows as in Remark `mg:rem:cap`(i).

**(C) $SU(3)$, $d=4$, Theorem `mg:thm:lower5` (`chars`), $\delta'/\beta^4$** (best of $\varpi=b^2,b^{5/2},b^3$; infinite volume).

| $\beta$ | $\theta_5$, general counts | $\theta_7$, general counts | $\theta_7$, $N^v_{k\le7}$, $Z_{n\le7,4}$ exact | also $N^v_8=2034$ | also $N^v_8$, $N^v_9=4344$ |
|---|---|---|---|---|---|
| $10^{-5}$ | $2.33\cdot10^8$ | $2.33\cdot10^8$ | $7.25\cdot10^6$ (majorant $7.07\cdot10^6$) | $3.94\cdot10^6$ | $2.78\cdot10^6$ |
| $10^{-6}$ | $2.33\cdot10^7$ | $2.33\cdot10^7$ | $7.25\cdot10^5$ | | $2.78\cdot10^5$ |
| $10^{-7}$ | | | | $3.93\cdot10^4$ | $2.78\cdot10^4$ |

The paper's values ($7.3\cdot10^6,7.3\cdot10^5,7.3\cdot10^4$; before $2.3\cdot10^8,2.3\cdot10^7,2.3\cdot10^6$) are reproduced with $N^v_k$ exact for $k\le7$ only, whereas
Remark `mg:rem:practical`(ii) now lists $N^v_8=2034$ among "the counts of (ii)".

**(D) threshold, $SU(3)$, $d=4$**: total loss $8(K+1)\theta+\log(\bar\Upsilon/\Upsilon)+\log(1+\delta')$ over $\beta^3/432$.
With the paper's counts: $0.125$, $0.282$, $0.501$, $0.783$, $1.127$ at $\beta=2,3,4,5,6\cdot10^{-8}$ (quadratic in $\beta$): crossing at $\beta\approx5.7\cdot10^{-8}$.
Before ($\theta_5$, general counts): $0.174$, $0.695$, $1.564$ at $1,2,3\cdot10^{-9}$: crossing at $\approx2.4\cdot10^{-9}$. With $N^v_8$: crossing $\approx7.7\cdot10^{-8}$; with
$N^v_8,N^v_9$: $0.769$, $0.973$, $1.201$ at $8,9,10\cdot10^{-8}$, crossing $\approx9.1\cdot10^{-8}$.

**(E) $SU(3)$, $d=3$** ($\theta_7$, $N^v_{k\le9}$, $Z_{n\le7,4}=0$): $\delta'/\beta^4=1.08\cdot10^4$, $3.23\cdot10^3$, $1.08\cdot10^3$ at $\beta=10^{-5}$, $3\cdot10^{-6}$, $10^{-6}$
(paper: $1.1\cdot10^4$, $3.2\cdot10^3$, $1.1\cdot10^3$); loss/$(\beta^3/432)=0.465$ at $10^{-6}$: crossing near $1.5\cdot10^{-6}$ (paper: resolved for $\beta\lesssim10^{-6}$).

**(G), (H), (I): Remark `mg:rem:practical`(iii) as revised at 22:34 UTC** ($\theta_7$, the counts of (ii), and $M_5$ of Lemma `mg:lem:jsums5` with $N_3,N_2,N'_2=876,320,312$
resp. $48,24,24$; `ref_numbers.py G`, `H`, `I`), $SU(3)$:

* $d=4$, $\beta=10^{-5}$: $\delta'/\beta^4=7.22\cdot10^4$ (paper $7.2\cdot10^4$); with $N^v_k$ exact only for $k\le7$ it would be $3.48\cdot10^5$ (the one-plaquette term with
  the crude $N^v_8$ then dominates).
* $d=4$, $\beta=1.4\cdot10^{-6}$: $\delta'/\beta^4=1.52\cdot10^3$ (paper $1.5\cdot10^3$) against $1/(432\beta)=1.65\cdot10^3$; total loss over $\beta^3/432$: $0.922$. Parts: majorant
  remainder $1.35\cdot10^3$, one-plaquette junctions $175$, the rest below $1$.
* $d=4$, $\beta=10^{-6}$: $\delta'/\beta^4=794$ (paper $7.9\cdot10^2$); loss ratio $0.343$. The ratio grows like $\beta^3$ (the remainder of $M_5$ is
  $O(\beta^6)$), so the crossing is at $\beta\approx1.44\cdot10^{-6}$, as claimed ($\beta\le1.4\cdot10^{-6}$). At $10^{-7}$: $\delta'/\beta^4=10.2$ (loss ratio
  $4.4\cdot10^{-4}$): $\delta'/\beta^4$ keeps decreasing with $\beta$, as the remark assumes.
* $d=3$, $\beta=1.7\cdot10^{-5}$: $\delta'/\beta^4=113$ (paper $1.1\cdot10^2$) against $136$; loss ratio $0.893$.
* With also $N^v_9=4344$ in $d=4$ (`I`): $\delta'/\beta^4=4.05\cdot10^4$ at $10^{-5}$ and $794$ at $1.4\cdot10^{-6}$ (loss ratio $0.48$; the one-plaquette term drops from
  $175$ to $1.2$), so the $d=4$ crossing would move to about $1.8\cdot10^{-6}$.
* $d=3$, $\beta=10^{-5}$ and $10^{-4}$: $\delta'/\beta^4=39.1$ and $3.91\cdot10^3$ (paper $39$, $3.9\cdot10^3$); loss ratio $0.173$ at $10^{-5}$, so the $d=3$ crossing is at
  $\beta\approx1.77\cdot10^{-5}$ (claimed: $\beta\le1.7\cdot10^{-5}$). In $d=3$ the whole of $\delta'$ is the remainder of $M_5$.
* Not evaluated: $d=4$ at $3\cdot10^{-6}$, $d=3$ at $3\cdot10^{-5}$ and $2\cdot10^{-5}$.

## Lemmas added to Section 16 during the review (`ref_three.py m case`)

Lemma `mg:lem:threelow` (junctions with exactly three front-pairs, one of them $p=(\partial q,\lambda)$), from the definitions (fronts as elements of a
deduplicated search, seeded at the vertices of $J$):

| case | $d=3$ | $d=4$ | paper |
|---|---|---|---|
| (i) $|P|=3$, fronts $4,4,4$ | $48$ (all caps) | $876$ ($864$ caps $=3\nu_2$, $12$ cube) | $N_3=48$, $876$ |
| (ii) $|P|=2$, fronts $4,4,6$ | $24$ | $320$ | $N_2=24$, $320$ |
| (iii) $|P|=2$, fronts $4,4,7$ | $24$ | $312$ | $N'_2=24$, $312$ |

Case (iii) was found here independently, as the class $(n,E)=(2,\frac72)$ missing from the first version of Lemma `mg:lem:jsums5`; the authors added it
($N'_2=K_S^2/2+K_SJ_v$) at the same time.

## Other checks

* `closed_sets.py` (authors) reproduced: $\{6{:}\,2\}$ ($\mathbb Z^3$), $\{6{:}\,4\}$ ($\mathbb Z^4$) up to 9 plaquettes.
* `junction_count.py` (authors) treats two front-pairs as never adjacent in $H_s$; Definition `mg:def:junction` makes them adjacent when their fronts share a vertex. For two four-edge
  front-pairs this never matters (the plaquettes covering the unshared edges at a common vertex share that vertex); `ref_junctions.py 3 6 nofrontadj` (that adjacency removed) gives the
  same $Z_{n,4}$ ($n\le6$, $d=4$).

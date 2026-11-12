# Independent checks of Section 13 of Part IV (the referee's code)

Written by the independent referee of Section 13 of `theory/mass-gap-strong-coupling.tex` ("Junctions of order $\beta^4$; the third-order term
with an effective constant"), 12 November 2026, and added with credit. A compass, not part of the proofs. The code works from the definitions
alone, in $\mathbb Z^n$ near the origin ($n=d-1$; a torus with all sides $\ge8$ looks like this locally):
* a front is any non-empty connected edge set in which every vertex has degree $\ge2$, found by a Redelmeier search (`lattice.py`), not just cycles;
* a junction (without support) must satisfy four conditions:
  * fronts of one layer are vertex-disjoint;
  * every front-pair touches a cluster;
  * the junction is connected in $H_s$;
  * $n(e)\ge2$ on $E(J)$.

| Script | What it checks | Output (in `outputs/`) |
|---|---|---|
| `check_i_iv.py 2 3 4` | Lemma 13.1(i) and (iv), by enumerating all multisets of plaquette boundaries that cover each of their links at least twice (a different search from ours). Controls: shifts and jumps $8$, $32$, $72$; $Z_{3,4}=0$. | `out_i_iv_n23.txt`, `out_i_iv_n4.txt`: (i) $0$, $12$, $24$ ($8(d-3)$ turns, $2(d-3)$ jumps, $2(d-3)$ U-turns); (iv) none |
| `check_ii_iii.py n K2 K3` | Lemma 13.1(ii) for second fronts with $\le$ K2 edges, and (iii) for further fronts with $\le$ K3 edges | `out_ii_iii_n2.txt` ($n=2$, $10$, $8$), `out_ii_iii_n3.txt` ($n=3$, $8$, $7$), `out_ii_n3_K9.txt` ($n=3$, $9$): exactly the cap and $K_S$, $K_S$, $4(P_v-1)-2K_S$ exits with $6$, $7$, $8$ edges; none with $5$, $9$, $10$; (iii) none with $|x|,|y|\le7$ ($d=4$), only $(4,4,8)$ for $d=3$ |
| `check_six.py n KMAX` | Lemma 13.3: the narrowings of six-edge fronts | `out_six_n2.txt`, `out_six_n3.txt`: $108$ dominos with exactly $2$ narrowings, $24$ skew hexagons with none ($n=3$) |
| `test_fronts.py`, `n7check.py` | the counts $N^v_4=P_v$, $N^v_6=12\binom{d-1}2+96\binom{d-1}3$, $N^v_7=\frac34P_vK_S$ | $N^v_6=12,132,456$, $N^v_7=12,108,360$ for $d=3,4,5$ |

Run from this directory (`lattice.py` is imported). Runtimes: seconds to a few minutes; `check_ii_iii.py 3 9 4` about 3 minutes.

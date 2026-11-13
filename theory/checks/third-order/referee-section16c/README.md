# Referee checks for Part IV, Section 16 (third pass): Lemma `mg:lem:corners`, Remarks `mg:rem:practical`(ii)–(iii), `mg:rem:nextorder`

Independent checks, made on 13 November, of these parts of `theory/mass-gap-strong-coupling.tex`:

* Remark `mg:rem:practical` (ii): the new counts.
* Remark `mg:rem:practical` (iii): the criterion with the whole deficit $D$.
* Lemma `mg:lem:corners` and the description of $S_6$ after it.
* Remark `mg:rem:nextorder` (i)–(iii). Parts (iv)–(v), the order after $S_6$ in $d=3$, were added during the review and were checked at the coordinator's request.

These checks are a compass, not part of the proofs. The enumerations were written from Definitions `mg:def:fronts`, `mg:def:junction` and Lemma `mg:lem:junction`. They use none of the
paper's lemmas and none of the authors' or earlier referees' code. The deficit is evaluated with the authors' evaluator (`../practical_beta.py`,
`../third_order_effective.py`). It is called through a wrapper that also reports which majorant is used and whether its series converge.

## Scripts

| Script | What it does | Result | Runtime |
|---|---|---|---|
| `fronts_redelmeier.c` | $N^v_k$ by Redelmeier enumeration of bond animals through the origin with a degree-deficiency prune (not a decision search); dumps the fronts for the junction search | $\mathbb Z^2$: $4,12,12,70,160,634,1728,6369$ ($k=4,6,\dots,12$); $\mathbb Z^3$: $12,132,108,2034,4344,43650$ ($k=4,6,\dots,10$) | 0.7 s, 5 s |
| `junction_search.py DUMP D n SIZES` | candidate junctions $(P,\Phi)\ni p=(\partial q,s-1)$ with $\lvert P\rvert=n$ and further fronts of the given sizes, satisfying (1)–(5) of the docstring (necessary conditions from the definitions). It adds one element per step: at the smallest link covered once, or by seeding at a vertex of a closed state. There are no exclusion sets, so memoising states is sound | see below | 0.1–36 s per class |
| `zexit_check.py D` | $Z_{2,k}$ by a second method: two plaquettes, exits as subsets of $\partial q\cup E(P)$ (from $n_J\ge2$), both layers tried, all conditions re-checked | $d=3$: $4,4,22,36,42,24,18$; $d=4$: $20,12,182,296,606,456,458$ ($k=6..12$); none with $k=5$; all exits in the other layer | 0.3 s |
| `closed_search.py n NMAX` | link-connected closed sets through a plaquette (one plaquette per step, no exclusions, memoised) | $\mathbb Z^3$: $\{6{:}2,10{:}10,11{:}11,12{:}24\}$ and none with 13; $\mathbb Z^4$: $\{6{:}4,10{:}60,11{:}66\}$; $\mathbb Z^5$ up to 9: $\{6{:}6\}$ | 2 s |
| `eval_thresholds.py MODE d beta..` | the authors' evaluator in the modes `s5` (default), `next`, `next7`, `next2`, plus `noM5` and `noM5sum` (without Lemma `mg:lem:jsums5`, with the maximum resp. the sum in $M_4$). Prints $D/\beta^4$, its parts, the majorant used and its $(\bar b,\bar a)$, and the ratios of the barred series to their radii of convergence | see below | 20–30 s per $\beta$ |

Junction counts (`outputs/junctions_*.txt`, `outputs/z_n4.txt`). Each pair of numbers is $d=3$, $d=4$, and each count agrees with the paper:

* Validation:
  * cap $1,1$;
  * one-plaquette exits with $6,7,8,9$ edges: $4,4,4,0$ and $12,12,20,0$;
  * $Z_{2,4}=8,32$;
  * $N_3=48,876$, $N_2=24,320$, $N'_2=24,312$;
  * $Z_{4,4}=0,12$ and $Z_{5,4}=0,26$;
  * $Z_{6,4}=302$ and $Z_{7,4}=288$ ($d=4$);
  * $Z_{n,4}=0$ for $3\le n\le10$ in $d=3$.
* $S_6\setminus S_5$:
  * $N_4=0,42$;
  * $N_{3,6}=24,680$;
  * $N_{3,7}=24,344$;
  * $N_{2,8}=48,992$;
  * $N_{2,444}=0,4$ (the cube with $p$ and its opposite face in layer $s-1$, a second pair of opposite faces in layer $s$, $P$ the third pair; $4(d-3)$);
  * $N_{2,9}=16,256$.
* $Z_{2,k}$, $k=6..9$: $4,4,22,36$ and $20,12,182,296$, by both methods.
* $S'_7$ ($d=3$, Remark (v)):
  * $0$ (five plaquettes, fronts $4,4,4$);
  * $276$ and $276$ (four plaquettes, fronts $4,4,6$ and $4,4,7$);
  * $252$, $0$, $424$ (three plaquettes);
  * $100,0,0$ and $212,0,0$ (two plaquettes).

Deficit, $SU(3)$ (`outputs/eval_*.txt`; $D/\beta^4$ against $\beta^3/432$ in units of $\beta^4$):

| mode | $d$ | values | paper |
|---|---|---|---|
| `s5` | 3 | $12.7$ ($10^{-5}$), $90.57<100.64$ ($2.3\cdot10^{-5}$), $101.4>96.5$ ($2.4\cdot10^{-5}$), $188$ ($3\cdot10^{-5}$), $11537$ ($10^{-4}$) | $13, 91, 1.9\cdot10^2, 1.2\cdot10^4$; threshold $2.3\cdot10^{-5}$ ✓ |
| `s5` | 4 | $128.9$ ($10^{-6}$), $871.4<890.3$ ($2.6\cdot10^{-6}$), $939.8>857.3$ ($2.7\cdot10^{-6}$), $1160$ ($3\cdot10^{-6}$), $12987$ ($10^{-5}$) | ✓ |
| `next` | 3 | $0.849$, $24.11<25.72$ ($9\cdot10^{-5}$), $28.7>24.4$ ($9.5\cdot10^{-5}$), $34.0$ ($10^{-4}$) | ✓ |
| `next` | 4 | $2.945$, $108.2$, $187.9<192.9$ ($1.2\cdot10^{-5}$), $239.6>178.1$ ($1.3\cdot10^{-5}$) | ✓ |
| `next7` | 3 | $73.1<79.8$ ($2.9\cdot10^{-5}$), $78.27<78.47$ ($2.95\cdot10^{-5}$), $83.7>77.2$ ($3\cdot10^{-5}$) | paper: threshold $2.8\cdot10^{-5}$ (conservative; it is $\approx2.95\cdot10^{-5}$) |
| `next7` | 4 | $206<231$ ($10^{-5}$), $288>210$ ($1.1\cdot10^{-5}$) | ✓ ($10^{-5}$) |
| `next2` | 3 | $5.21$ ($1.6\cdot10^{-4}$), $8.29<12.86$ ($1.8\cdot10^{-4}$), $12.58>11.57$ ($2\cdot10^{-4}$) | ✓ |
| `noM5sum` / `noM5` | 4 | $23760>23148$ at $10^{-7}$ (crossing $\approx0.99\cdot10^{-7}$) / crossing $\approx1.4\cdot10^{-7}$ | "about $10^{-7}$" ✓ |
| `noM5sum` / `noM5` | 3 | $1081<1543$ at $1.5\cdot10^{-6}$ (crossing $\approx1.8\cdot10^{-6}$) / $\approx2.2\cdot10^{-6}$ | paper: "about $1.5\cdot10^{-6}$" (low) |

The θ-loss is $1.73\beta^4$ at $9\cdot10^{-5}$ ($d=3$) and $2.6\cdot10^{-3}\beta^4$ at $1.2\cdot10^{-5}$ ($d=4$), as Remark `mg:rem:nextorder`(ii) states. At every selected point the
majorant is the corner (or $M_5$) majorant, at $(\bar b,\bar a)=(2^{-7},2^{-12})$ ($d=3$) or $(2^{-8},2^{-13})$ ($d=4$), resp. $(2^{-9},2^{-14})$ for `s5` in $d=4$. There the
barred series are well inside their disks of convergence (ratios $0.07$–$0.16$ for $\bar Y$ and $0.21$–$0.31$ for $\bar A$), so the truncations at `NMAX`, `KMAX` are harmless.

## Commands

```
gcc -O2 -o $S/fronts_redelmeier fronts_redelmeier.c          # S: a scratch directory
$S/fronts_redelmeier 2 12; $S/fronts_redelmeier 3 10
$S/fronts_redelmeier 2 11 $S/fronts_D2_k11.txt 4; $S/fronts_redelmeier 3 9 $S/fronts_D3_k9.txt 4
python3 junction_search.py $S/fronts_D2_k11.txt 2 4 4,4      # etc.; one call per n, several size classes allowed
python3 junction_search.py $S/fronts_D3_k9.txt 3 2 4,7 4,8 4,4,4 4,9
python3 zexit_check.py 2; python3 zexit_check.py 3
python3 closed_search.py 3 13; python3 closed_search.py 4 11; python3 closed_search.py 5 9
python3 eval_thresholds.py s5 3 2.3e-5 2.4e-5                 # and the other modes and beta of the table
```

Total CPU time of all runs: about 20 minutes.

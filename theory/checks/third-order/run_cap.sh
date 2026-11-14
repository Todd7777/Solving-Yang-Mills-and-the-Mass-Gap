#!/bin/bash
# Part IV, Section 16 (Theorem mg:thm:thirdcap, Corollaries mg:cor:sectorscap, mg:cor:classic): every evaluation of the computer-assisted proofs, with the exact commands, written
# to orders/ with the command line as the first line of each output. Run from theory/checks/third-order:  bash run_cap.sh [JOBS] [theorem|corollary|classic]
# (default all three; see the README for the run times; the d=3 runs of the corollaries sum the fronts to 100 edges and take longer than those of the theorem)
# beta_* = beta_0(d)(1 + 1e-30) >= beta_0(d), as decimal strings (the ninth referee's run_reproductions.sh; checked >= beta_0(d) in Arb by referee-section16i/arb_eval.py).
set -e
J=${1:-2}
B3=0.00040262143352281974401285396895067832618294529962141   # beta_0(3)(1+1e-30)
B4=0.00015376644804254827379172376907999414111262406766301   # beta_0(4)(1+1e-30)
CL="--closed 10:60,11:66,12:320,13:0,14:1316,15:2520,16:12552 --closed-from 17"
PART=${2:-all}
run() { out=$1; shift; { echo "# python3 $*"; python3 "$@"; } > "orders/$out" 2>&1; echo "orders/$out: $(tail -1 orders/$out | cut -c1-160)"; }
# Theorem mg:thm:thirdcap: the lower bound (scaling from beta_*), the error budget, the upper bound, the group bounds
if [ "$PART" = all ] || [ "$PART" = theorem ]; then
run cap_su3_d3_scaling.out scaling_check.py --bounds --certified $B3 $J
run cap_su5_d3_scaling.out scaling_check.py --bounds --certified --group 5 $B3 $J
run cap_su3_d4_scaling.out scaling_check.py --bounds --certified --d4 orders/d4_s10x_92.json $CL $B4 $J
run cap_su5_d4_scaling.out scaling_check.py --bounds --certified --group 5 --d4 orders/d4_s10x_92.json $CL $B4 $J
run error_budget.out error_budget.py
run cap_upper.out certify_upper.py
run rigorous_bounds.out rigorous_bounds.py
fi
# Corollary mg:cor:sectorscap: d=3 with the empty kinds of Lemma mg:lem:planar(ii) (orders/d3_s9z.json, from orders/assemble_d3_s9z.py) and the fronts of Z^2 beyond
# the exact counts by Lemma mg:lem:frontsdecision (--nv-decision-d3, the eleventh referee's Issue 1; then the optimal abar is small and the fronts are summed to 100 edges),
# SU(4) in d=4, the error budget of the new cases, the sector upper bounds, and the series and splitting (rigorous_bounds.py, above)
if [ "$PART" = all ] || [ "$PART" = corollary ]; then
run capz_su3_d3_scaling.out scaling_check.py --bounds --certified --nv-decision-d3 --kmax 100 --d3 orders/d3_s9z.json $B3 $J
run capz_su4_d3_scaling.out scaling_check.py --bounds --certified --group 4 --nv-decision-d3 --kmax 100 --d3 orders/d3_s9z.json $B3 $J
run capz_su5_d3_scaling.out scaling_check.py --bounds --certified --group 5 --nv-decision-d3 --kmax 100 --d3 orders/d3_s9z.json $B3 $J
run cap_su4_d4_scaling.out scaling_check.py --bounds --certified --group 4 --d4 orders/d4_s10x_92.json $CL $B4 $J
run error_budget_sectors.out error_budget.py --groups 3,4,5 --dims 4
run error_budget_sectors_d3.out error_budget.py --d3 orders/d3_s9z.json --groups 3,4,5 --dims 3 --nv-decision-d3 --kmax 100
run sector_upper.out sector_upper.py
fi
# Corollary mg:cor:classic: SU(2) (--group 2, iota = 1) and U(1) (--group 1) in bound mode with u_2 <= e^beta - 1 - beta. For U(1) in d=3 the constants F, F_6
# are raised by 1e-7 (needed with the binomial front bound, where delta'_* ~ 5e-9; with Lemma mg:lem:frontsdecision delta'_* ~ 5e-12 and either would do; the
# raise keeps the stronger form of Lemma mg:lem:scaling with room); the fronts are summed to 100 edges (the optimal abar is small,
# 2^-8 in d=4, and the barred front sums converge slowly: at 60 edges their tails would be 1e-13 of the sums, above the room of the six inequalities); in d=3 the
# fronts of Z^2 beyond the exact counts by Lemma mg:lem:frontsdecision, as for Corollary mg:cor:sectorscap
if [ "$PART" = all ] || [ "$PART" = classic ]; then
run capz_su2_d3_scaling.out scaling_check.py --bounds --certified --group 2 --nv-decision-d3 --kmax 100 --d3 orders/d3_s9z.json $B3 $J
run cap_su2_d4_scaling.out scaling_check.py --bounds --certified --group 2 --kmax 100 --d4 orders/d4_s10x_92.json $CL $B4 $J
run capz_u1_d3_scaling.out scaling_check.py --bounds --certified --group 1 --f-raise 1e-7 --nv-decision-d3 --kmax 100 --d3 orders/d3_s9z.json $B3 $J
run cap_u1_d4_scaling.out scaling_check.py --bounds --certified --group 1 --kmax 100 --d4 orders/d4_s10x_92.json $CL $B4 $J
run error_budget_su2_d3.out error_budget.py --d3 orders/d3_s9z.json --groups 2 --dims 3 --nv-decision-d3 --kmax 100
run error_budget_su2_d4.out error_budget.py --groups 2 --dims 4 --kmax 100
run error_budget_u1_d3.out error_budget.py --d3 orders/d3_s9z.json --groups 1 --dims 3 --f-raise 1e-7 --nv-decision-d3 --kmax 100
run error_budget_u1_d4.out error_budget.py --groups 1 --dims 4 --kmax 100
run sector_upper.out sector_upper.py
run rigorous_bounds.out rigorous_bounds.py
fi

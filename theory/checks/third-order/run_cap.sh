#!/bin/bash
# Part IV, Section 16 (Theorem mg:thm:thirdcap): every evaluation of the computer-assisted proof, with the exact commands, written
# to orders/ with the command line as the first line of each output. Run from theory/checks/third-order:  bash run_cap.sh [JOBS]   (about 15 minutes with 2 jobs)
# beta_* = beta_0(d)(1 + 1e-30) >= beta_0(d), as decimal strings (the ninth referee's run_reproductions.sh; checked >= beta_0(d) in Arb by referee-section16i/arb_eval.py).
set -e
J=${1:-2}
B3=0.00040262143352281974401285396895067832618294529962141   # beta_0(3)(1+1e-30)
B4=0.00015376644804254827379172376907999414111262406766301   # beta_0(4)(1+1e-30)
CL="--closed 10:60,11:66,12:320,13:0,14:1316,15:2520,16:12552 --closed-from 17"
run() { out=$1; shift; { echo "# python3 $*"; python3 "$@"; } > "orders/$out" 2>&1; echo "orders/$out: $(tail -1 orders/$out | cut -c1-160)"; }
# Theorem mg:thm:thirdcap: the lower bound (scaling from beta_*), the error budget, the upper bound, the group bounds
run cap_su3_d3_scaling.out scaling_check.py --bounds --certified $B3 $J
run cap_su5_d3_scaling.out scaling_check.py --bounds --certified --group 5 $B3 $J
run cap_su3_d4_scaling.out scaling_check.py --bounds --certified --d4 orders/d4_s10x_92.json $CL $B4 $J
run cap_su5_d4_scaling.out scaling_check.py --bounds --certified --group 5 --d4 orders/d4_s10x_92.json $CL $B4 $J
run error_budget.out error_budget.py
run cap_upper.out certify_upper.py
run rigorous_bounds.out rigorous_bounds.py

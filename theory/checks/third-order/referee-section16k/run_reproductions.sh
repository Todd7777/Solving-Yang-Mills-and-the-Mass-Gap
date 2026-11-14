#!/bin/bash
# Referee (eleventh pass on Section 16, Corollary mg:cor:classic): rerun the "classic" part of run_cap.sh, with every output written to
# referee-section16k/outputs/ instead of orders/ (so that nothing outside this folder is changed), and without writing bytecode caches.
# Run from theory/checks/third-order:  bash referee-section16k/run_reproductions.sh [JOBS]
# The four scaling checks run two at a time (each with JOBS workers), then the error budgets in parallel.
export PYTHONDONTWRITEBYTECODE=1
J=${1:-2}
O=referee-section16k/outputs
B3=0.00040262143352281974401285396895067832618294529962141   # beta_0(3)(1+1e-30)
B4=0.00015376644804254827379172376907999414111262406766301   # beta_0(4)(1+1e-30)
CL="--closed 10:60,11:66,12:320,13:0,14:1316,15:2520,16:12552 --closed-from 17"
run() { out=$1; shift; { echo "# python3 $*"; ( time python3 "$@" ) 2>&1; } > "$O/$out"; echo "$O/$out: $(grep -v '^real\|^user\|^sys\|^$' $O/$out | tail -1 | cut -c1-160)"; }
run repro_capz_su2_d3_scaling.out scaling_check.py --bounds --certified --group 2 --d3 orders/d3_s9z.json $B3 $J &
run repro_capz_u1_d3_scaling.out scaling_check.py --bounds --certified --group 1 --f-raise 1e-7 --d3 orders/d3_s9z.json $B3 $J &
wait
run repro_cap_su2_d4_scaling.out scaling_check.py --bounds --certified --group 2 --kmax 100 --d4 orders/d4_s10x_92.json $CL $B4 $J &
run repro_cap_u1_d4_scaling.out scaling_check.py --bounds --certified --group 1 --kmax 100 --d4 orders/d4_s10x_92.json $CL $B4 $J &
wait
run repro_error_budget_su2_d3.out error_budget.py --d3 orders/d3_s9z.json --groups 2 --dims 3 &
run repro_error_budget_su2_d4.out error_budget.py --groups 2 --dims 4 --kmax 100 &
run repro_error_budget_u1_d3.out error_budget.py --d3 orders/d3_s9z.json --groups 1 --dims 3 --f-raise 1e-7 &
run repro_error_budget_u1_d4.out error_budget.py --groups 1 --dims 4 --kmax 100 &
wait
run repro_sector_upper.out sector_upper.py
run repro_rigorous_bounds.out rigorous_bounds.py
# addendum (the coordinator's follow-up): the d = 3 runs with --nv-decision-d3 (Lemma mg:lem:frontsdecision for Z^2) for SU(3), SU(4), SU(5), and their error budget
run repro_nvdec_capz_su3_d3_scaling.out scaling_check.py --bounds --certified --nv-decision-d3 --kmax 100 --d3 orders/d3_s9z.json $B3 $J &
run repro_nvdec_capz_su4_d3_scaling.out scaling_check.py --bounds --certified --group 4 --nv-decision-d3 --kmax 100 --d3 orders/d3_s9z.json $B3 $J &
wait
run repro_nvdec_capz_su5_d3_scaling.out scaling_check.py --bounds --certified --group 5 --nv-decision-d3 --kmax 100 --d3 orders/d3_s9z.json $B3 $J &
run repro_nvdec_error_budget_sectors_d3.out error_budget.py --d3 orders/d3_s9z.json --groups 3,4,5 --dims 3 --nv-decision-d3 --kmax 100 &
wait

#!/bin/bash
# Referee (section 16, ninth pass): rerun the four bound-mode scaling checks, the error budget and certify_upper (Theorem mg:thm:thirdcap).
# Run from theory/checks/third-order:  bash referee-section16i/run_reproductions.sh [JOBS]
J=${1:-2}
O=referee-section16i/outputs
B3=0.00040262143352281974401285396895067832618294529962141   # beta_0(3)(1+1e-30)
B4=0.00015376644804254827379172376907999414111262406766301   # beta_0(4)(1+1e-30)
CL="--closed 10:60,11:66,12:320,13:0,14:1316,15:2520,16:12552 --closed-from 17"
( time python3 scaling_check.py --bounds --certified $B3 $J ) > $O/repro_cap_su3_d3_scaling.out 2>&1
( time python3 scaling_check.py --bounds --certified --group 5 $B3 $J ) > $O/repro_cap_su5_d3_scaling.out 2>&1
( time python3 scaling_check.py --bounds --certified --d4 orders/d4_s10x_92.json $CL $B4 $J ) > $O/repro_cap_su3_d4_scaling.out 2>&1
( time python3 scaling_check.py --bounds --certified --group 5 --d4 orders/d4_s10x_92.json $CL $B4 $J ) > $O/repro_cap_su5_d4_scaling.out 2>&1
( time python3 certify_upper.py ) > $O/repro_cap_upper.out 2>&1
( time python3 error_budget.py ) > $O/repro_error_budget.out 2>&1

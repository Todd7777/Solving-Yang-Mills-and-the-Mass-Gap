#!/bin/bash
# Referee (tenth pass on Section 16, Corollary mg:cor:sectorscap): rerun the corollary part of run_cap.sh, with every output written to
# referee-section16j/outputs/ instead of orders/ (so that nothing outside this folder is changed), and without writing bytecode caches.
# Run from theory/checks/third-order:  bash referee-section16j/run_reproductions.sh [JOBS]
export PYTHONDONTWRITEBYTECODE=1
J=${1:-2}
O=referee-section16j/outputs
B3=0.00040262143352281974401285396895067832618294529962141   # beta_0(3)(1+1e-30)
B4=0.00015376644804254827379172376907999414111262406766301   # beta_0(4)(1+1e-30)
CL="--closed 10:60,11:66,12:320,13:0,14:1316,15:2520,16:12552 --closed-from 17"
run() { out=$1; shift; { echo "# python3 $*"; ( time python3 "$@" ) 2>&1; } > "$O/$out"; echo "$O/$out: $(grep -v '^real\|^user\|^sys\|^$' $O/$out | tail -1 | cut -c1-160)"; }
run repro_capz_su3_d3_scaling.out scaling_check.py --bounds --certified --d3 orders/d3_s9z.json $B3 $J
run repro_capz_su4_d3_scaling.out scaling_check.py --bounds --certified --group 4 --d3 orders/d3_s9z.json $B3 $J
run repro_capz_su5_d3_scaling.out scaling_check.py --bounds --certified --group 5 --d3 orders/d3_s9z.json $B3 $J
run repro_cap_su4_d4_scaling.out scaling_check.py --bounds --certified --group 4 --d4 orders/d4_s10x_92.json $CL $B4 $J
run repro_error_budget_sectors.out error_budget.py --d3 orders/d3_s9z.json --groups 3,4,5
run repro_sector_upper.out sector_upper.py
run repro_rigorous_bounds.out rigorous_bounds.py

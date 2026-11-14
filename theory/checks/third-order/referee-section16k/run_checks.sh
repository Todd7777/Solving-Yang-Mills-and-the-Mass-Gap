#!/bin/bash
# Referee (eleventh pass on Section 16, Corollary mg:cor:classic): the referee's own checks. Run from theory/checks/third-order/referee-section16k:
#   bash run_checks.sh          (about 1 minute on 4 cores; run_reproductions.sh, from theory/checks/third-order, reruns the authors' commands: about 3 minutes)
# Every output goes to outputs/ with the command line as its first line. No file outside this folder is written (no bytecode caches either).
export PYTHONDONTWRITEBYTECODE=1
set -e
run() { out=$1; shift; { echo "# python3 $*"; ( time python3 "$@" ) 2>&1; } > "outputs/$out"; echo "outputs/$out: $(grep -v '^real\|^user\|^sys\|^$' outputs/$out | tail -1 | cut -c1-150)"; }
# the lower bound in ball arithmetic, the corollary's four cases (GROUP d c_-), and the d = 3 suggestion (Lemma mg:lem:frontsdecision for Z^2)
run arb_eval_classic_su2_d3.out arb_eval_classic.py 2 3 1968 &
run arb_eval_classic_su2_d4.out arb_eval_classic.py 2 4 570 &
run arb_eval_classic_u1_d3.out arb_eval_classic.py 1 3 177200 &
run arb_eval_classic_u1_d4.out arb_eval_classic.py 1 4 8043 &
wait
run arb_eval_decision_d3_g2.out arb_eval_classic.py 2 3 12.7 --decision-d3 &
run arb_eval_decision_d3_g1.out arb_eval_classic.py 1 3 181 --decision-d3 &
wait
# addendum (the coordinator's follow-up): SU(3), SU(4), SU(5) in d = 3 with Lemma mg:lem:frontsdecision for Z^2 (Corollary mg:cor:sectorscap), and, as a check of the
# extension to N >= 3, the same groups without the option against the tenth referee's certified values (0.5993437096, 0.1332401769, 0.07288809269)
run arb_eval_decision_d3_g3.out arb_eval_classic.py 3 3 0.185 --decision-d3 &
run arb_eval_decision_d3_g4.out arb_eval_classic.py 4 3 0.094 --decision-d3 &
run arb_eval_decision_d3_g5.out arb_eval_classic.py 5 3 0.067 --decision-d3 &
wait
run arb_eval_validate_su3_d3_s9z.out arb_eval_classic.py 3 3 0.600 &
run arb_eval_validate_su4_d3_s9z.out arb_eval_classic.py 4 3 0.134 &
run arb_eval_validate_su5_d3_s9z.out arb_eval_classic.py 5 3 0.073 &
wait
# the upper bounds in ball arithmetic; the series, the Turan determinant and the splitting; the inputs of bound mode with the exact group data
run upper_arb_classic.out upper_arb_classic.py
run series_check_classic.out series_check_classic.py
run applicability_check.out applicability_check.py
# exact group data against bound mode (Numerical); the d = 3 probe and scaling with the decision bound (suggestion)
for c in "2 3" "2 4" "1 3" "1 4"; do set -- $c; run exact_vs_bound_g$1_d$2.out exact_vs_bound.py $c & done; wait
for c in "2 on" "1 on" "3 on" "3 off"; do set -- $c; run nv_probe_g$1_$2.out nv_decision_d3_probe.py $c & done; wait
run scaling_decision_d3_g2.out scaling_decision_d3.py 2 2 &
run scaling_decision_d3_g1.out scaling_decision_d3.py 1 2 &
wait
# the two technical points of the proof, tested by the authors' error budget with the defaults (both must FAIL): kmax 60 in d = 4, the F raise 1e-9 in d = 3
cd ..
{ echo "# python3 error_budget.py --groups 2 --dims 4 --kmax 60   (referee: d=4 at the default truncation)"; ( time python3 error_budget.py --groups 2 --dims 4 --kmax 60 ) 2>&1; } > referee-section16k/outputs/error_budget_kmax60_g2.out &
{ echo "# python3 error_budget.py --groups 1 --dims 4 --kmax 60   (referee: d=4 at the default truncation)"; ( time python3 error_budget.py --groups 1 --dims 4 --kmax 60 ) 2>&1; } > referee-section16k/outputs/error_budget_kmax60_g1.out &
{ echo "# python3 error_budget.py --d3 orders/d3_s9z.json --groups 1 --dims 3   (referee: U(1), d=3 with the default F raise 1e-9)"; ( time python3 error_budget.py --d3 orders/d3_s9z.json --groups 1 --dims 3 ) 2>&1; } > referee-section16k/outputs/error_budget_u1_d3_default_raise.out &
wait
cd referee-section16k
grep -h 'smallest room' outputs/error_budget_kmax60_g*.out outputs/error_budget_u1_d3_default_raise.out
# the quoted numbers against the outputs above
run numbers_check_classic.out numbers_check_classic.py

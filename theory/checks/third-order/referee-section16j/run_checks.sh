#!/bin/bash
# Referee (tenth pass on Section 16, Corollary mg:cor:sectorscap): the referee's own checks. Run from theory/checks/third-order/referee-section16j:
#   bash run_checks.sh          (about 2 minutes on 4 cores; run_reproductions.sh, from theory/checks/third-order, reruns the authors' commands: about 5 minutes)
# Every output goes to outputs/ with the command line as its first line. No file outside this folder is written (no bytecode caches either).
export PYTHONDONTWRITEBYTECODE=1
set -e
run() { out=$1; shift; { echo "# python3 $*"; ( time python3 "$@" ) 2>&1; } > "outputs/$out"; echo "outputs/$out: $(grep -v '^real\|^user\|^sys\|^$' outputs/$out | tail -1 | cut -c1-150)"; }
run corners_check.out corners_check.py
run series_check.out series_check.py
run split_check_arb.out split_check_arb.py
run weyl_check.out weyl_check.py
run upper_arb_sectors.out upper_arb_sectors.py
run sector_upper_compare.out sector_upper_compare.py
# the lower bound in ball arithmetic: the corollary's four new cases, and (for the comparison of the majorant, and the stronger form of Lemma mg:lem:scaling)
# the four cases of Theorem mg:thm:thirdcap
for c in "3 3 s9z 0.600" "4 3 s9z 0.134" "5 3 s9z 0.073" "4 4 x92 1.083" "3 3 s9 2.558" "5 3 s9 0.346" "3 4 x92 3.937" "5 4 x92 0.520"; do
  set -- $c; run arb_eval_cor_su$1_d$2_$3.out arb_eval_cor.py $c &
done; wait
run numbers_check.out numbers_check.py

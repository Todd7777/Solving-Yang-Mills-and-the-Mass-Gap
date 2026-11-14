#!/bin/bash
# Referee (Section 16h): rerun the authors' evaluators, two lanes (at most 2 cores), logs in $SCR/logs.
SCR=/tmp/claude-0/-home-user/e2720f80-9d25-5995-bdae-a9bca5a2da7f/scratchpad/referee16h
TO=/home/user/Solving-Yang-Mills-and-the-Mass-Gap/theory/checks/third-order
export PYTHONPYCACHEPREFIX=$SCR/pycache
mkdir -p $SCR/logs
B3=0.00040262143352281974401
B4=0.00015376644804254827379
C3="--closed 10:10,11:11,12:24,13:0,14:70 --closed-from 15"
C4="--closed 10:60,11:66,12:320,13:0,14:1316,15:2520,16:12552 --closed-from 17"
run() { name=$1; shift; echo "# $*" > $SCR/logs/$name.out; ( time nice -n 5 python3 "$@" ) >> $SCR/logs/$name.out 2>&1; echo "done $name" >> $SCR/logs/lanes.log; }
laneA() {
  cd $SCR
  run su5_d3_beta0     $TO/practical_beta.py --group 5 --terms $TO/orders/d3_s9.json $C3 3 $B3
  run su5_d4_92_beta0  $TO/practical_beta.py --group 5 --terms $TO/orders/d4_s10x_92.json $C4 4 $B4
  run su5_d4_s10x_beta0 $TO/practical_beta.py --group 5 --terms $TO/orders/d4_s10x.json $C4 4 $B4
  run su5_monotone_d3  $TO/monotone_check.py --group 5
  run su5_monotone_d4  $TO/monotone_check.py --group 5 --d4
  run su5_scaling_d3   $TO/scaling_check.py --group 5 $B3 1
  run su5_scaling_d4   $TO/scaling_check.py --group 5 --d4 $TO/orders/d4_s10x_92.json $C4 $B4 1
  run su5_upper_d3     $TO/upper_exact.py --group 5 3
  run sun_group        $TO/sun_group.py
  echo "laneA finished" >> $SCR/logs/lanes.log
}
laneB() {
  cd $SCR
  run su3_d4_s10x_beta0 $TO/practical_beta.py --terms $TO/orders/d4_s10x.json $C4 4 $B4
  run su3_d3_s9_beta0   $TO/practical_beta.py --terms $TO/orders/d3_s9.json $C3 3 $B3
  run su3_d4_92_beta0   $TO/practical_beta.py --terms $TO/orders/d4_s10x_92.json $C4 4 $B4
  run su3_monotone_d3   $TO/monotone_check.py
  run su3_monotone_d4   $TO/monotone_check.py --d4
  run su3_scaling_d4_92 $TO/scaling_check.py --d4 $TO/orders/d4_s10x_92.json $C4 $B4 1
  run su3_upper_d3      $TO/upper_exact.py 3
  run su3_upper_d4      $TO/upper_exact.py 4
  run su3_scaling_d3    $TO/scaling_check.py $B3 1
  echo "laneB finished" >> $SCR/logs/lanes.log
}
laneA &
laneB &
wait
echo ALL_LANES_DONE >> $SCR/logs/lanes.log

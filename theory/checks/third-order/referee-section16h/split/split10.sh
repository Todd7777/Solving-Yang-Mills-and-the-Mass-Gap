#!/bin/bash
cd /tmp/claude-0/-home-user/e2720f80-9d25-5995-bdae-a9bca5a2da7f/scratchpad/referee16h/split
for K in A Q; do
  s=$(date +%s); ./admissible_split 4 $K 10 2 > un_${K}_10.txt; e=$(date +%s); echo "unsplit $K 10: $((e-s)) s" >> times.txt
  if [ $K = A ]; then J=4; else J=16; fi
  s=$(date +%s); for k in $(seq 0 $((J-1))); do ./admissible_split 4 $K 10 2 $k $J > pc_${K}_10_${k}_${J}.txt; done; e=$(date +%s)
  echo "split $K 10 J=$J done: $((e-s)) s" >> times.txt
done
echo SPLIT10_DONE >> times.txt

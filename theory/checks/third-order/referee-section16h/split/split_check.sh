#!/bin/bash
# Referee: the pieces of admissible_split.c add up to the unsplit counts (Z^4, S_A and S_Q, nmax 8 and 9, several J).
cd /tmp/claude-0/-home-user/e2720f80-9d25-5995-bdae-a9bca5a2da7f/scratchpad/referee16h/split
for NM in 8 9; do
  for K in A Q; do
    s=$(date +%s); ./admissible_split 4 $K $NM 2 > un_${K}_${NM}.txt; e=$(date +%s); echo "unsplit $K $NM: $((e-s)) s" >> times.txt
    if [ $K = A ]; then JS="2 3 4"; else JS="2 3 16"; fi
    for J in $JS; do
      for k in $(seq 0 $((J-1))); do ./admissible_split 4 $K $NM 2 $k $J > pc_${K}_${NM}_${k}_${J}.txt; done
      echo "split $K $NM J=$J done" >> times.txt
    done
  done
done
echo SPLIT_DONE >> times.txt

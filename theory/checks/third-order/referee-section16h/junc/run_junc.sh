#!/bin/bash
# Referee: reproduce the class 4,4,4 count to nine plaquettes in Z^3 (Remark mg:rem:su5(iii)) with the two existing programs.
J=/tmp/claude-0/-home-user/e2720f80-9d25-5995-bdae-a9bca5a2da7f/scratchpad/referee16h/junc
TO=/home/user/Solving-Yang-Mills-and-the-Mass-Gap/theory/checks/third-order
cd $J
gcc -O2 -o junctions $TO/junctions.c 2> cc1.log
gcc -O2 -o fronts_enum $TO/referee-section16f/fronts_enum.c 2> cc2.log
gcc -O2 -o junc_partition_big $TO/referee-section16g/junc_partition_big.c 2> cc3.log
s=$(date +%s); ./junctions 3 9 4,4 > junctions_3_9_44.out 2>&1; e=$(date +%s); echo "junctions.c: $((e-s)) s" >> times.txt
./fronts_enum 3 4 $J/fronts_D3_k4.txt > fronts_enum.out 2>&1
s=$(date +%s); ./junc_partition_big $J/fronts_D3_k4.txt 3 9 4,4 > partition_3_9_44.out 2>&1; e=$(date +%s); echo "partition: $((e-s)) s" >> times.txt
echo JUNC_DONE >> times.txt

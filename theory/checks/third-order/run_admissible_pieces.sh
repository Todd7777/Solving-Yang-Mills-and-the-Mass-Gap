#!/bin/bash
# Remark mg:rem:su5(iv): the admissible sets of 11 plaquettes in Z^4, for (S_A,{p0}) in 4 pieces and (S_Q,Q) in 16, 4 at a time, resumably
# (each piece writes OUT/KIND_K_J.txt; finished pieces are skipped).  Usage: ./run_admissible_pieces.sh OUTDIR     (then: python3 sum_admissible_pieces.py OUTDIR)
OUT=${1:?usage: $0 OUTDIR}; mkdir -p "$OUT"; HERE=$(cd "$(dirname "$0")" && pwd)
[ -x "$OUT/admissible_split" ] || cc -O2 -o "$OUT/admissible_split" "$HERE/admissible_split.c"
run_one() { kind=$1; k=$2; J=$3; f=$OUT/${kind}_${k}_${J}.txt
  if grep -q "admissible sets by size" "$f" 2>/dev/null; then return; fi
  ( time "$OUT/admissible_split" 4 $kind 11 2 $k $J ) > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"; }
export -f run_one; export OUT
{ for k in 0 1 2 3; do echo "A $k 4"; done; for k in $(seq 0 15); do echo "Q $k 16"; done; } | xargs -P 4 -L 1 bash -c 'run_one $0 $1 $2'
echo ALL_PIECES_DONE

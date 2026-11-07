#!/bin/bash
# Slurm array job for the glueball / Wilson-loop production runs (one Markov chain per task).
#
#   cd $REPO/lattice/slurm
#   N=$(grep -vc '^#' params_su2.txt)
#   sbatch --array=1-$N glueball_array.sbatch            # uses params_su2.txt
#   PARAMS=params_su3.txt sbatch --array=1-6 glueball_array.sbatch
#
# Every run checkpoints (configuration + RNG state) every 50 measurements and resumes
# bit-for-bit when restarted, so a preemptible/requeue partition is safe: a requeued task
# simply continues.  Check partition names, time limits and the Python module against the
# current FASRC documentation for your account before submitting.
#SBATCH -J ym-glueball
#SBATCH -p serial_requeue
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 2-00:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH -o logs/%x_%A_%a.out
set -euo pipefail

SLURM_DIR=${SLURM_SUBMIT_DIR:-$(pwd)}
REPO=${REPO:-$(cd "$SLURM_DIR/../.." && pwd)}
PARAMS=${PARAMS:-params_su2.txt}
OUTDIR=${OUTDIR:-$REPO/results/fasrc}
mkdir -p "$OUTDIR" "$SLURM_DIR/logs"

line=$(grep -v '^#' "$SLURM_DIR/$PARAMS" | grep -v '^[[:space:]]*$' | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "$line" ]; then echo "no parameter line ${SLURM_ARRAY_TASK_ID} in $PARAMS"; exit 1; fi
read -r GROUP BETA L0 L NMEAS NOR LEVELS SEED MAXLOOP WEVERY <<< "$line"

# Python: a virtual environment with the package installed (see slurm/README.md).
if [ -n "${VENV:-}" ]; then source "$VENV/bin/activate"; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

OUT="$OUTDIR/gb_${GROUP}_b${BETA}_${L}x${L0}_s${SEED}.npz"
echo "task ${SLURM_ARRAY_TASK_ID}: $line -> $OUT ($(date))"
cd "$REPO/lattice"
python scripts/glueball.py --group "$GROUP" --beta "$BETA" --shape "$L0" "$L" "$L" "$L" \
    --nmeas "$NMEAS" --nor "$NOR" --levels ${LEVELS//,/ } --seed "$SEED" \
    --maxloop "$MAXLOOP" --wloop-every "$WEVERY" --out "$OUT"
echo "done $(date)"

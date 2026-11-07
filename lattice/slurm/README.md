# Running the production Monte Carlo on a Slurm cluster (e.g. FASRC Cannon)

The cloud container used for development has 4 cores, which is enough for validation but not for
the volumes needed by milestone T1.3 (the 0++ mass gap versus asymptotic scaling). These scripts
run the same code on a cluster, one Markov chain per array task.

## One-time setup (on a login node)

```bash
git clone https://github.com/Todd7777/Solving-Yang-Mills-and-the-Mass-Gap.git
cd Solving-Yang-Mills-and-the-Mass-Gap
git checkout claude/blissful-wozniak-ytnciq
module load python            # or the specific python/3.x module your cluster provides
python -m venv ~/venvs/ymlattice
source ~/venvs/ymlattice/bin/activate
pip install -e "lattice[test]"
cd lattice && pytest -q && cd ..   # about 2-3 minutes; all tests must pass before production
```

## Submit

```bash
cd lattice/slurm
VENV=~/venvs/ymlattice sbatch --export=ALL --array=1-12 glueball_array.sbatch
VENV=~/venvs/ymlattice PARAMS=params_su3.txt sbatch --export=ALL --array=1-6 glueball_array.sbatch
```

* Each line of `params_su2.txt` / `params_su3.txt` is one task (group, beta, L0, L, measurements,
  overrelaxation sweeps, smearing levels, seed, Wilson-loop size, Wilson-loop frequency).
* Tasks checkpoint every 50 measurements and resume exactly when requeued, so the preemptible
  `serial_requeue` partition is safe. Use `-p shared` (or your lab partition) to avoid preemption.
* Expected cost (pure NumPy, one core): roughly 1 s per measurement on 8^3 x 16, 5 s on 12^3 x 24,
  15 s on 16^3 x 32, 40 s on 20^3 x 40 (SU(3) about 2-3 times slower). The largest tasks may need
  more than one 2-day allocation; requeue or resubmission continues from the checkpoint.

## Collect and analyse

Raw `.npz` files land in `results/fasrc/` (not committed: they are large). Then:

```bash
cd lattice
python scripts/glueball_summary.py ../results/fasrc --md ../results/glueball-su2-fasrc.md --png ../results/glueball-su2-fasrc.png
```

Commit the summary (`.md`, `.png`) and the job logs' parameter lines, not the raw data. Independent
chains with the same coupling, lattice and smearing levels are pooled automatically.

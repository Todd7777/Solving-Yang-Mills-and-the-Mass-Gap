# Referee (Section 16, sixth pass): the capacity prune of order_counts.count_junctions, tested against the same search without it.
# Usage: python3 prune_ab.py MODE D K1:NMAX1 [K2:NMAX2 ...]     MODE in {both, on, off};  K as 4,4,6;  e.g. python3 prune_ab.py both 2 4,4:6 4,8:4
"""Runs the authors' count_junctions (../order_counts.py) for each class (sizes K of the further fronts, n <= NMAX) with PRUNE = True and/or
PRUNE = False and prints the counts by n and the times.  The two settings must give identical counts, since the prune is claimed to be a
necessary condition.  The search code is the authors'; only the flag is toggled."""
import os, sys, time
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
import order_counts as OC

def run(D, K, nmax, prune):
    OC.PRUNE = prune
    t0 = time.time()
    c = OC.count_junctions(D, nmax, K)
    return dict(sorted(c.items())), time.time() - t0

if __name__ == "__main__":
    mode, D = sys.argv[1], int(sys.argv[2])
    for spec in sys.argv[3:]:
        ks, nm = spec.split(':'); K = tuple(sorted(int(x) for x in ks.split(','))); nmax = int(nm)
        res = {}
        for prune in ([True, False] if mode == 'both' else [mode == 'on']):
            res[prune] = run(D, K, nmax, prune)
            print(f"D={D} K={K} n<={nmax} PRUNE={prune}: {res[prune][0]}  ({res[prune][1]:.1f} s)"); sys.stdout.flush()
        if mode == 'both':
            print(f"   -> {'IDENTICAL' if res[True][0] == res[False][0] else 'MISMATCH'}"); sys.stdout.flush()

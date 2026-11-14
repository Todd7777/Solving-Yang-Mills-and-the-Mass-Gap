# Referee probe (eleventh pass, Corollary mg:cor:classic; a SUGGESTION, not the corollary's setting): the authors' scaling_check.py in d = 3 with the fronts beyond
# the exact counts bounded by Lemma mg:lem:frontsdecision for Z^2 (third_order_effective.NV_DECISION = True; practical_beta.load_terms sets it False in d = 3) and
# fronts summed to 100 edges (abar = 2^-6 then). Nothing in the repository is modified: load_terms is wrapped in this process only.
# Usage: python3 scaling_decision_d3.py GROUP [JOBS]     (GROUP = 2 for SU(2), 1 for U(1))
import os, sys, runpy
here = os.path.dirname(os.path.abspath(__file__)); parent = os.path.join(here, '..'); sys.path.insert(0, parent)
G = sys.argv[1]; J = sys.argv[2] if len(sys.argv) > 2 else '2'
import practical_beta as PB, third_order_effective as T
_orig = PB.load_terms
def load_terms(d, *a, **k):
    _orig(d, *a, **k)
    if d == 3: T.NV_DECISION = True
PB.load_terms = load_terms
os.chdir(parent)
sys.argv = ['scaling_check.py', '--bounds', '--certified', '--group', G, '--kmax', '100', '--d3', 'orders/d3_s9z.json',
            '0.00040262143352281974401285396895067832618294529962141', J]
runpy.run_path(os.path.join(parent, 'scaling_check.py'), run_name='__main__')

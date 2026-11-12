# Written by the independent referee of Section 13 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs. Code written from the definitions alone (no lemma of the paper is used).
import time, sys
from collections import Counter
from lattice import *
for n, K in [(2, 8), (3, 7), (3, 8)]:
    t = time.time()
    F = fronts_through(tuple([0]*n), K, n)
    assert all(is_front(f) for f in F)
    c = Counter(len(f) for f in F)
    # same with a larger ball, to see that the ball does not cut anything
    t2 = time.time()
    print(f"n={n} K={K}: {dict(sorted(c.items()))}  ({t2-t:.1f}s)")
    sys.stdout.flush()
    if K <= 7:
        F2 = fronts_through(tuple([0]*n), K, n, R=K//2+3)
        print("   larger ball identical:", F2 == F)

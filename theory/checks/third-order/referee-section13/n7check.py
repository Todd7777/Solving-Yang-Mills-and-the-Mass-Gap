# Written by the independent referee of Section 13 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs. Code written from the definitions alone (no lemma of the paper is used).
import time
from collections import Counter
from lattice import *
for n in (2, 3, 4):
    t = time.time(); F = fronts_through(tuple([0]*n), 7, n)
    c = Counter(len(f) for f in F); d = n + 1; Pv = 2*n*(n-1); KS = 4*(2*n-3)
    print(f"d={d}: N4={c[4]} (P_v={Pv}), N6={c[6]} (12C(d-1,2)+96C(d-1,3)={12*n*(n-1)//2+96*n*(n-1)*(n-2)//6}), N7={c[7]} ((3/4)P_v K_S={3*Pv*KS//4})  ({time.time()-t:.0f}s)", flush=True)

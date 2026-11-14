# Referee check (tenth pass): the authors' sector_upper.E_sector (mpmath, their tails upper_exact.tails) against the referee's Arb formula
# (upper_arb_sectors.E_sector) at the same point beta_0(d), to many digits. Usage: python3 sector_upper_compare.py
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..'))
for sub in ('second-order', 'upper-bound'): sys.path.insert(0, os.path.join(here, '..', '..', sub))
from mpmath import mp, mpf, binomial, nstr
from flint import arb
import sector_upper as SU, upper_exact as UE, upper_arb_sectors as UA
from constants_upper import beta0
mp.dps = 60
worst = 0
for d in (3, 4):
    b0 = beta0(d); K = 4*(2*d - 3)
    if d == 3: cn = {3: UE.counts_to(3, 12)}
    else:
        cn = {4: {'Q': ({**{n: c for n, c in UE.COUNTS[4]['Q'][0].items() if n <= 10}, 11: binomial(48 + 11*K, 11)}, 18),
                  'A': ({**{n: c for n, c in UE.COUNTS[4]['A'][0].items() if n <= 10}, 11: binomial(20 + 11*K, 11)}, 4)}}
    bA = arb(nstr(b0, 70))
    for N in (3, 4, 5):
        for C in ('+', '-'):
            ea = SU.E_sector(b0, d, N, 0, C, cn)/b0**4
            eb, _ = UA.E_sector(bA, d, N, C); eb = eb/bA**4
            rel = abs(mpf(eb.mid().str(40, radius=False)) - ea)/ea
            worst = max(worst, rel)
            print(f"SU({N}) d={d} C={C}: authors {nstr(ea, 15)}, referee (Arb) {eb.mid().str(15, radius=False)}, relative difference {nstr(rel, 3)}")
print("largest relative difference:", nstr(worst, 3))

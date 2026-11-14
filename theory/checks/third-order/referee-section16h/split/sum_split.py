import re, glob, os, sys
from collections import Counter
here = os.path.dirname(os.path.abspath(__file__))
def parse(f):
    t = open(f).read()
    m = re.search(r'admissible sets by size:(.*)', t); m2 = re.search(r'of which with two components:(.*)', t); bh = re.search(r'boundary_hit=(\d)', t)
    a = Counter({int(x.split(':')[0]): int(x.split(':')[1]) for x in m.group(1).split()})
    b = Counter({int(x.split(':')[0]): int(x.split(':')[1]) for x in m2.group(1).split()})
    return a, b, int(bh.group(1))
for NM in (8, 9, 10):
    for K in 'AQ':
        f = os.path.join(here, f'un_{K}_{NM}.txt')
        if not os.path.exists(f): continue
        un, un2, bh = parse(f)
        line = f"Z^4 {K} nmax={NM}: unsplit {dict(sorted(un.items()))} (two comp. {dict(sorted(un2.items()))}, boundary_hit {bh})"
        for J in sorted({int(x.split('_')[-1][:-4]) for x in glob.glob(os.path.join(here, f'pc_{K}_{NM}_*_*.txt'))}):
            pcs = [os.path.join(here, f'pc_{K}_{NM}_{k}_{J}.txt') for k in range(J)]
            if not all(os.path.exists(p) and 'admissible sets by size' in open(p).read() for p in pcs): line += f"; J={J}: incomplete"; continue
            s, s2, hit = Counter(), Counter(), 0
            for p in pcs:
                a, b, h = parse(p); s += a; s2 += b; hit |= h
            line += f"; J={J}: sum of pieces equal: {s == un and s2 == un2}"
        print(line)

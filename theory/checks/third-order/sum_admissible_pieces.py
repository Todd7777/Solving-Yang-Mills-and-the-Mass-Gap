# Sums the pieces of run_admissible_pieces.sh: the admissible sets by size for (S_A,{p0}) and (S_Q,Q) in Z^4, and checks that every piece is present.
import sys, re, glob, os
out = sys.argv[1]
for kind, J in (('A', 4), ('Q', 16)):
    tot, tot2, missing = {}, {}, []
    for k in range(J):
        f = os.path.join(out, f'{kind}_{k}_{J}.txt')
        if not os.path.exists(f): missing.append(k); continue
        txt = open(f).read()
        m = re.search(r'admissible sets by size:(.*)', txt); m2 = re.search(r'of which with two components:(.*)', txt)
        if not m: missing.append(k); continue
        for part in m.group(1).split(): n, c = map(int, part.split(':')); tot[n] = tot.get(n, 0) + c
        for part in m2.group(1).split(): n, c = map(int, part.split(':')); tot2[n] = tot2.get(n, 0) + c
    print(f"{kind}: {J - len(missing)}/{J} pieces; admissible sets by size {dict(sorted(tot.items()))}; two components {dict(sorted(tot2.items()))}"
          + (f"; MISSING pieces {missing}" if missing else ""))

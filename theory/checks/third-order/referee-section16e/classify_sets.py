# Referee check (Remark mg:rem:upperexact (i)): the shapes of the admissible sets found by admissible_partition.py.
# Checks: every 5-set is a cap, every 6-set a unit cube sharing a link with p0 (an edge with Q) and containing no plaquette of C, every 9-set the
# boundary of a 2x1x1 box minus one face F in C (with exactly one face of the box in C), every two-component 10-set two caps; also lists the d=4 8-sets.
# Usage: python3 classify_sets.py d nmax
import sys
from collections import Counter
import admissible_partition as AP

def cube_faces(base, dirs, add, E):
    dirs = sorted(dirs); out = set()
    for i in range(3):
        for j in range(i + 1, 3):
            k = 3 - i - j
            out.add((base, dirs[i], dirs[j])); out.add((add(base, E[dirs[k]]), dirs[i], dirs[j]))
    return out

def all_cubes_through(p, d, add, E):
    x, a, b = p; res = []
    for c in range(d):
        if c in (a, b): continue
        for base in (x, tuple(xi - (1 if i == c else 0) for i, xi in enumerate(x))):
            res.append(frozenset(cube_faces(base, (a, b, c), add, E)))
    return res

if __name__ == "__main__":
    d, nmax = int(sys.argv[1]), int(sys.argv[2])
    for kind in ('A', 'Q'):
        conn, allsets, nodes, geo = AP.all_admissible(d, kind, nmax)
        links, plaqs, C, S, PS, E, add = geo
        bysize = Counter(); report = Counter()
        for s in allsets:
            Y = frozenset(p for c in s for p in c); n = len(Y)
            if len(s) == 2:
                report[(n, 'two components: two caps' if all(AP.is_cap(c, geo) for c in s) else 'two components: OTHER')] += 1
                continue
            if n == 5:
                report[(n, 'cap' if AP.is_cap(Y, geo) else 'OTHER')] += 1
            elif n == 6:
                cube = AP.is_cube(Y, geo)
                shares = bool(set(l for p in Y for l in links(p)) & S)
                report[(n, 'unit cube sharing a link with S, no plaquette of C' if cube and shares and not (Y & C) else 'OTHER')] += 1
            elif n == 9:
                ok = False
                for F in C:
                    Z = Y | {F}
                    if len(Z) != 10: continue
                    # Z is the boundary of a box Gamma1 u Gamma2 (two unit cubes sharing a face) iff Z = faces(G1) xor faces(G2) for cubes through F...
                    for G1 in all_cubes_through(F, d, add, E):
                        for P in G1:
                            for G2 in all_cubes_through(P, d, add, E):
                                if G2 == G1: continue
                                if (G1 ^ G2) == Z and len(Z & C) == 1:
                                    ok = True
                report[(n, 'box 2x1x1 minus its one face in C' if ok else 'OTHER')] += 1
            elif n in (7, 8):
                # in d = 4, (S_Q, Q): describe as (faces of G1 xor faces of G2) minus two faces of C?
                ok = False
                FF = [F for F in C]
                for i, F1 in enumerate(FF):
                    for F2 in FF[i + 1:]:
                        Z = Y | {F1, F2}
                        if len(Z) != n + 2: continue
                        for G1 in all_cubes_through(F1, d, add, E):
                            for G2 in all_cubes_through(F2, d, add, E):
                                if G1 != G2 and len(G1 & G2) == 1 and (G1 ^ G2) == Z:
                                    ok = True
                report[(n, 'two cubes on adjacent faces of Q sharing a face, minus the two faces of Q' if ok else 'OTHER')] += 1
            else:
                report[(n, 'not classified')] += 1
        print(f"d={d} (S,C)=({kind}):")
        for k in sorted(report): print(f"   n={k[0]:2d}: {report[k]:5d}  {k[1]}")

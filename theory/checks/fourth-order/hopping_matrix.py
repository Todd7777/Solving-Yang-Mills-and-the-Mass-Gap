# Fourth order (Part IV, Section 12): the hopping matrices H_+ (C-even) and H_- (C-odd) of the plaquette state on the spatial
# lattice Z^n (n = d-1), from their definition by orientations (no formula of the paper is assumed):
#   H_C(q,q) = K_S (number of plaquettes sharing a link with q: the bulges);
#   H_C(q,c) += s for c sharing a link with q (the hops), s = 1 (C = +) or sigma(q,c) (C = -), sigma(q,c) the sign with
#     bd(q) + sigma bd(c) of six edges;
#   H_C(q,q') += s for each unit cube with faces q != q' (turns and jumps), s = 1 (C = +) or sigma_G(q,q') = -eps_G(q) eps_G(q'),
#     eps_G(p) = +1 if the standard orientation of p is the outward one.
# Checked: (i) coplanar hops and jumps have sign +1, and on adjacent faces sigma_G = -sigma (turn and hop cancel for C = -);
# (ii) row sums h_+ = 26d - 70, h_- = 10d - 22, and these are the largest eigenvalues on tori;
# (iii) the k = 0 block over the n(n-1)/2 orientations: eigenvalues (A_1, E, T_1 in d = 4);
# (iv) the Rayleigh quotient of the uniform vector on the plaquettes of a box of side L tends to h_C (the trial sets of the theorem).
import itertools, sys
import numpy as np
sys.path.insert(0, '.')
from zn_block_spectra import geometry, chain_add, sign_to

def build(n, L, periodic=True):
    nl, plaqs, cubes = geometry((L,)*n, periodic)
    keys = []
    rng = [range(L) if periodic else range(L + 1)]*n
    for x in itertools.product(*rng):
        for i, j in itertools.combinations(range(n), 2):
            if periodic or (x[i] < L and x[j] < L): keys.append((x, i, j))
    if len(keys) != len(plaqs):        # open blocks: rebuild keys in the order of geometry()
        keys = None
    P = len(plaqs)
    Hp = np.zeros((P, P)); Hm = np.zeros((P, P)); stats = {'coplanar hop signs': set(), 'jump signs': set(), 'turn+hop': set()}
    hop = {}
    for q in range(P):
        for c in range(P):
            if c != q and set(plaqs[q]) & set(plaqs[c]):
                s = [a for a in (1, -1) if len(chain_add((1, plaqs[q]), (a, plaqs[c]))) == 6][0]
                hop[(q, c)] = s
                Hp[q, q] += 1; Hm[q, q] += 1; Hp[c, q] += 1; Hm[c, q] += s
                if keys and keys[q][1:] == keys[c][1:]: stats['coplanar hop signs'].add(s)
    for fs in cubes:
        for q in fs:
            for q2 in fs:
                if q2 == q: continue
                S = [f for f in fs if f not in (q, q2)]
                sg = sign_to(plaqs, q, S, q2)
                Hp[q2, q] += 1; Hm[q2, q] += sg
                if (q, q2) in hop: stats['turn+hop'].add(sg + hop[(q, q2)])
                else: stats['jump signs'].add(sg)
    return Hp, Hm, keys, stats

def k0_block(H, keys, n):
    planes = list(itertools.combinations(range(n), 2))
    idx = {pl: k for k, pl in enumerate(planes)}
    B = np.zeros((len(planes), len(planes)))
    q0 = {pl: [i for i, kk in enumerate(keys) if kk[1:] == pl and all(c == 0 for c in kk[0])][0] for pl in planes}
    for pl in planes:
        col = H[:, q0[pl]]
        for i, kk in enumerate(keys):
            B[idx[kk[1:]], idx[pl]] += col[i]
    return np.round(np.linalg.eigvalsh((B + B.T)/2), 6)

if __name__ == "__main__":
    for n in (2, 3, 4):
        d = n + 1; L = 4 if n < 4 else 3
        Hp, Hm, keys, stats = build(n, L)
        rs = (set(np.round(Hp.sum(axis=0), 9)), set(np.round(Hm.sum(axis=0), 9)))
        print(f"d={d} (torus {L}^{n}): {dict((k, sorted(v)) for k, v in stats.items())}")
        print(f"   row sums: H_+ {sorted(rs[0])} (26d-70 = {26*d-70}), H_- {sorted(rs[1])} (10d-22 = {10*d-22});"
              f" top eigenvalues {np.linalg.eigvalsh(Hp)[-1]:.6f}, {np.linalg.eigvalsh(Hm)[-1]:.6f}")
        print(f"   k=0 orientation blocks: H_+ {k0_block(Hp, keys, n).tolist()}, H_- {k0_block(Hm, keys, n).tolist()}")
    # (iv) uniform trial vectors on boxes inside a large torus (n = 3, d = 4): plaquettes with all vertices in [0, Lb]^3
    n, Lt = 3, 12
    Hp, Hm, keys, _ = build(n, Lt)
    for Lb in (1, 2, 3, 4, 6, 8):
        box = [i for i, (x, a, b) in enumerate(keys) if all(c + (1 if k in (a, b) else 0) <= Lb for k, c in enumerate(x))]
        for name, H, h in (('+', Hp, 34), ('-', Hm, 18)):
            sub = H[np.ix_(box, box)]
            rq = sub.sum()/len(box); top = np.linalg.eigvalsh(sub)[-1]
            print(f"d=4, box of side {Lb} ({len(box)} plaquettes), C{name}: uniform Rayleigh quotient {rq:.4f}, top eigenvalue {top:.4f} (h = {h})")

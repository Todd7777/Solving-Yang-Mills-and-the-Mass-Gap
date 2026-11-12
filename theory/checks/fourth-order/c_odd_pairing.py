"""The hop and the turn between perpendicular neighbours, as two-slice blocks (a compass, not part of the proofs).

For a plaquette q and a neighbour c sharing the link l, the block of the slices s-1, s is
    K^{1/2}_{s-2} M_{s-1} K_{s-1} M_s K^{1/2}_s,
with the front dq below and dc above.
  hop:  c at slice s-1, the six temporal plaquettes over dq xor dc in layer s-1, q at slice s;
  turn: the four temporal plaquettes over dq in layer s-1 (the tube), the four faces of the cube other than q, c at slice s.
We compute the single-copy integrals J(a,b,S) = E[ conj chi_a(U_q at s-1) chi_b(U_c at s) prod_{p in S} f(U_p) ] by exact group sums
over all link variables of the local three-dimensional complex (temporal links included, gauge fixed on a spanning tree), for every
subset S of the eight plaquettes. From them we assemble the two-copy block
    B[(a,a'),(b,b')] = uh_{aa'}^2 uh_{bb'}^2 sum_{S1 u S2 = P} J(a,b,S1) J(a',b',S2)
and check:
  (1) only S = empty (a = b = 0) and S = P contribute, with J = u_a^8 [b = a^sigma] for the hop and [b = a^sigma_Gamma] for the turn;
  (2) B_hop = uh^12 on the map (a,a') -> (a^s, a'^s), and B_turn = B_hop o C (C: a -> conj a), so B_hop + B_turn vanishes on the
      C-odd sector, for both signs sigma(q,c) = -1, +1;
  (3) the coplanar hop has sigma = +1 (no turn exists there).
Usage: python3 c_odd_pairing.py
"""
import itertools
import numpy as np
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0] + '/referee-section14')
from groups import Zn  # noqa: E402

E = [np.array(v) for v in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]


def plaq(x, i, j):
    """Standard orientation of the plaquette at x spanned by axes i<j: x -> x+e_i -> x+e_i+e_j -> x+e_j."""
    x = np.array(x)
    return [tuple(x), tuple(x + E[i]), tuple(x + E[i] + E[j]), tuple(x + E[j])]


def cube_faces(x):
    return [plaq(np.array(x) + s * E[k], i, j) for k, (i, j) in ((2, (0, 1)), (1, (0, 2)), (0, (1, 2))) for s in (0, 1)]


def links_of(cyc):
    return [(cyc[t], cyc[(t + 1) % len(cyc)]) for t in range(len(cyc))]


def build(kind, q, c):
    """Return (spatial cycles at time s-1 and s, temporal plaquettes, insertion cycles) as cycles of (vertex, time)."""
    ql, cl = set(map(frozenset, links_of(q))), set(map(frozenset, links_of(c)))
    plaqs = []
    if kind == 'hop':
        plaqs.append([(v, 0) for v in c])            # c at slice s-1 (time 0)
        plaqs.append([(v, 1) for v in q])            # q at slice s (time 1)
        front = ql ^ cl                              # six-edge front in layer s-1
    else:
        Gam = cube_faces(tuple(np.min(np.array(q + c), axis=0)))
        rest = [f for f in Gam if set(map(frozenset, links_of(f))) not in (ql, cl)]
        assert len(rest) == 4
        plaqs += [[(v, 1) for v in f] for f in rest]  # four faces at slice s
        front = ql                                   # the tube in layer s-1
    for e in front:
        v, w = sorted(e)
        plaqs.append([(v, 0), (w, 0), (w, 1), (v, 1)])  # temporal plaquette over e in layer s-1
    ins_in = [(v, 0) for v in q]
    ins_out = [(v, 1) for v in c]
    return plaqs, ins_in, ins_out


def run(G, irho, beta, kind, q, c):
    plaqs, ins_in, ins_out = build(kind, q, c)
    cycles = plaqs + [ins_in, ins_out]
    links = sorted({tuple(sorted(l)) for cyc in cycles for l in links_of(cyc)})
    verts = sorted({v for l in links for v in l})
    # spanning tree
    tree, seen, changed = [], {verts[0]}, True
    while changed:
        changed = False
        for (v, w) in links:
            if (v in seen) != (w in seen):
                tree.append((v, w)); seen |= {v, w}; changed = True
    assert len(seen) == len(verts)
    free = [l for l in links if l not in tree]
    n = G.n
    grids = np.meshgrid(*[np.arange(n)] * len(free), indexing='ij')
    U = {l: g.ravel() for l, g in zip(free, grids)}
    size = n ** len(free)
    for l in tree:
        U[l] = np.zeros(size, dtype=np.int64)

    def hol(cyc):
        g = None
        for (v, w) in links_of(cyc):
            x = U[(v, w)] if (v, w) in U else G.inv[U[(w, v)]]
            g = x if g is None else G.mul[g, x]
        return g

    f, u = G.setup(irho, beta)
    fP = [f[hol(p)] for p in plaqs]
    hq, hc = hol(ins_in), hol(ins_out)
    idx = hq * n + hc
    R = G.R
    J = np.zeros((2 ** len(plaqs), R, R), dtype=complex)
    for S in range(2 ** len(plaqs)):
        w = np.ones(size)
        for i in range(len(plaqs)):
            if S >> i & 1:
                w = w * fP[i]
        W = np.bincount(idx, weights=w, minlength=n * n).reshape(n, n) / size
        J[S] = G.chars.conj() @ W @ G.chars.T  # J[a,b] = sum_{x,y} conj chi_a(x) chi_b(y) W(x,y)
    return J, u, len(plaqs)


def two_copy(J, u, nP):
    R = J.shape[1]
    full = 2 ** nP - 1
    uh = np.array([[0.0 if a == 0 and b == 0 else (u[b] if a == 0 else (u[a] if b == 0 else u[a] * u[b]))
                    for b in range(R)] for a in range(R)])
    B = np.zeros((R, R, R, R), dtype=complex)
    for S1 in range(2 ** nP):
        for S2 in range(2 ** nP):
            if S1 | S2 == full:
                B += np.einsum('ab,cd->acbd', J[S1], J[S2])
    B *= (uh ** 2)[:, :, None, None] * (uh ** 2)[None, None, :, :]
    return B.reshape(R * R, R * R), uh


def check(G, irho, beta, q, c, label, with_turn=True):
    out = {}
    for kind in (('hop', 'turn') if with_turn else ('hop',)):
        J, u, nP = run(G, irho, beta, kind, q, c)
        full = 2 ** nP - 1
        # (1) only S = empty and S = P
        mid = max(np.abs(J[S]).max() for S in range(1, full))
        e0 = np.abs(J[0] - np.diag([1.0] + [0.0] * (G.R - 1))).max()
        perm = np.argmax(np.abs(J[full][1:, 1:]), axis=1) + 1
        ok = np.allclose(np.abs(J[full][np.arange(1, G.R), perm]), u[1:] ** 8, rtol=1e-9, atol=1e-15)
        twist = 'id' if np.all(perm == np.arange(1, G.R)) else ('conj' if np.all(perm == G.conjperm[1:]) else 'other')
        B, uh = two_copy(J, u, nP)
        out[kind] = (B, twist)
        print(f'  {label} {kind}: max|J(S)| over 0<S<P = {mid:.1e}; |J(empty)-delta| = {e0:.1e}; '
              f'J(P) = u^8 on b = {twist}(a): {ok}')
    R = G.R
    Cp = np.zeros((R * R, R * R))
    for a in range(R):
        for b in range(R):
            Cp[a * R + b, G.conjperm[a] * R + G.conjperm[b]] = 1
    Qm = (np.eye(R * R) - Cp) / 2
    Bh = out['hop'][0]
    tw = G.conjperm if out['hop'][1] == 'conj' else np.arange(R)
    J, u, nP = run(G, irho, beta, 'hop', q, c)
    uh = np.array([[0.0 if a == 0 and b == 0 else (u[b] if a == 0 else (u[a] if b == 0 else u[a] * u[b]))
                    for b in range(R)] for a in range(R)])
    pred = np.zeros((R * R, R * R))
    for a in range(R):
        for b in range(R):
            if a or b:
                pred[a * R + b, tw[a] * R + tw[b]] = uh[a, b] ** 12
    print(f'  {label}: |B_hop - uh^12 on the twist|/|B_hop| = {np.abs(Bh - pred).max() / np.abs(Bh).max():.1e}')
    if with_turn:
        Bt = out['turn'][0]
        r1 = np.abs(Bt - Bh @ Cp).max() / np.abs(Bh).max()
        r2 = np.abs(Qm @ (Bh + Bt) @ Qm).max() / np.abs(Bh).max()
        Qp = np.eye(R * R) - Qm
        r3 = np.linalg.norm(Qp @ (Bh + Bt) @ Qp, 2) / np.linalg.norm(Bh, 2)
        print(f'  {label}: |B_turn - B_hop C|/|B_hop| = {r1:.1e};  C-odd part of B_hop+B_turn: {r2:.1e};  '
              f'norm of the C-even part / norm of B_hop: {r3:.6f} (expect 2)')
        return r1, r2
    return None


if __name__ == '__main__':
    beta = 0.3
    q = plaq((0, 0, 0), 0, 1)
    c_minus = plaq((0, 0, 0), 0, 2)    # shares l = {(0,0,0),(1,0,0)}, traversed the same way: sigma = -1
    c_plus = plaq((0, 0, -1), 0, 2)    # shares l from below, traversed oppositely: sigma = +1
    c_flat = plaq((0, -1, 0), 0, 1)    # coplanar neighbour: sigma = +1, no common cube
    worst = 0.0
    for n in (3, 4, 5):
        G = Zn(n)
        print(f'Z{n}, rho = 1, beta = {beta}')
        for c, lab in ((c_minus, 'sigma=-1'), (c_plus, 'sigma=+1')):
            r1, r2 = check(G, 1, beta, q, c, lab)
            worst = max(worst, r1, r2)
        check(G, 1, beta, q, c_flat, 'coplanar', with_turn=False)
    print(f'worst relative residual of the identities: {worst:.1e}')
    assert worst < 1e-10

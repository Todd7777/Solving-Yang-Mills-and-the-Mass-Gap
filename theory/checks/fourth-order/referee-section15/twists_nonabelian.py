# Written by the independent referee of Section 15 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs.
# Referee of Section 15 (2026-11-12). A compass, not part of the proofs.
# Lemma 15.4 for non-abelian groups with non-self-conjugate irreducible representations (F21 = Z7 x| Z3, with its 3-dimensional complex
# irreps; A4 with its complex 1-dimensional ones), from the lattice geometry: holonomies are products of link variables along oriented
# cycles (the hexagon is the cycle of the 1-chain dq + sigma dc, sigma of Definition 12.1), integrals are exact group sums over all links of
# one slice (gauge fixed on a spanning tree). Single-slice junction integrals I(a,b,S) (Lemma 14.1) for the widening (q -> hexagon, plaquette c),
# the narrowing (hexagon -> c, plaquette q) and the cube junction of the turn (q -> c, the four faces). Checks:
#  - widening and narrowing are diagonal in the irreps (twist id) for the hexagon oriented as d(q + sigma c); so the hop has the twist a -> a^sigma;
#  - the turn has the twist a -> a^{sigma_Gamma} with sigma_Gamma(q,c) = -sigma(q,c), eps_Gamma from d(dGamma)=0;
#  - assembled two-copy blocks: Y_hop = sum uhat^12 |e^q_aa'><e^c_{a^s a'^s}| and Y0_turn = Y_hop C; Q-(Y_hop + Y0_turn)Q- = 0.
import sys, itertools
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 2)[0] + '/referee-section14')
from groups import F21, A4, Zn
from geom import add, unit, boundary, chain_add, sigma, cube_faces, eps_gamma

def cycle_order(chain):
    """Ordered traversal [(link, sign)] of an oriented 1-cycle given as dict link -> +-1 (a single cycle)."""
    steps = []
    for l, s in chain.items():
        b, d = l; h = add(b, unit(d))
        steps.append((l, s, b if s > 0 else h, h if s > 0 else b))
    order = [steps.pop(0)]
    while steps:
        k = [i for i, st in enumerate(steps) if st[2] == order[-1][3]][0]
        order.append(steps.pop(k))
    assert order[-1][3] == order[0][2]
    return [(l, s) for l, s, _, _ in order]

def setup_links(cycles):
    links = sorted({l for cyc in cycles for l in cyc})
    vs = sorted({v for (b, d) in links for v in (b, add(b, unit(d)))})
    tree, seen, ch = [], {vs[0]}, True
    while ch:
        ch = False
        for (b, d) in links:
            h = add(b, unit(d))
            if (b in seen) != (h in seen):
                tree.append((b, d)); seen |= {b, h}; ch = True
    assert len(seen) == len(vs)
    return links, [l for l in links if l not in tree]

def integrals(G, f, cyc_in, cyc_out, plaqs):
    """I[S][a,b] = E[ conj chi_a(U_in) chi_b(U_out) prod_{p in S} f(U_p) ] for all subsets S of plaqs (single copy)."""
    cycles = [cyc_in, cyc_out] + [cycle_order(boundary(p)) for p in plaqs]
    links, free = setup_links([[l for l, _ in cyc] for cyc in cycles])
    nf = len(free)
    grids = np.indices((G.n,) * nf).reshape(nf, -1) if nf else np.zeros((0, 1), dtype=int)
    U = {l: grids[free.index(l)] if l in free else np.zeros(grids.shape[1], dtype=int) for l in links}
    def hol(cyc):
        g = None
        for l, s in cyc:
            x = U[l] if s > 0 else G.inv[U[l]]
            g = x if g is None else G.mul[g, x]
        return g
    hin, hout = hol(cyc_in), hol(cyc_out)
    fP = [f[hol(cycle_order(boundary(p)))] for p in plaqs]
    out = {}
    for S in range(2 ** len(plaqs)):
        w = np.ones(grids.shape[1])
        for i in range(len(plaqs)):
            if S >> i & 1:
                w = w * fP[i]
        out[S] = (G.chars[:, hin].conj() * w) @ G.chars[:, hout].T / grids.shape[1]
    return out

def two_copy(I, nP, R):
    full = 2 ** nP - 1
    B = np.zeros((R, R, R, R), dtype=complex)
    for S1 in range(2 ** nP):
        for S2 in range(2 ** nP):
            if S1 | S2 == full:
                B += np.einsum('ab,cd->acbd', I[S1], I[S2])
    return B.reshape(R * R, R * R)

def run(G, irho, beta, q, c):
    f, u = G.setup(irho, beta)
    R = G.R
    uh = np.array([[0.0 if a == 0 and b == 0 else (u[b] if a == 0 else (u[a] if b == 0 else u[a] * u[b])) for b in range(R)] for a in range(R)]).ravel()
    sg = sigma(q, c)
    hexc = cycle_order(chain_add(boundary(q), boundary(c), sg))
    cq, cc = cycle_order(boundary(q)), cycle_order(boundary(c))
    Gam = [F for F in (cube_faces(add(q[0], d)) for d in itertools.product((-1, 0), repeat=3)) if q in F and c in F][0]
    eps = eps_gamma(Gam); sgG = -eps[q] * eps[c]
    rest = [p for p in Gam if p not in (q, c)]
    Iw = integrals(G, f, cq, hexc, [c]); In = integrals(G, f, hexc, cc, [q]); It = integrals(G, f, cq, cc, rest)
    perm = lambda M: [int(np.argmax(np.abs(M[a]))) for a in range(R)]
    ids, cj = list(range(R)), list(G.conjperm)
    live = [a for a in range(1, R) if u[a] > 1e-12]          # rows with u_a = 0 (e.g. the 3-dim irrep of A4 when f factors through Z3) carry no twist
    tw = lambda M: 'id' if all(perm(M)[a] == a for a in live) else ('conj' if all(perm(M)[a] == cj[a] for a in live) else 'other')
    mag = lambda M, k: max(abs(abs(M[a, perm(M)[a]]) - u[a] ** k) for a in live)
    off = lambda Is, full: max(np.abs(Is[S]).max() for S in Is if S not in (0, full))
    print(f'  {G.name} rho={irho} (dim {G.dims[irho]}, self-conj: {G.conjperm[irho] == irho}), q={q}, c={c}: sigma={sg:+d}, sigma_Gamma={sgG:+d}')
    print(f'    widening twist {tw(Iw[1])} (|I|-u_a: {mag(Iw[1],1):.1e}); narrowing twist {tw(In[1])} ({mag(In[1],1):.1e}); '
          f'turn twist {tw(It[15])} (|I|-u_a^4: {mag(It[15],4):.1e}; 0<S<P: {off(It,15):.1e})')
    hop_tw = 'id' if (tw(Iw[1]), tw(In[1])) in (('id', 'id'), ('conj', 'conj')) else 'conj'
    print(f'    hop twist = {hop_tw} (expect {"id" if sg > 0 else "conj"});  turn twist = {tw(It[15])} (expect {"id" if sgG > 0 else "conj"})')
    # two-copy blocks with the K factors: widening O1 = uh^2 [M_gc] uh^3, narrowing O2 = uh^3 [M_gq] uh^2, hexagon layer split in halves
    O1 = (uh ** 2)[:, None] * two_copy(Iw, 1, R) * (uh ** 3)[None, :]
    O2 = (uh ** 3)[:, None] * two_copy(In, 1, R) * (uh ** 2)[None, :]
    Yh = O1 @ O2
    Y0 = (uh ** 4)[:, None] * ((uh ** 2)[:, None] * two_copy(It, 4, R) * (uh ** 2)[None, :])   # plain tube on q (layer s-1), then the cube junction
    Cm = np.zeros((R * R, R * R))
    for a in range(R):
        for b in range(R):
            Cm[G.conjperm[a] * R + G.conjperm[b], a * R + b] = 1
    Qm = (np.eye(R * R) - Cm) / 2
    nrm = lambda X: np.linalg.norm(X, 2)
    print(f'    |Y_hop|/u_*^12 = {nrm(Yh)/max(u[1:])**12:.10f};  |Y0_turn - Y_hop C|/|Y_hop| = {nrm(Y0 - Yh @ Cm)/nrm(Yh):.1e};'
          f'  |Q-(Y_hop+Y0_turn)Q-|/|Y_hop| = {nrm(Qm @ (Yh + Y0) @ Qm)/nrm(Yh):.1e}')
    return nrm(Y0 - Yh @ Cm) / nrm(Yh)

if __name__ == '__main__':
    q = ((0, 0, 0), 0, 1)
    cs = [((0, 0, 0), 0, 2), ((0, 0, -1), 0, 2), ((0, 0, 0), 1, 2), ((1, 0, 0), 1, 2)]
    worst = 0
    for G, irho in ((F21(), None), (A4(), None), (Zn(5), 1)):
        if irho is None:
            irho = [r for r in range(G.R) if G.conjperm[r] != r][-1]
        for c in cs:
            worst = max(worst, run(G, irho, 0.8, q, c))
    print(f'worst |Y0_turn - Y_hop C|/|Y_hop| = {worst:.1e}')

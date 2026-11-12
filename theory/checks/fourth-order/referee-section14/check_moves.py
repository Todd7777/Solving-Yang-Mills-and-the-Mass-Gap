# Written by the independent referee of Section 14 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs.
# Independent check of Lemma 14.1 (moves by characters) on finite groups, from the lattice geometry:
#  - widening/narrowing on a FLAT and a BENT domino (q, c sharing a link), fronts dq and dq xor dc,
#  - turn and jump through a unit cube, fronts dq and dq', P = the four other faces.
# For every irreducible a, b and EVERY subset S of P we compute I(a,b,S) by exact group sums (gauge fixed on a spanning tree; for the
# domino also without gauge fixing for small groups), and then the two-copy matrix <e^-_{aa'}, M_J e^+_{bb'}> both from the
# S1 u S2 = P expansion with ALL subsets, and (domino) by a direct brute-force sum over both copies. Norms vs u_*^k.
import numpy as np, itertools, sys
from groups import Q8, A4, F21, S3, Zn

def lattice(vertices, faces_vertices):
    """vertices: list of int tuples; faces: dict name -> oriented vertex cycle. Links = unordered nearest-neighbour pairs used by faces."""
    links = set()
    for cyc in faces_vertices.values():
        for t in range(len(cyc)):
            v, w = cyc[t], cyc[(t + 1) % len(cyc)]
            links.add((min(v, w), max(v, w)))
    return sorted(links)

def spanning_tree(vertices, links):
    tree = []; seen = {vertices[0]}; changed = True
    while changed:
        changed = False
        for (v, w) in links:
            if (v in seen) != (w in seen):
                tree.append((v, w)); seen |= {v, w}; changed = True
    assert len(seen) == len(vertices)
    return tree

def holonomy(G, cyc, U):
    """U: dict link -> array of element indices (vectorised). Path-ordered product along the oriented cycle."""
    g = None
    for t in range(len(cyc)):
        v, w = cyc[t], cyc[(t + 1) % len(cyc)]
        if (v, w) in U: x = U[(v, w)]
        else: x = G.inv[U[(w, v)]]
        g = x if g is None else G.mul[g, x]
    return g

def configs(G, links, tree, gauge_fix=True):
    free = [l for l in links if not (gauge_fix and l in tree)]
    grids = np.meshgrid(*[np.arange(G.n)] * len(free), indexing='ij')
    U = {l: grid.ravel() for l, grid in zip(free, grids)}
    size = G.n ** len(free)
    if gauge_fix:
        for l in tree: U[l] = np.zeros(size, dtype=np.int64)
    return U, size

def I_table(G, f, cyc_minus, cyc_plus, Pcycles, U):
    """I[a,b,S] = mean conj(chi_a(U_-)) chi_b(U_+) prod_{p in S} f(U_p), S as bitmask over Pcycles."""
    hm = holonomy(G, cyc_minus, U); hp = holonomy(G, cyc_plus, U)
    fP = [f[holonomy(G, c, U)] for c in Pcycles]
    nP = len(Pcycles); R = G.R
    I = np.zeros((R, R, 2 ** nP), dtype=complex)
    Cm = G.chars[:, hm].conj(); Cp = G.chars[:, hp]
    for S in range(2 ** nP):
        w = np.ones(len(hm))
        for i in range(nP):
            if S >> i & 1: w = w * fP[i]
        I[:, :, S] = (Cm * w) @ Cp.T / len(hm)
    return I

def check_I(G, u, I, nP, tag):
    R = G.R; full = 2 ** nP - 1; ok = True; sigma = {}
    for a in range(R):
        for b in range(R):
            for S in range(2 ** nP):
                v = I[a, b, S]
                if a == 0 and b == 0 and S == 0: exp = 1.0
                elif S == full and a != 0 and abs(v) > 1e-9: exp = u[a] ** nP; sigma.setdefault(a, set()).add(b)
                else: exp = 0.0
                if abs(v - exp) > 1e-9:
                    ok = False; print('   MISMATCH', tag, a, b, S, v, exp)
    # sigma must be a bijection of the non-trivial irreps with sigma(a) in {a, conj(a)}
    sig_ok = all(len(sigma.get(a, ())) == 1 for a in range(1, R))
    sig = {a: next(iter(sigma[a])) for a in sigma}
    sig_ok &= sorted(sig.values()) == list(range(1, R)) and all(sig[a] in (a, G.conjperm[a]) for a in sig)
    return ok and sig_ok, sig

def uhat(u, a, a2):
    if a == 0 and a2 == 0: return 0.0
    if a2 == 0: return u[a]
    if a == 0: return u[a2]
    return u[a] * u[a2]

def two_copy_from_I(G, u, I, nP, km, kp):
    """Matrix of P^J O_J P^J on class functions: sum over ALL (S1,S2) with S1|S2 = P of I(a,b,S1) I(a',b',S2), times uhat factors."""
    R = G.R; full = 2 ** nP - 1
    pairs = [(S1, S2) for S1 in range(full + 1) for S2 in range(full + 1) if (S1 | S2) == full]
    M = np.zeros((R * R, R * R), dtype=complex); Mraw = np.zeros_like(M)
    for a, a2, b, b2 in itertools.product(range(R), repeat=4):
        raw = sum(I[a, b, S1] * I[a2, b2, S2] for S1, S2 in pairs)
        Mraw[a * R + a2, b * R + b2] = raw
        M[a * R + a2, b * R + b2] = uhat(u, a, a2) ** (km / 2) * raw * uhat(u, b, b2) ** (kp / 2)
    return M, Mraw

def domino(kind):
    if kind == 'flat':
        q = [(0, 0), (1, 0), (1, 1), (0, 1)]; c = [(1, 0), (2, 0), (2, 1), (1, 1)]
        hexa = [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (0, 1)]
    else:
        q = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]; c = [(0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)]
        hexa = [(1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 0), (0, 0, 1), (1, 0, 1)]
    verts = sorted(set(q) | set(c))
    links = lattice(verts, {'q': q, 'c': c})
    # sanity: hexagon links = (dq xor dc)
    hl = {(min(hexa[t], hexa[(t+1) % 6]), max(hexa[t], hexa[(t+1) % 6])) for t in range(6)}
    ql = {(min(q[t], q[(t+1) % 4]), max(q[t], q[(t+1) % 4])) for t in range(4)}
    cl = {(min(c[t], c[(t+1) % 4]), max(c[t], c[(t+1) % 4])) for t in range(4)}
    assert hl == ql ^ cl and len(ql & cl) == 1
    return verts, links, q, c, hexa

def cube():
    V = list(itertools.product((0, 1), repeat=3))
    faces = {}
    for k in range(3):
        for s in (0, 1):
            i, j = [m for m in range(3) if m != k]
            def P(x, y): return tuple(s if m == k else (x if m == i else y) for m in range(3))
            faces[(k, s)] = [P(0, 0), P(1, 0), P(1, 1), P(0, 1)]
    links = lattice(V, faces)
    assert len(links) == 12
    return V, links, faces

def run(G, irho, beta, full_domino=False):
    f, u = G.setup(irho, beta)
    us = max(abs(u[1:]))
    out = []
    for kind in ('flat', 'bent'):
        verts, links, q, c, hexa = domino(kind)
        tree = spanning_tree(verts, links)
        for gf in ([True, False] if full_domino else [True]):
            U, size = configs(G, links, tree, gauge_fix=gf)
            for name, (cm, cp) in (('widening', (q, hexa)), ('narrowing', (hexa, q))):
                I = I_table(G, f, cm, cp, [c], U)
                ok, sig = check_I(G, u, I, 1, f'{kind} {name}')
                km, kp = len(cm), len(cp)
                M, Mraw = two_copy_from_I(G, u, I, 1, km, kp)
                nrm = np.linalg.svd(M, compute_uv=False).max()
                out.append((f'{kind:4s} domino {name:9s} {"gauge-fixed" if gf else "all links  "}', ok, sig, nrm, us ** 6))
        # direct two-copy brute force (gauge fixed in each copy) of <e^-_{aa'}, M_J e^+_{bb'}> for the widening
        U, size = configs(G, links, tree, True)
        hm = holonomy(G, q, U); hp = holonomy(G, hexa, U); hc = holonomy(G, c, U)
        fc = f[hc]
        # g_c(U,U') = (1+f(Uc))(1+f(U'c)) - 1 ; two independent copies -> outer products over config indices
        R = G.R
        Cm = G.chars[:, hm]; Cp = G.chars[:, hp]
        g = np.outer(1 + fc, 1 + fc) - 1                 # [cfg, cfg']
        direct = np.einsum('ax,ry,xy,bx,sy->arbs', Cm.conj(), Cm.conj(), g, Cp, Cp) / size ** 2
        I = I_table(G, f, q, hexa, [c], U)
        _, Mraw = two_copy_from_I(G, u, I, 1, 4, 6)
        err = abs(direct.reshape(R * R, R * R) - Mraw).max()
        out.append((f'{kind:4s} domino widening: direct two-copy sum vs expansion, max diff {err:.1e}', err < 1e-10, None, None, None))
    V, links, faces = cube()
    tree = spanning_tree(V, links)
    U, size = configs(G, links, tree, True)
    q = (2, 0)
    for kind, qp in (('turn', (0, 0)), ('turn2', (1, 1)), ('jump', (2, 1))):
        Pn = [F for F in faces if F not in (q, qp)]
        for orient in ('std', 'rev'):
            cq = faces[q]; cqp = faces[qp] if orient == 'std' else faces[qp][::-1]
            I = I_table(G, f, cq, cqp, [faces[F] for F in Pn], U)
            ok, sig = check_I(G, u, I, 4, f'cube {kind}')
            M, Mraw = two_copy_from_I(G, u, I, 4, 4, 4)
            nrm = np.linalg.svd(M, compute_uv=False).max()
            out.append((f'cube {kind:5s} (front orientation {orient})', ok, sig, nrm, us ** 8))
    return u, out

if __name__ == '__main__':
    for G, irho, full in ((Q8(), 4, True), (A4(), 3, True), (F21(), 3, False), (S3(), 2, True), (Zn(4), 1, True)):
        for beta in (0.4, 1.3):
            u, out = run(G, irho, beta, full)
            print(f'== {G.name} (rho dim {G.dims[irho]}), beta={beta}: u_r = {np.round(u, 6)}')
            for tag, ok, sig, nrm, ref in out:
                s = f'   {tag}: lemma values {"OK" if ok else "FAIL"}'
                if sig is not None: s += f', sigma={sig}'
                if nrm is not None: s += f', norm={nrm:.10e}, u_*^k={ref:.10e}, rel.diff={abs(nrm/ref-1):.1e}'
                print(s)

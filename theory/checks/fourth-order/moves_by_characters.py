# Lemma 14.1 of Part IV (the moves by characters), checked by brute force on finite groups (a compass, not part of the proof).
# For a widening (one plaquette c, fronts bd q and bd q xor bd c) and for a turn and a jump through a unit cube (fronts bd q, bd q'), we compute
#   I(a,b,S) = int conj(chi_a(U_-)) chi_b(U_+) prod_{p in S} f(U_p) dU    (single copy, exact group sums)
# for all irreducible a, b and all subsets S of the junction's plaquettes, and compare with the lemma: 1 if a=b=0, S=empty; u_a^|P| if S=P, b=sigma(a);
# 0 otherwise. Then we assemble the two-copy operator P^J O_J P^J on the class functions of the front holonomies and compare its norm with u_*^k,
# k = |P| + (|phi_-|+|phi_+|)/2 (6 for the widening, 8 for the cube moves).
# Groups: Z_3 (faithful character) and S_3 (two-dimensional irreducible representation). Usage: python3 moves_by_characters.py  (about a minute)
import itertools, numpy as np

def group_Zn(n):
    els = list(range(n)); mul = lambda x, y: (x + y) % n; inv = lambda x: (-x) % n
    chars = [lambda x, k=k: np.exp(2j*np.pi*k*x/n) for k in range(n)]
    return els, mul, inv, chars, [1]*n, 1          # rho = the faithful character k=1

def group_S3():
    els = list(itertools.permutations(range(3)))
    mul = lambda x, y: tuple(x[y[i]] for i in range(3)); inv = lambda x: tuple(sorted(range(3), key=lambda i: x[i]))
    def sign(p):
        s = 1
        for i in range(3):
            for j in range(i+1, 3):
                if p[i] > p[j]: s = -s
        return s
    fix = lambda p: sum(p[i] == i for i in range(3))
    chars = [lambda x: 1.0, lambda x: float(sign(x)), lambda x: float(fix(x) - 1)]
    return els, mul, inv, chars, [1, 1, 2], 2      # rho = the two-dimensional irrep

def setup(G, beta):
    els, mul, inv, chars, dims, irho = G
    Nrho = dims[irho]; w = {x: np.exp(beta/Nrho*np.real(chars[irho](x))) for x in els}
    c0 = np.mean(list(w.values())); f = {x: w[x]/c0 - 1 for x in els}
    u = [np.real(np.mean([f[x]*np.conj(chars[r](x)) for x in els]))/dims[r] for r in range(len(chars))]   # f = sum_r d_r u_r chi_r
    return els, mul, inv, chars, dims, f, u

def domino_I(G, beta):
    els, mul, inv, chars, dims, f, u = setup(G, beta)
    R = len(chars); out = {}
    # links: V (x->y, shared), and W_q, W_c (products of three links each; Haar-distributed, independent) -> U_q = V W_q, U_c = V W_c, U_+ ~ W_q W_c^-1
    for a, b in itertools.product(range(R), repeat=2):
        for S in ((), ('c',)):
            tot = 0
            for V, Wq, Wc in itertools.product(els, repeat=3):
                Uq = mul(V, Wq); Uc = mul(V, Wc); U6 = mul(Wq, inv(Wc))
                h = f[Uc] if S else 1.0
                tot += np.conj(chars[a](Uq))*chars[b](U6)*h
            out[(a, b, S)] = tot/len(els)**3
    return out, u, dims

def cube_I(G, beta, kind):
    """kind 'turn': fronts = bottom face q and a side face q'; 'jump': bottom q and top q'. P = the other four faces. Gauge fixed on a spanning tree."""
    els, mul, inv, chars, dims, f, u = setup(G, beta)
    R = len(chars)
    V = list(itertools.product((0, 1), repeat=3)); idx = {v: i for i, v in enumerate(V)}
    E = [(v, w) for v in V for w in V if v < w and sum(abs(v[k]-w[k]) for k in range(3)) == 1]
    faces = {}
    for k in range(3):
        for s in (0, 1):
            vs = sorted([v for v in V if v[k] == s]); i, j = [m for m in range(3) if m != k]
            a0 = tuple(s if m == k else 0 for m in range(3))
            def P(x, y, k=k, s=s, i=i, j=j): return tuple(s if m == k else (x if m == i else y) for m in range(3))
            faces[(k, s)] = [P(0, 0), P(1, 0), P(1, 1), P(0, 1)]                 # a cycle of four vertices
    q = (2, 0); qp = (0, 0) if kind == 'turn' else (2, 1)
    P_ = [F for F in faces if F not in (q, qp)]
    # spanning tree: BFS from (0,0,0); its links are set to the identity
    tree = set(); seen = {V[0]}; queue = [V[0]]
    while queue:
        v = queue.pop(0)
        for (x, y) in E:
            w = y if x == v else (x if y == v else None)
            if w is not None and w not in seen: seen.add(w); queue.append(w); tree.add((x, y))
    free = [e for e in E if e not in tree]; e_id = els[[i for i, x in enumerate(els) if all(mul(x, z) == z for z in els)][0]]
    def hol(cycle, U):
        g = e_id
        for t in range(4):
            v, w = cycle[t], cycle[(t+1) % 4]
            g = mul(g, U[(v, w)] if (v, w) in U else inv(U[(w, v)]))
        return g
    out = {}
    for vals in itertools.product(els, repeat=len(free)):
        U = {e: e_id for e in tree}; U.update(dict(zip(free, vals)))
        hq = hol(faces[q], U); hqp = hol(faces[qp], U); fP = [f[hol(faces[F], U)] for F in P_]
        for a, b in itertools.product(range(R), repeat=2):
            base = np.conj(chars[a](hq))*chars[b](hqp)
            for S in range(16):                                # all subsets S of the four plaquettes (bit i: the i-th plaquette of P)
                term = base*np.prod([fP[i] for i in range(4) if S >> i & 1])
                out[(a, b, S)] = out.get((a, b, S), 0) + term
    for k in out: out[k] /= len(els)**len(free)
    return out, u, dims

def two_copy_norm(Ifun, u, dims, k_minus, k_plus):
    R = len(u); uh = lambda a, a2: (0 if a == a2 == 0 else (u[a] if a2 == 0 else (u[a2] if a == 0 else u[a]*u[a2])))
    M = np.zeros((R*R, R*R), complex)
    for (a, a2), (b, b2) in itertools.product(itertools.product(range(R), repeat=2), repeat=2):
        tot = 0
        for S1, S2 in (('P', ''), ('', 'P'), ('P', 'P')):
            tot += Ifun(a, b, S1)*Ifun(a2, b2, S2)
        M[a*R+a2, b*R+b2] = uh(a, a2)**(k_minus/2)*tot*uh(b, b2)**(k_plus/2)      # (G_g^{1/2}) on the links of phi_-, resp. phi_+
    return np.linalg.svd(M, compute_uv=False).max(), max(abs(x) for x in u[1:])

for gname, G in (('Z_3', group_Zn(3)), ('S_3', group_S3())):
    for beta in (0.3, 1.0):
        I, u, dims = domino_I(G, beta)
        R = len(u); ok = True
        for a, b in itertools.product(range(R), repeat=2):
            exp0 = 1.0 if a == b == 0 else 0.0; expc = (u[a] if (a == b and a != 0) else 0.0)
            ok &= abs(I[(a, b, ())] - exp0) < 1e-10 and abs(I[(a, b, ('c',))] - expc) < 1e-10
        Ifun = lambda a, b, S: I[(a, b, ('c',) if S == 'P' else ())]
        nrm, us = two_copy_norm(Ifun, u, dims, 4, 6)
        print(f"{gname}, beta={beta}: widening: I(a,b,S) as in the lemma: {ok};  norm {nrm:.6e}  u_*^6 = {us**6:.6e}")
        for kind in ('turn', 'jump'):
            I, u, dims = cube_I(G, beta, kind)
            ok = True; perm = {}
            for a, b in itertools.product(range(R), repeat=2):
                exp0 = 1.0 if a == b == 0 else 0.0
                ok &= abs(I[(a, b, 0)] - exp0) < 1e-10
                ok &= all(abs(I[(a, b, S)]) < 1e-10 for S in range(1, 15))     # the other subsets vanish
                v = I[(a, b, 15)]
                if abs(v) > 1e-10:
                    ok &= a != 0 and abs(v - u[a]**4) < 1e-10 and a not in perm; perm[a] = b
            ok &= sorted(perm) == list(range(1, R))
            Ifun = lambda a, b, S, I=I: I[(a, b, 15 if S == 'P' else 0)]
            nrm, us = two_copy_norm(Ifun, u, dims, 4, 4)
            print(f"{gname}, beta={beta}: {kind}: I(a,b,S) as in the lemma (sigma = {perm}): {ok};  norm {nrm:.6e}  u_*^8 = {us**8:.6e}")

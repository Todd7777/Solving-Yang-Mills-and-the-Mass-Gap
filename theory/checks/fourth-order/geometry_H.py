# geometry_H.py -- written by the independent referee of Section 12 of Part IV (2026-11-12), added to the repository with credit;
# unchanged apart from this header. A compass, not part of the proofs.
# Referee code (independent): oriented plaquettes, the signs sigma(q,c) and sigma_Gamma(q,q'), and the hopping matrices H_+-
# of Definition 12.1, built from scratch on a spatial lattice given as a list of vertices and oriented links.
# Conventions (Definition 12.1): standard orientation of the plaquette (x; i<j): x -> x+e_i -> x+e_i+e_j -> x+e_j -> x.
# sigma(q,c): bd q + sigma bd c has six links.  eps_Gamma(p) = +1 iff the standard orientation of p equals the boundary
# orientation induced by a fixed orientation of the 3-cube Gamma; sigma_Gamma(q,q') = -eps(q) eps(q').
# The boundary orientation of the 3-cube is obtained WITHOUT outward normals: we solve sum_p eps(p) bd p = 0 over eps in {+-1}^6
# (boundary of a boundary), which has exactly the two solutions +-eps.
import itertools
import numpy as np

def unit(n, i):
    return tuple(1 if k == i else 0 for k in range(n))

def vadd(x, y, L=None):
    z = tuple(a + b for a, b in zip(x, y))
    if L is not None:
        z = tuple(a % l for a, l in zip(z, L))
    return z

class Lattice:
    """n-dimensional block [0,L1]x...x[0,Ln] (open) or torus Z_L1 x ... x Z_Ln (periodic)."""
    def __init__(self, dims, periodic):
        self.n = n = len(dims); self.dims = dims; self.periodic = periodic
        Lw = dims if periodic else None
        rng = [range(L) if periodic else range(L + 1) for L in dims]
        self.verts = list(itertools.product(*rng))
        vset = set(self.verts)
        self.links = []          # (x, i): from x to x+e_i
        for x in self.verts:
            for i in range(n):
                y = vadd(x, unit(n, i), Lw)
                if y in vset:
                    self.links.append((x, i))
        self.lidx = {l: k for k, l in enumerate(self.links)}
        self.Lw = Lw
        self.plaqs = []          # (x, i, j), i<j
        self.bd = []             # dict link index -> +-1 (standard orientation)
        for x in self.verts:
            for i, j in itertools.combinations(range(n), 2):
                xi = vadd(x, unit(n, i), Lw); xj = vadd(x, unit(n, j), Lw)
                if xi in vset and xj in vset:
                    d = {}
                    for l, s in (((x, i), 1), ((xi, j), 1), ((xj, i), -1), ((x, j), -1)):
                        k = self.lidx[l]; d[k] = d.get(k, 0) + s
                    self.plaqs.append((x, i, j)); self.bd.append(d)
        self.pidx = {p: k for k, p in enumerate(self.plaqs)}
        self.cubes = []          # lists of 6 plaquette indices
        for x in self.verts:
            for i, j, k in itertools.combinations(range(n), 3):
                ok = all(vadd(x, unit(n, a), Lw) in vset for a in (i, j, k))
                if not ok: continue
                xi, xj, xk = (vadd(x, unit(n, a), Lw) for a in (i, j, k))
                faces = [(x, i, j), (xk, i, j), (x, i, k), (xj, i, k), (x, j, k), (xi, j, k)]
                self.cubes.append([self.pidx[f] for f in faces])

    def P(self): return len(self.plaqs)

def chain(*terms):
    out = {}
    for s, d in terms:
        for l, c in d.items():
            out[l] = out.get(l, 0) + s * c
    return {l: c for l, c in out.items() if c}

def cube_eps(lat, cube):
    sols = []
    for eps in itertools.product((1, -1), repeat=6):
        if not chain(*[(e, lat.bd[p]) for e, p in zip(eps, cube)]):
            sols.append(eps)
    assert len(sols) == 2 and all(a == -b for a, b in zip(*sols)), sols
    return dict(zip(cube, sols[0]))

def build_H(lat, block_KS=True):
    """returns H_+, H_-, and the separate pieces (KS diag, hop matrix A_signed, A_unsigned, cube matrix C_signed, C_unsigned)."""
    P = lat.P()
    KS = np.zeros(P); Ap = np.zeros((P, P)); Am = np.zeros((P, P)); Cp = np.zeros((P, P)); Cm = np.zeros((P, P))
    sig = {}
    for q in range(P):
        for c in range(P):
            if c == q: continue
            common = set(lat.bd[q]) & set(lat.bd[c])
            if not common: continue
            assert len(common) == 1
            s = [a for a in (1, -1) if len(chain((1, lat.bd[q]), (a, lat.bd[c]))) == 6]
            assert len(s) == 1
            sig[(q, c)] = s[0]
            KS[q] += 1; Ap[q, c] += 1; Am[q, c] += s[0]
    checks = {'sigmaG=-sigma (adjacent)': set(), 'sigmaG (opposite)': set(), 'sigma coplanar': set(), 'cubes per pair': set()}
    paircount = {}
    for cube in lat.cubes:
        eps = cube_eps(lat, cube)
        for q in cube:
            for q2 in cube:
                if q2 == q: continue
                sg = -eps[q] * eps[q2]
                Cp[q, q2] += 1; Cm[q, q2] += sg
                paircount[(q, q2)] = paircount.get((q, q2), 0) + 1
                if (q, q2) in sig: checks['sigmaG=-sigma (adjacent)'].add(sg == -sig[(q, q2)])
                else: checks['sigmaG (opposite)'].add(sg)
    for (q, c), s in sig.items():
        if lat.plaqs[q][1:] == lat.plaqs[c][1:]:
            checks['sigma coplanar'].add(s)
    checks['cubes per pair'] = set(paircount.values())
    Hp = np.diag(KS) + Ap + Cp
    Hm = np.diag(KS) + Am + Cm
    return Hp, Hm, dict(KS=KS, Ap=Ap, Am=Am, Cp=Cp, Cm=Cm, sig=sig), checks

if __name__ == '__main__':
    for dims, per in (((5, 5), True), ((4, 4, 4), True), ((5, 5, 5), True), ((1, 1, 1), False), ((3, 3), True), ((4, 4), True)):
        lat = Lattice(dims, per)
        Hp, Hm, parts, checks = build_H(lat)
        d = len(dims) + 1
        print(f"{'torus' if per else 'block'} {dims} (d={d}): P={lat.P()}, cubes={len(lat.cubes)}; checks: "
              + ', '.join(f'{k}: {sorted(v)}' for k, v in checks.items()))
        rsp = sorted(set(np.round(Hp.sum(1), 9))); rsm = sorted(set(np.round(Hm.sum(1), 9)))
        ep = np.sort(np.linalg.eigvalsh(Hp))[::-1]; em = np.sort(np.linalg.eigvalsh(Hm))[::-1]
        print(f"   row sums H+ {rsp}, H- {rsm}; 26d-70={26*d-70}, 10d-22={10*d-22}; top eigs {ep[0]:.6f}, {em[0]:.6f}")
        print(f"   K_S values {sorted(set(parts['KS']))} (4(2d-5) = {4*(2*d-5)}); min entries H+ {Hp.min():.0f}, H- {Hm.min():.0f};"
              f" H- has negative entries: {(Hm < -1e-12).any()}")
        if lat.P() <= 20:
            print('   spec H+ =', np.round(ep, 6).tolist()); print('   spec H- =', np.round(em, 6).tolist())
        # k = 0 block in d = 4 on tori: sum H over plaquettes of each orientation
        if per and len(dims) == 3:
            planes = sorted(set(p[1:] for p in lat.plaqs))
            for name, H in (('+', Hp), ('-', Hm)):
                B = np.zeros((3, 3))
                for a, pa in enumerate(planes):
                    q0 = lat.pidx[((0, 0, 0),) + pa]
                    for b, pb in enumerate(planes):
                        B[a, b] = sum(H[q0, q] for q in range(lat.P()) if lat.plaqs[q][1:] == pb)
                print(f'   k=0 block H{name}:', B.tolist(), ' eigenvalues', np.round(np.linalg.eigvalsh(B), 6).tolist())
            # second moments sum_q' H(q,q')^2 (relevant to the free-energy coefficient of N_t^2 u^8)
            print('   sum_q\' H+(q,q\')^2 =', sorted(set(np.round((Hp ** 2).sum(1), 6))), '; H-:', sorted(set(np.round((Hm ** 2).sum(1), 6))))

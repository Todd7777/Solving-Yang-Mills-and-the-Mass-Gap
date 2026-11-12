# Written by the independent referee of Section 15 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs.
# Referee of Section 15 (2026-11-12). A compass, not part of the proofs.
# Lattice geometry from the definitions: oriented links (base, dir), plaquettes (base, i, j) with i<j and the standard orientation
# x -> x+e_i -> x+e_i+e_j -> x+e_j -> x, boundary 1-chains, sigma(q,c) of Definition 12.1, and eps_Gamma from d(d Gamma)=0.
import itertools
import numpy as np

def add(x, y):
    return tuple(a + b for a, b in zip(x, y))

def unit(i, n=3):
    return tuple(1 if k == i else 0 for k in range(n))

def boundary(p):
    """Oriented boundary of the plaquette p=(x,i,j): dict link -> +-1, link = (base, dir) oriented base -> base+e_dir."""
    x, i, j = p
    return {(x, i): 1, (add(x, unit(i)), j): 1, (add(x, unit(j)), i): -1, (x, j): -1}

def chain_add(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + s * v
        if out[k] == 0:
            del out[k]
    return out

def sigma(q, c):
    """Definition 12.1: the sign for which dq + sigma dc has six links."""
    res = [s for s in (1, -1) if len(chain_add(boundary(q), boundary(c), s)) == 6]
    assert len(res) == 1
    return res[0]

def cube_faces(x):
    """The six faces of the unit 3-cube at x spanned by e1,e2,e3 (axes 0,1,2)."""
    out = []
    for (i, j), k in (((0, 1), 2), ((0, 2), 1), ((1, 2), 0)):
        out += [(x, i, j), (add(x, unit(k)), i, j)]
    return out

def eps_gamma(faces):
    """Signs eps(p) with sum eps(p) dp = 0 (the induced boundary orientation, up to a global sign)."""
    sols = []
    for signs in itertools.product((1, -1), repeat=len(faces)):
        tot = {}
        for s, p in zip(signs, faces):
            tot = chain_add(tot, boundary(p), s)
        if not tot:
            sols.append(signs)
    assert len(sols) == 2
    return dict(zip(faces, sols[0]))

def holonomy(p, U, links, n):
    """Holonomy of p (abelian Z_n, additive) as an array expression over link variables U[k] (k = index in links)."""
    h = 0
    for l, s in boundary(p).items():
        h = h + s * U[links.index(l)]
    return h % n

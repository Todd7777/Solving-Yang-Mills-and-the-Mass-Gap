# zn_links_mp.py -- written by the independent referee of Section 12 of Part IV (2026-11-12), added to the repository with credit;
# unchanged apart from this header. A compass, not part of the proofs.
# High-precision (mpmath) version of the covariance-matrix test of zn_links_tm.py for Z_3 on one cube (d = 4):
# entries of T are exact integer combinations of w0^(12-m) w1^m, evaluated at 50 digits; Omega by power iteration.
# Prints T1_qq'/(u^8 Ups^2) and G_qq'/(u^4 Ups^2) for adjacent and opposite faces, next to A_s + 2 C_s and C_s.
import sys
import numpy as np
import mpmath as mp
from geometry_H import Lattice, build_H
from zn_links_tm import build, T_matrix

mp.mp.dps = 50

def run(n=3, betas=('0.1', '0.05', '0.025', '0.0125')):
    lat = Lattice((1, 1, 1), False)
    data = build(lat, n)
    Hp, Hm, parts, _ = build_H(lat)
    Ncnt = data['Ncnt']; nl = data['nl']; norb = data['norb']; B = data['B']; P = lat.P()
    for bs in betas:
        beta = mp.mpf(bs)
        th = [2 * mp.pi * k / n for k in range(n)]
        w = [mp.e ** (beta * mp.cos(t)) for t in th]
        w0, w1 = w[0] / n, w[1] / n
        kern = [w0 ** (nl - m) * w1 ** m for m in range(nl + 1)]
        S = [sum(mp.cos(2 * mp.pi * int(B[p, o]) / n) for p in range(P)) for o in range(norb)]
        sv = [mp.e ** (beta * s / 2) for s in S]
        T = mp.matrix(norb, norb)
        for o in range(norb):
            for o2 in range(o, norb):
                val = mp.fsum(int(Ncnt[o, o2, m]) * kern[m] for m in range(nl + 1) if Ncnt[o, o2, m])
                T[o, o2] = T[o2, o] = sv[o] * sv[o2] * val
        # constants
        c0 = sum(w) / n
        u = sum(wk * mp.cos(t) for wk, t in zip(w, th)) / n / c0
        pr = [wk / sum(w) for wk in w]
        def var(f): return sum(p * f(t) ** 2 for p, t in zip(pr, th)) - sum(p * f(t) for p, t in zip(pr, th)) ** 2
        Ups = {'+': var(lambda t: mp.sqrt(2) * mp.cos(t)), '-': var(lambda t: mp.sqrt(2) * mp.sin(t))}
        # ground state by power iteration from the double-precision guess
        Td = T_matrix(data, float(beta)); ev, evec = np.linalg.eigh(Td)
        x = mp.matrix([mp.mpf(float(v)) for v in evec[:, -1]])
        for _ in range(12):
            y = T * x; nrm = mp.sqrt(mp.fsum(v ** 2 for v in y)); x = y / nrm
        lam0 = mp.fsum(x[i] * (T * x)[i] for i in range(norb))
        print(f'beta = {bs}: u = {mp.nstr(u, 8)}')
        for C in ('+', '-'):
            fn = (lambda b: mp.sqrt(2) * mp.cos(2 * mp.pi * b / n)) if C == '+' else (lambda b: mp.sqrt(2) * mp.sin(2 * mp.pi * b / n))
            vecs = []
            for q in range(P):
                v = mp.matrix([fn(int(B[q, o])) * x[o] for o in range(norb)])
                mean = mp.fsum(v[i] * x[i] for i in range(norb))
                vecs.append(v - mean * x)
            Tv = [T * v for v in vecs]
            G = [[mp.fsum(vecs[a][i] * vecs[b][i] for i in range(norb)) for b in range(P)] for a in range(P)]
            T1 = [[mp.fsum(vecs[a][i] * Tv[b][i] for i in range(norb)) / lam0 for b in range(P)] for a in range(P)]
            Y = Ups[C]
            As = parts['Ap'] if C == '+' else parts['Am']; Cs = parts['Cp'] if C == '+' else parts['Cm']
            rows = []
            for a in range(P):
                for b in range(P):
                    if a == b: continue
                    rows.append((a, b, float(T1[a][b] / (u ** 8 * Y ** 2)), As[a, b] + 2 * Cs[a, b], float(G[a][b] / (u ** 4 * Y ** 2)), Cs[a, b]))
            worst = max(abs(r[2] - r[3]) for r in rows); worstG = max(abs(r[4] - r[5]) for r in rows)
            diag = [float((T1[a][a] / (u ** 4 * Y ** 2) - 1) / u ** 4) for a in range(P)]
            print(f'  C{C}: Ups = {mp.nstr(Y, 10)}; max |T1/(u^8 Ups^2) - (A+2C)| = {worst:.4e}; max |G/(u^4 Ups^2) - C| = {worstG:.3e};'
                  f' diag (T1/(u^4Ups^2)-1)/u^4 in [{min(diag):.5f}, {max(diag):.5f}] (K_S = 4)')
            seen = set()
            for r in rows:
                key = (round(r[3]), round(r[5]))
                if key in seen: continue
                seen.add(key)
                print(f'     e.g. faces {r[0]},{r[1]}: T1/(u^8 Ups^2) = {r[2]:+.6f}  vs A+2C = {r[3]:+.0f};   G/(u^4 Ups^2) = {r[4]:+.6f} vs C = {r[5]:+.0f}')
        sys.stdout.flush()

if __name__ == '__main__':
    run()

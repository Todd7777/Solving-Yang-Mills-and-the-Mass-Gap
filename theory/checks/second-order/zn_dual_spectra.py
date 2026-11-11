# zn_dual_spectra.py -- written by the independent referee of Section 10 of Part IV (2026-11-11), added to the repository with
# credit; unchanged apart from this header and the run modes at the end. A compass, not part of the proofs.
# Runtime: z2, z3L3, z4L3 seconds; z2L5 about 12 min; z3 (L=4) about 30 min.
# Independent referee check of Corollary 10.9 by exact transfer-matrix spectra.
# Z_n lattice gauge theory in d = 3 (spatial L x L torus), Wilson weight w(m) = exp(beta cos(2 pi m/n)) (N = 1, w = exp(beta Re chi)).
# Zero-winding (electric flux) sector: states = functions of the plaquette variables b_p in Z_n with sum_p b_p = 0.
#   T = V^{1/2} K V^{1/2},  V = prod_p w(b_p),  K = prod_links k_l,  k_l = (1/n) sum_m w(m) tau_l^m,
#   where tau_l^m shifts the link variable by m, i.e. b_{p+} -> b_{p+} + m, b_{p-} -> b_{p-} - m for the two plaquettes containing l.
# The last plaquette variable is eliminated (b_last = -sum of the others).
# Lowest states in the sectors (C = +-, momentum 0 or (pi,pi)) by power iteration with deflation of the vacuum.
# Prediction of Section 10:  m = 4 log(1/u) - log Upsilon_C + O(beta^3), Upsilon_C the C-even/odd plaquette variance.
import numpy as np, sys, time

def setup(n, L):
    P = L*L
    idx = lambda x, y: (x % L)*L + (y % L)
    links = []
    for x in range(L):
        for y in range(L):
            links.append((idx(x, y), idx(x, y-1)))   # x-link at (x,y): + in p(x,y), - in p(x,y-1)
            links.append((idx(x-1, y), idx(x, y)))   # y-link at (x,y): + in p(x-1,y), - in p(x,y)
    return P, links, idx

class TM:
    def __init__(self, n, L, beta):
        self.n, self.L, self.beta = n, L, beta
        self.P, self.links, self.idx = setup(n, L)
        self.shape = (n,)*(self.P-1)
        m = np.arange(n)
        self.w = np.exp(beta*np.cos(2*np.pi*m/n))
        # plaquette values on the grid
        grids = np.indices(self.shape, dtype=np.int64) if self.P-1 <= 16 else None
        # build log V incrementally to save memory
        logV = np.zeros(self.shape)
        s = np.zeros(self.shape, dtype=np.int64)
        for ax in range(self.P-1):
            sh = [1]*(self.P-1); sh[ax] = n
            vals = np.arange(n).reshape(sh)
            logV = logV + np.log(self.w)[vals]
            s = s + vals
        last = (-s) % n
        logV = logV + np.log(self.w)[last]
        self.sqrtV = np.exp(0.5*(logV - logV.max()))
        self.last = last
        del s
        # character coefficients: u_r = c_r/c_0, c_r = (1/n) sum_m w(m) e^{-2 pi i r m/n}
        c = np.array([np.mean(self.w*np.exp(-2j*np.pi*r*m/n)) for r in range(n)])
        self.u = (c[1]/c[0]).real

    def b_of(self, p):
        """array of b_p over the grid (as a broadcastable view or full array for the last plaquette)"""
        if p == self.P-1:
            return self.last
        sh = [1]*(self.P-1); sh[p] = self.n
        return np.arange(self.n).reshape(sh)

    def K(self, psi):
        n = self.n
        for (pp, pm) in self.links:
            out = self.w[0]*psi
            for mm in range(1, n):
                shifts, axes = [], []
                if pp != self.P-1: shifts.append(mm); axes.append(pp)
                if pm != self.P-1: shifts.append(-mm); axes.append(pm)
                out += self.w[mm]*np.roll(psi, tuple(shifts), axis=tuple(axes))
            psi = out/n
        return psi

    def T(self, psi):
        return self.sqrtV*self.K(self.sqrtV*psi)

def upsilons(n, beta):
    m = np.arange(n); w = np.exp(beta*np.cos(2*np.pi*m/n)); p = w/w.sum()
    if n == 2:
        chi = np.cos(np.pi*m)
        return {'+': np.sum(p*chi**2) - np.sum(p*chi)**2}
    phis = {'+': np.sqrt(2)*np.cos(2*np.pi*m/n), '-': np.sqrt(2)*np.sin(2*np.pi*m/n)}
    # covariance matrix in the basis (sqrt2 cos, sqrt2 sin); for Z_n it is diagonal by C-symmetry
    return {k: np.sum(p*v**2) - np.sum(p*v)**2 for k, v in phis.items()}

def run(n, L, beta, sectors, iters=40, tol=1e-15, verbose=True):
    t0 = time.time()
    tm = TM(n, L, beta)
    u = tm.u
    # vacuum
    Om = tm.sqrtV.copy(); Om /= np.linalg.norm(Om)
    lam0_old = 0
    for it in range(60):
        v = tm.T(Om); lam0 = np.vdot(Om, v).real
        Om = v/np.linalg.norm(v)
        if abs(lam0-lam0_old) < 1e-15*abs(lam0): break
        lam0_old = lam0
    res = {}
    for (C, mom) in sectors:
        # start vector: sum_p phase(p) phi_C(b_p), times sqrt V (roughly the dressed plaquette state)
        psi = np.zeros(tm.shape)
        for x in range(L):
            for y in range(L):
                p = tm.idx(x, y)
                ph = 1.0 if mom == '0' else (-1.0)**(x+y)
                bp = tm.b_of(p)
                f = np.cos(2*np.pi*bp/n) if C == '+' else np.sin(2*np.pi*bp/n)
                psi = psi + ph*np.broadcast_to(f, tm.shape)
        psi = tm.sqrtV*psi
        psi -= np.vdot(Om, psi)*Om; psi /= np.linalg.norm(psi)
        lam_old = 0
        for it in range(iters):
            v = tm.T(psi)
            lam = np.vdot(psi, v).real
            v -= np.vdot(Om, v)*Om
            psi = v/np.linalg.norm(v)
            if abs(lam-lam_old) < tol*abs(lam): break
            lam_old = lam
        res[(C, mom)] = (lam, it)
    ups = upsilons(n, beta)
    out = []
    for (C, mom), (lam, it) in res.items():
        m = -np.log(lam/lam0)
        pred = 4*np.log(1/u) - np.log(ups[C])
        out.append((C, mom, m, pred, m - pred, it))
        if verbose:
            print(f"  Z{n} L={L} beta={beta:<8g} u={u:.6e}  C={C} k={mom}: m-4log(1/u)={m-4*np.log(1/u):+.10e}  "
                  f"-logUps_C={-np.log(ups[C]):+.10e}  diff={m-pred:+.3e}  diff/u^3={(m-pred)/u**3:+.4f} diff/u^4={(m-pred)/u**4:+.3f} (it={it})")
    if verbose: print(f"     [{time.time()-t0:.1f}s]"); sys.stdout.flush()
    return u, out

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else 'z2'
    if which == 'z2':
        for L in (4,):
            for beta in (0.2, 0.1, 0.05, 0.025):
                run(2, L, beta, [('+', '0'), ('+', 'pi')])
    elif which == 'z2L5':      # the 5x5 torus satisfies the hypotheses of Theorem mg:thm:lower2 (spatial sides >= 5)
        for beta in (0.1, 0.05):
            run(2, 5, beta, [('+', '0')], iters=25, tol=1e-12)
    elif which == 'z3':
        for beta in (0.1, 0.05, 0.025):
            run(3, 4, beta, [('+', '0'), ('-', '0'), ('+', 'pi'), ('-', 'pi')], iters=25, tol=1e-12)
    elif which == 'z3L3':
        for beta in (0.1, 0.05, 0.02, 0.01):
            run(3, 3, beta, [('+', '0'), ('-', '0')])
    elif which == 'z4L3':
        for beta in (0.1, 0.05, 0.02):
            run(4, 3, beta, [('+', '0'), ('-', '0')])

# Non-abelian check of the positivity decomposition: G = S_3 with its 2-dim irrep (real, iota=1),
# spatial graph = a single square (4 links), Lambda = (Z/L0) x square.
# (1) transfer matrix: T = V K1 P V, K1^(F) with k' = k - c0 Pi_0 on the links of F, c0 Pi_0 on E0\F;
#     check Tr(T^(L0-1) A~ T^(F) A~) >= 0 for all F and L0, and that they sum to the full trace.
# (2) brute force (L0 = 2, gauge-fixed on a spanning tree): the F-terms as expectations of
#     (A(0)-a)(A(1)-a) prod_{s in F}(1-c0/w_s) prod_{s in E0\F}(c0/w_s) agree with (1).
import numpy as np, itertools, math
# S_3 as 2x2 orthogonal matrices
c,s=math.cos(2*math.pi/3),math.sin(2*math.pi/3)
R=np.array([[c,-s],[s,c]]); Fl=np.array([[1,0],[0,-1]])
els=[np.eye(2),R,R@R,Fl,R@Fl,R@R@Fl]
nG=len(els)
def find(M):
    for k,E in enumerate(els):
        if np.allclose(M,E): return k
    raise ValueError
mul=np.array([[find(els[i]@els[j]) for j in range(nG)] for i in range(nG)])
inv=np.array([find(els[i].T) for i in range(nG)])
tr=np.array([np.trace(E) for E in els])
beta=0.7; N=2
w=np.exp(beta/N*tr)                 # w(U) = exp((beta/N) Re Tr rho(U))
c0=w.mean()
# --- transfer matrix on G^4 (links l_i from site i to i+1 mod 4)
states=list(itertools.product(range(nG),repeat=4)); ns=len(states)
sidx={st:i for i,st in enumerate(states)}
def plaq(st): return mul[mul[mul[st[0],st[1]],st[2]],st[3]]
Vd=np.array([math.sqrt(w[plaq(st)]) for st in states])
# one-link kernels as matrices w.r.t. normalised counting measure
k=np.array([[w[mul[a,inv[b]]] for b in range(nG)] for a in range(nG)])/nG
Pi0=np.ones((nG,nG))/nG
kp=k-c0*Pi0
assert np.linalg.eigvalsh(kp).min()>-1e-12, np.linalg.eigvalsh(kp)
def kron4(ms):
    out=ms[0]
    for m in ms[1:]: out=np.kron(out,m)
    return out
# gauge projection: average over g in G^4 of R_g, (R_g psi)(U) = psi(g^{-1}.U), (g.U)_i = g_i U_i g_{i+1}^{-1}
Pm=np.zeros((ns,ns))
for g in itertools.product(range(nG),repeat=4):
    gi=[inv[x] for x in g]
    for i,st in enumerate(states):
        # (g^{-1}.U)_l = g_l^{-1} U_l g_{l+1}
        t=tuple(mul[mul[gi[l],st[l]],g[(l+1)%4]] for l in range(4))
        Pm[i,sidx[t]]+=1
Pm/=nG**4
A=np.array([tr[plaq(st)] for st in states])      # A = Tr rho(U_p) (real)
def Tmat(Fset):
    if Fset is None: K1=kron4([k]*4)
    else: K1=kron4([kp if i in Fset else c0*Pi0 for i in range(4)])
    return (Vd[:,None]*(K1@Pm))*Vd[None,:]
T=Tmat(None)
TF={F:Tmat(set(F)) for r in range(5) for F in itertools.combinations(range(4),r)}
assert np.allclose(sum(TF.values()),T)
for F,M in TF.items(): assert np.abs(M-M.T).max()<1e-14 and np.linalg.eigvalsh((M+M.T)/2).min()>-1e-12*np.abs(T).max()
ev,EV=np.linalg.eigh(T)
print("T eigenvalues (top 4):",ev[-4:])
for L0 in [2,3,4,6,10]:
    TL=np.linalg.matrix_power(T,L0); Z=np.trace(TL)
    T1=np.linalg.matrix_power(T,L0-1)
    a=np.trace(TL*A[None,:])/Z if False else np.trace(T1@T@np.diag(A))/Z  # <A(0)>
    At=np.diag(A-a)
    vals={F:np.trace(T1@At@M@At)/Z for F,M in TF.items()}
    tot=np.trace(T1@At@T@At)/Z
    assert abs(sum(vals.values())-tot)<1e-10*max(1,abs(tot))
    print(f"L0={L0}: Cov={tot:.6e}, F=E0 term={vals[(0,1,2,3)]:.6e}, min over F={min(vals.values()):+.3e}")
# --- brute force, L0 = 2, gauge fixed: V0_i = 1 (temporal links of layer 0) and U0_0=U0_1=U0_2=1
TL=T@T; Z=np.trace(TL); a=np.trace(TL@np.diag(A))/Z
idxs=np.array(list(itertools.product(range(nG),repeat=9)),dtype=np.int64)  # U0_3, U1_0..3, V1_0..3
e=0
U0=[np.zeros(len(idxs),dtype=np.int64)]*3+[idxs[:,0]]
U1=[idxs[:,1+i] for i in range(4)]
V1=[idxs[:,5+i] for i in range(4)]    # temporal link from (i,1) to (i,0) (layer 1, time 1 -> 2=0)
def prod(*xs):
    out=xs[0]
    for x in xs[1:]: out=mul[out,x]
    return out
p_s0=prod(U0[0],U0[1],U0[2],U0[3]); p_s1=prod(U1[0],U1[1],U1[2],U1[3])
# layer 0 sides (V0 = 1): holonomy U0_i * 1 * U1_i^{-1} * 1
side0=[prod(U0[i],inv[U1[i]]) for i in range(4)]
# layer 1 sides: U1_i V1_{i+1} U0_i^{-1} V1_i^{-1}
side1=[prod(U1[i],V1[(i+1)%4],inv[U0[i]],inv[V1[i]]) for i in range(4)]
Wrest=w[p_s0]*w[p_s1]*np.prod([w[x] for x in side1],axis=0)
Wside=np.prod([w[x] for x in side0],axis=0)
Zb=np.sum(Wrest*Wside)
ab=np.sum(tr[p_s0]*Wrest*Wside)/Zb
print(f"brute force: a={ab:.10f} vs transfer {a:.10f}")
base=(tr[p_s0]-ab)*(tr[p_s1]-ab)*Wrest
for F in [(),(0,),(0,1),(0,2),(0,1,2),(0,1,2,3)]:
    fac=np.ones(len(idxs))
    for i in range(4):
        ws=w[side0[i]]
        fac*= (ws-c0) if i in F else c0
    bf=np.sum(base*fac)/Zb
    T1=T; At=np.diag(A-a)
    tm=np.trace(T1@At@TF[F]@At)/Z
    print(f"F={F}: brute force {bf:+.10e}  transfer {tm:+.10e}")

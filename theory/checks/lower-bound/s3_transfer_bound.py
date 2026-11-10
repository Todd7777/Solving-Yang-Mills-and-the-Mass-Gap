# Check of Lemma mg:lem:transfer (operator-norm bound for duplicated terms) for S_3 (2-dim irrep) on
# Lambda = (Z/2) x (one square).  (1) Validate the single-copy transfer representation against brute force
# (gauge-fixed, 6^9 configurations) on random (H, Q).  (2) Assemble Phi_F(X) from single-copy terms via the
# transfer representation and compare with the bound.
import numpy as np, itertools, math, random
c,s=math.cos(2*math.pi/3),math.sin(2*math.pi/3)
R=np.array([[c,-s],[s,c]]); Fl=np.array([[1,0],[0,-1]])
els=[np.eye(2),R,R@R,Fl,R@Fl,R@R@Fl]; nG=6
def find(M):
    for k,E in enumerate(els):
        if np.allclose(M,E): return k
mul=np.array([[find(els[i]@els[j]) for j in range(nG)] for i in range(nG)])
inv=np.array([find(els[i].T) for i in range(nG)])
tr=np.array([np.trace(E) for E in els])
sign=np.array([1,1,1,-1,-1,-1.]); chis={'sign':sign,'two':tr}; dims={'sign':1,'two':2}
beta=0.9; N=2
w=np.exp(beta/N*tr); c0=w.mean(); f=w/c0-1
u={k:np.mean(chis[k]*w)/(dims[k]*c0) for k in chis}; ustar=max(u.values())
gam=np.exp(4*beta)-1; TrF=np.exp(beta)/c0-1; TrG=(1+TrF)**2-1
print(f"u_sign={u['sign']:.5f} u_two={u['two']:.5f} u*={ustar:.5f}  TrG={TrG:.5f} <= e^(2b)-1={np.exp(2*beta)-1:.5f}")
# transfer side
states=list(itertools.product(range(nG),repeat=4)); ns=len(states); sidx={st:i for i,st in enumerate(states)}
def plq(st): return mul[mul[mul[st[0],st[1]],st[2]],st[3]]
Pm=np.zeros((ns,ns))
for g in itertools.product(range(nG),repeat=4):
    gi=[inv[x] for x in g]
    for i,st in enumerate(states):
        t=tuple(mul[mul[gi[l],st[l]],g[(l+1)%4]] for l in range(4)); Pm[i,sidx[t]]+=1
Pm/=nG**4
Pi0=np.ones((nG,nG))/nG
Ff=np.array([[f[mul[a,inv[b]]] for b in range(nG)] for a in range(nG)])/nG
ev=np.linalg.eigvalsh(Ff); assert abs(ev.max()-ustar)<1e-12 and ev.min()>-1e-12
def kron4(ms):
    out=ms[0]
    for m in ms[1:]: out=np.kron(out,m)
    return out
LAY={}
for r in range(5):
    for T in itertools.combinations(range(4),r):
        LAY[frozenset(T)]=kron4([Ff if i in T else Pi0 for i in range(4)])@Pm
fS=np.array([f[plq(st)] for st in states]); Ad=np.array([tr[plq(st)] for st in states]).astype(float)
_tc={}
def I_tm(H0,H1,Q):
    key=(H0,H1,frozenset(Q))
    if key in _tc: return _tc[key]
    T0=frozenset(q[2] for q in Q if q[0]=='T' and q[1]==0); T1=frozenset(q[2] for q in Q if q[0]=='T' and q[1]==1)
    m0=np.ones(ns); m1=np.ones(ns)
    if ('S',0) in Q: m0=m0*fS
    if ('S',1) in Q: m1=m1*fS
    if H0: m0=m0*Ad
    if H1: m1=m1*Ad
    val=np.sum(m0[:,None]*LAY[T0]*m1[None,:]*LAY[T1].T)
    _tc[key]=val; return val
# brute force side (single copy), gauge fixed
idx=np.array(list(itertools.product(range(nG),repeat=9)),dtype=np.int64)
U0=[np.zeros(len(idx),dtype=np.int64)]*3+[idx[:,0]]; U1=[idx[:,1+i] for i in range(4)]; V1=[idx[:,5+i] for i in range(4)]
def prod(*xs):
    out=xs[0]
    for x in xs[1:]: out=mul[out,x]
    return out
plaq={('S',0):prod(U0[0],U0[1],U0[2],U0[3]),('S',1):prod(U1[0],U1[1],U1[2],U1[3])}
for i in range(4):
    plaq[('T',0,i)]=prod(U0[i],inv[U1[i]]); plaq[('T',1,i)]=prod(U1[i],V1[(i+1)%4],inv[U0[i]],inv[V1[i]])
A0=tr[plaq[('S',0)]].astype(float); B1=tr[plaq[('S',1)]].astype(float)
def I_bf(H0,H1,Q):
    val=np.ones(len(idx))
    if H0: val=val*A0
    if H1: val=val*B1
    for q in Q: val=val*f[plaq[q]]
    return val.mean()
allp=[('S',0),('S',1)]+[('T',0,i) for i in range(4)]+[('T',1,i) for i in range(4)]
random.seed(1); maxerr=0
for trial in range(40):
    Q=[p for p in allp if random.random()<0.45]; H0=random.random()<0.5; H1=random.random()<0.5
    a=I_bf(H0,H1,Q); b=I_tm(H0,H1,Q); maxerr=max(maxerr,abs(a-b))
print(f"transfer representation vs brute force, 40 random (H,Q): max |difference| = {maxerr:.2e}")
tests=[[('T',0,i) for i in range(4)],
       [('T',0,i) for i in range(4)]+[('S',0)],
       [('T',0,i) for i in range(4)]+[('T',1,i) for i in range(4)],
       [('T',0,i) for i in range(4)]+[('T',1,0),('T',1,1)]+[('S',1)],
       [('T',0,i) for i in range(4)]+[('S',0),('S',1),('T',1,2)],
       [('T',0,0),('T',0,1),('T',0,2),('T',1,0),('T',1,3),('S',0),('S',1)]]
for X in tests:
    phi=0.0
    for assign in itertools.product(range(3),repeat=len(X)):
        Q1=[X[j] for j in range(len(X)) if assign[j] in (0,2)]; Q2=[X[j] for j in range(len(X)) if assign[j] in (1,2)]
        phi+=(I_tm(1,1,Q1)*I_tm(0,0,Q2)-I_tm(1,0,Q1)*I_tm(0,1,Q2)-I_tm(0,1,Q1)*I_tm(1,0,Q2)+I_tm(0,0,Q1)*I_tm(1,1,Q2))
    TX=[q for q in X if q[0]=='T']; SX=[q for q in X if q[0]=='S']
    mins=min(sum(1 for q in TX if q[1]==l) for l in (0,1))
    bound=4*2*2*gam**len(SX)*(np.exp(2*beta)-1)**mins*ustar**(len(TX)-mins)   # 4||A|| ||B||, ||A|| = N = 2
    print(f"|X|={len(X)} |T|={len(TX)} |S|={len(SX)} m(X)={mins}: Phi_F={phi:+.6e}  bound={bound:.3e}  |Phi|/bound={abs(phi)/bound:.2e}")

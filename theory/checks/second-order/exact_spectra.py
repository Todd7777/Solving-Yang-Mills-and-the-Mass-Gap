# Section 10 of Part IV, Conjecture mg:conj:second and Remark mg:rem:second: exact transfer-matrix spectra against
#   m = 4 log(1/u) - log Lambda_C + O(beta^3),  Lambda_C = Var_{w/c0}(phi_C)  (phi_C the unit plaquette function of the C-sector).
# (1) Z_n, d=3, LxL spatial torus, winding-zero sector, both C-sectors (Z_2: 4x4; Z_3, Z_4, Z_5: 3x3).
# (2) Spatial lattice = one square (class functions): SU(2) (characters), S_3 and Q_8 (two-dimensional representations).
# (3) Spatial lattice = two squares sharing an edge (full link space G^7 with the gauge projection): S_3, Q_8 -- this
#     includes the shift of the tube between neighbouring plaquettes, which cancels at second order.
# A compass, not part of the proofs. Usage: python3 exact_spectra.py [zn|square|ladder|all]
import numpy as np, math, itertools, sys
from scipy.sparse.linalg import LinearOperator, eigsh
def run(nG, L, beta, kev=2):
    nl=2*L*L; NP=L*L; ng=NP-1; shape=(nG,)*ng
    def lid(x,y,i): return ((x%L)*L+(y%L))*2+i
    def bd(x,y):
        v=np.zeros(nl,dtype=np.int64); v[lid(x,y,0)]+=1; v[lid(x+1,y,1)]+=1; v[lid(x,y+1,0)]-=1; v[lid(x,y,1)]-=1; return v
    P=[bd(x,y) for x in range(L) for y in range(L)]
    labels=(np.indices(shape).reshape(ng,-1).T@np.array(P[:ng]))%nG
    k=np.arange(nG); w=np.exp(beta*np.cos(2*np.pi*k/nG)); c0=w.mean()
    cr=np.array([np.real(np.mean(w*np.exp(-2j*np.pi*r*k/nG))) for r in range(nG)])
    wh=np.exp(beta*np.cos(2*np.pi*k/nG)/2); br=np.array([np.real(np.mean(wh*np.exp(-2j*np.pi*r*k/nG))) for r in range(nG)])
    Kd=np.prod((cr/cr[0])[labels],axis=1)
    def shift(v,p,s):
        if p<ng: return np.roll(v,s,axis=p)
        return np.roll(v,(-s,)*ng,axis=tuple(range(ng)))
    def Vop(v):
        v=v.reshape(shape)
        for p in range(NP):
            acc=br[0]*v
            for r in range(1,nG): acc=acc+br[r]*shift(v,p,r)
            v=acc/br[0]
        return v.reshape(-1)
    idxC=np.ix_(*[np.mod(-np.arange(nG),nG)]*ng)
    def Cop(v): return np.ascontiguousarray(v.reshape(shape)[idxC]).reshape(-1)
    n=nG**ng; out={}
    for sgn in ((+1,-1) if nG>2 else (+1,)):
        Pc=lambda v: 0.5*(v+sgn*Cop(v)) if nG>2 else v
        op=LinearOperator((n,n),matvec=lambda v: Pc(Vop(Kd*Vop(Pc(v)))),dtype=float)
        out[sgn]=np.sort(eigsh(op,k=kev,which='LA',tol=1e-14,return_eigenvectors=False))[::-1]
    u=cr[1]/cr[0]; lam0=out[1][0]
    res={}
    for sgn in out:
        lam=out[sgn][1] if sgn==1 else out[sgn][0]
        res[sgn]=math.log(lam0/lam)-4*math.log(1/u)
    # conjecture: Var under the single-plaquette measure p = w/(n c0)
    p=w/w.sum(); th=2*np.pi*k/nG
    def var(phi): return np.sum(p*phi**2)-np.sum(p*phi)**2
    if nG>2: conj={+1:-math.log(var(math.sqrt(2)*np.cos(th))), -1:-math.log(var(math.sqrt(2)*np.sin(th)))}
    else: conj={+1:-math.log(var(np.cos(th)))}
    return u,res,conj
def class_data(G):
    els, mul, inv, irreps = G
    n=len(els); idx={g:i for i,g in enumerate(els)}
    # conjugacy classes
    seen=set(); classes=[]
    for g in els:
        if g in seen: continue
        cl={mul(mul(h,g),inv(h)) for h in els}; seen|=cl; classes.append(sorted(cl,key=lambda x: idx[x]))
    w=np.array([len(c)/n for c in classes])
    X={r:np.array([np.trace(R(c[0])) for c in classes]) for r,R in irreps.items()}
    d={r:irreps[r](els[0]).shape[0] for r in irreps}
    return w,X,d
def S3():
    els=list(itertools.permutations(range(3))); mul=lambda a,b: tuple(a[b[i]] for i in range(3)); inv=lambda a: tuple(np.argsort(a))
    B=np.array([[1,-1,0],[1,1,-2]],float); B[0]/=np.sqrt(2); B[1]/=np.sqrt(6)
    def std(p):
        Pm=np.zeros((3,3))
        for i in range(3): Pm[p[i],i]=1
        return B@Pm@B.T
    sgn=lambda p: round(np.linalg.det(np.eye(3)[list(p)]))
    return els,mul,inv,{'0':lambda g:np.eye(1),'sgn':lambda g:np.array([[float(sgn(g))]]),'2':std}
def Q8():
    # quaternion units as 2x2 complex matrices
    I2=np.eye(2); i=np.array([[1j,0],[0,-1j]]); j=np.array([[0,1],[-1,0]]); k=i@j
    mats=[I2,-I2,i,-i,j,-j,k,-k]; keys=[tuple(np.round(m,6).ravel()) for m in mats]
    els=list(range(8)); mul=lambda a,b: keys.index(tuple(np.round(mats[a]@mats[b],6).ravel())); inv=lambda a: keys.index(tuple(np.round(np.linalg.inv(mats[a]),6).ravel()))
    G=(els,mul,inv,{})
    # irreps: trivial, three 1-dim sign reps, 2-dim
    def one(sig): return lambda g: np.array([[ {0:1,1:1,2:sig[0],3:sig[0],4:sig[1],5:sig[1],6:sig[0]*sig[1],7:sig[0]*sig[1]}[g] ]])
    irreps={'0':one((1,1)),'a':one((1,-1)),'b':one((-1,1)),'c':one((-1,-1)),'2':lambda g: mats[g]}
    return els,mul,inv,irreps
def gap_classfun(w,X,d,rho,beta):
    N=d[rho]; reps=list(X); chi=X[rho]
    wt=np.exp(beta/N*np.real(chi)); c0=np.sum(w*wt)
    c={r:np.real(np.sum(w*np.conj(X[r])*wt)) for r in reps}; u={r:c[r]/(d[r]*c0) for r in reps}
    # basis: orthonormal class functions chi_r (L2 w.r.t. class weights w)
    M=np.array([[np.sum(w*np.conj(X[a])*(wt/c0)*X[b]) for b in reps] for a in reps])
    Kh=np.diag([u[r]**2 if r!='0' else 1.0 for r in reps])
    T=Kh@M@Kh; ev=np.sort(np.linalg.eigvalsh((T+T.conj().T)/2).real)[::-1]
    m=math.log(ev[0]/ev[1]); uu=u[rho]
    f=wt/c0-1; p=w*(1+f)
    iota=1 if np.allclose(chi,np.conj(chi)) else 0
    if iota: phis=[chi.real]
    else: phis=[np.sqrt(2)*chi.real, np.sqrt(2)*chi.imag]
    Lam=max(np.sum(p*ph**2)-np.sum(p*ph)**2 for ph in phis)
    return uu, m-4*math.log(1/uu), -math.log(Lam)
def su2_classfun(beta,J=12):
    # Weyl integration on theta in [0,pi]: dmu = (2/pi) sin^2
    M_=4000; th=(np.arange(M_)+0.5)*np.pi/M_; w=(2/np.pi)*np.sin(th)**2*np.pi/M_
    X={str(jj/2):np.sin((jj+1)*th)/np.sin(th) for jj in range(0,J)}; d={str(jj/2):jj+1 for jj in range(J)}
    X={('0' if k=='0.0' else k):v for k,v in X.items()}; d={('0' if k=='0.0' else k):v for k,v in d.items()}
    return gap_classfun(w,X,d,'0.5',beta)
V_=[(x,y) for y in range(2) for x in range(3)]; vid={v:i for i,v in enumerate(V_)}
links=[((0,0),(1,0)),((1,0),(2,0)),((0,1),(1,1)),((1,1),(2,1)),((0,0),(0,1)),((1,0),(1,1)),((2,0),(2,1))]
# plaquettes as oriented link sequences: left square (0,0)->(1,0)->(1,1)->(0,1)->(0,0): +L0 +L5 -L2 -L4
plaq=[[(0,+1),(5,+1),(2,-1),(4,-1)],[(1,+1),(6,+1),(3,-1),(5,-1)]]
def build(G,rho,beta):
    els,mul,inv,irreps=G; n=len(els); idx={g:i for i,g in enumerate(els)}
    M=np.array([[idx[mul(a,b)] for b in els] for a in els]); I=np.array([idx[inv(a)] for a in els])
    chi=np.array([np.trace(irreps[rho](g)) for g in els]); N=irreps[rho](els[0]).shape[0]
    wt=np.exp(beta/N*np.real(chi)); c0=wt.mean(); onef=wt/c0
    # one-link kernel (1+f)(U V^{-1})/|G|
    Kl=np.array([[onef[M[a,I[b]]] for b in range(n)] for a in range(n)])/n
    L=len(links); shape=(n,)*L; grid=np.indices(shape).reshape(L,-1)
    def hol(p):
        h=np.zeros(grid.shape[1],dtype=np.int64)  # identity index
        h[:]=idx[els[0]] if False else 0
        e0=idx[mul(els[0],inv(els[0]))]; h[:]=e0
        for (l,s) in p:
            g=grid[l] if s>0 else I[grid[l]]
            h=M[h,g]
        return h
    Vh=np.ones(grid.shape[1])
    for p in plaq: Vh*=np.sqrt(onef[hol(p)])
    # gauge transformations at vertex v by group element g: U_l -> g U_l if l starts at v ; U_l -> U_l g^{-1} if l ends at v
    vert_perms=[]
    for v in V_:
        perms=[]
        for g in range(n):
            new=grid.copy()
            for li,(a,b) in enumerate(links):
                if a==v: new[li]=M[g,new[li]]
                if b==v: new[li]=M[new[li],I[g]]
            perms.append(np.ravel_multi_index(new,shape))
        vert_perms.append(perms)
    def P(x):
        for perms in vert_perms:
            acc=np.zeros_like(x)
            for pm in perms: acc+=x[pm]
            x=acc/n
        return x
    def Kop(x):
        x=x.reshape(shape)
        for l in range(L):
            x=np.moveaxis(np.tensordot(Kl,x,axes=([1],[l])),0,l)
        return x.reshape(-1)
    def T(x): return P(Vh*Kop(Vh*P(x)))
    Ntot=n**L
    op=LinearOperator((Ntot,Ntot),matvec=T,dtype=float)
    ev=np.sort(eigsh(op,k=4,which='LA',tol=1e-13,return_eigenvectors=False))[::-1]
    w,X,d=class_data(G); reps=list(X)
    c={r:np.real(np.sum(w*np.conj(X[r])*np.exp(beta/N*np.real(X[rho])))) for r in reps}
    u=c[rho]/(N*c['0'] if '0' in c else 1)
    f=np.exp(beta/N*np.real(X[rho]))/np.sum(w*np.exp(beta/N*np.real(X[rho])))-1; pw=w*(1+f); ch=X[rho]
    iota=1 if np.allclose(ch,np.conj(ch)) else 0
    phis=[ch.real] if iota else [np.sqrt(2)*ch.real,np.sqrt(2)*ch.imag]
    Lam=max(np.sum(pw*ph**2)-np.sum(pw*ph)**2 for ph in phis)
    return u,ev,-math.log(Lam)
if __name__=="__main__":
    what=sys.argv[1] if len(sys.argv)>1 else 'all'
    if what in ('zn','all'):
        print("(1) Z_n tori, winding zero: m - 4log(1/u) exact vs -log Lambda_C")
        for nG,L in ((2,4),(3,3),(4,3),(5,3)):
            for beta in (0.02,0.005):
                u,res,conj=run(nG,L,beta)
                s='  '.join(f"C{'+' if c>0 else '-'}: exact {res[c]:+.6e} conj {conj[c]:+.6e} diff/u^3 {(res[c]-conj[c])/u**3:+.2f}" for c in res)
                print(f"  Z_{nG} {L}x{L} beta={beta}: u={u:.4e}  {s}")
    if what in ('square','all'):
        print("(2) one square")
        for beta in (0.04,0.01):
            u,ex,cj=su2_classfun(beta); print(f"  SU(2) beta={beta}: exact {ex:+.6e} conj {cj:+.6e} diff/beta^3 {(ex-cj)/beta**3:+.4f}")
        for name,G in (('S3',S3()),('Q8',Q8())):
            w,X,d=class_data(G)
            for beta in (0.04,0.01):
                u,ex,cj=gap_classfun(w,X,d,'2',beta); print(f"  {name} beta={beta}: exact {ex:+.6e} conj {cj:+.6e} diff/beta^3 {(ex-cj)/beta**3:+.4f}")
    if what in ('ladder','all'):
        print("(3) two squares sharing an edge (the two lowest plaquette states, symmetric and antisymmetric)")
        for name,G in (('S3',S3()),('Q8',Q8())):
            for beta in (0.04,0.01):
                u,ev,cj=build(G,'2',beta); ms=[math.log(ev[0]/e)-4*math.log(1/u) for e in ev[1:3]]
                print(f"  {name} beta={beta}: exact {ms[0]:+.6e}, {ms[1]:+.6e}; conj {cj:+.6e}; diff/beta^3 {(ms[0]-cj)/beta**3:+.4f}")

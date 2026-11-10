# Lemma mg:lem:capchar (Section 9 of Part IV), by brute force: on L^2(G^{E(q)} x G^{E(q)}) (the four links of a square, two
# copies) O_cap = (G_g^{1/2})^{(x)4} M_{g_q} (G_g^{1/2})^{(x)4} and P^q = the gauge projection at the four corners (each copy).
# ||P^q O_cap P^q|| (matrix-free Lanczos on the full space) against the class-function formula ||Lambda M_g Lambda|| of the proof
# and the bound u^4 kbar_cap of the lemma, for Z_3 and for S_3 with its two-dimensional representation (where kappa ~ beta/2 != 0).
import numpy as np, itertools, math
from scipy.sparse.linalg import LinearOperator, eigsh

def group_Z3():
    els=list(range(3)); mul=lambda a,b:(a+b)%3; inv=lambda a:(-a)%3
    om=np.exp(2j*np.pi/3)
    irreps={'0':lambda g:np.array([[1]]),'1':lambda g:np.array([[om**g]]),'2':lambda g:np.array([[om**(2*g)]])}
    return els,mul,inv,irreps,'1'
def group_S3():
    els=list(itertools.permutations(range(3)))
    mul=lambda a,b: tuple(a[b[i]] for i in range(3))
    inv=lambda a: tuple(np.argsort(a))
    def sgn(p): return round(np.linalg.det(np.eye(3)[list(p)]))
    # standard 2-dim rep: restriction of permutation rep to sum-zero plane, orthonormal basis
    B=np.array([[1,-1,0],[1,1,-2]],float); B[0]/=np.sqrt(2); B[1]/=np.sqrt(6)
    def std(p):
        Pm=np.zeros((3,3)); 
        for i in range(3): Pm[p[i],i]=1
        return B@Pm@B.T
    irreps={'0':lambda g:np.array([[1.]]),'sgn':lambda g:np.array([[float(sgn(g))]]),'2':std}
    return els,mul,inv,irreps,'2'

def run(G, beta, twocopy=True):
    els,mul,inv,irreps,rho=G; n=len(els); idx={g:i for i,g in enumerate(els)}
    chi={r:np.array([np.trace(R(g)) for g in els]) for r,R in irreps.items()}
    d={r:irreps[r](els[0]).shape[0] for r in irreps}
    N=d[rho]; iota=1 if np.allclose(chi[rho],np.conj(chi[rho])) else 0; tau0=(1+iota)/2
    w=np.exp(beta/N*np.real(chi[rho])); c0=w.mean(); f=w/c0-1
    c={r:np.real(np.mean(np.conj(chi[r])*w)) for r in irreps}; u={r:c[r]/(d[r]*c0) for r in irreps}
    # one-link operator F (kernel f(U V^{-1})), Pi0
    F=np.array([[f[idx[mul(els[i],inv(els[j]))]] for j in range(n)] for i in range(n)])/n
    Pi0=np.ones((n,n))/n
    Gg=np.kron(Pi0+F,Pi0+F)-np.kron(Pi0,Pi0)
    ev,evec=np.linalg.eigh(Gg); ev=np.clip(ev,0,None); Gh=(evec*np.sqrt(ev))@evec.T
    # holonomy of the square: U_q = U1 U2 U3^{-1} U4^{-1}
    def hol(t): return mul(mul(els[t[0]],els[t[1]]),inv(mul(els[t[3]],els[t[2]])))
    confs=list(itertools.product(range(n),repeat=4))
    fq=np.array([f[idx[hol(t)]] for t in confs])
    # gauge transformations at corners a,b,c,d: U1: a->b, U2: b->c, U3: d->c, U4: a->d  (U_l -> g_start U_l g_end^{-1})
    def gauge_perm(ga,gb,gc,gd):
        perm=np.empty(len(confs),dtype=np.int64)
        for k,t in enumerate(confs):
            U=[els[i] for i in t]
            V=[mul(mul(ga,U[0]),inv(gb)),mul(mul(gb,U[1]),inv(gc)),mul(mul(gd,U[2]),inv(gc)),mul(mul(ga,U[3]),inv(gd))]
            perm[k]=np.ravel_multi_index([idx[x] for x in V],(n,)*4)
        return perm
    # corner projections: average over one corner's gauge group (generated one corner at a time)
    e0=els[0]
    corner_perms=[[gauge_perm(*(g if j==k else e0 for j in range(4))) for g in els] for k in range(4)]
    def Pcorner1(v,k,axis_copy):  # v shape (n^4, n^4): apply corner k projection on copy axis
        out=np.zeros_like(v)
        for p in corner_perms[k]:
            out+= v[p,:] if axis_copy==0 else v[:,p]
        return out/n
    hol_ok=True
    # class-function side: Lambda M_g Lambda
    reps=list(irreps); 
    lam1={r:(u[r] if r!='0' else 0.0) for r in reps}
    def lam(a,b):
        if a=='0' and b=='0': return 0.0
        if b=='0': return u[a]
        if a=='0': return u[b]
        return u[a]*u[b]
    m=np.array([[np.mean(np.conj(chi[a])*(1+f)*chi[b]) for b in reps] for a in reps])
    pairs=[(a,b) for a in reps for b in reps]
    Mg=np.array([[m[reps.index(a),reps.index(b)]*m[reps.index(a2),reps.index(b2)]-(a==b and a2==b2) for (b,b2) in pairs] for (a,a2) in pairs])
    Lam=np.diag([lam(a,a2)**2 for (a,a2) in pairs])
    norm_class=np.max(np.abs(np.linalg.eigvalsh((Lam@Mg@Lam+ (Lam@Mg@Lam).conj().T)/2)))
    # kappa and the bound
    phi=lambda a,b: np.real(np.mean(np.conj(chi[a])*f*chi[b]))
    kappa=np.real(np.mean(np.real(chi[rho])**2*f))/tau0
    eps=math.expm1(2*beta); others=[r for r in reps if r not in ('0',rho) and not np.allclose(chi[r],np.conj(chi[rho]))]
    u2=max([u[r] for r in others],default=0.0); ur=u[rho]
    gam=np.max(np.abs(np.add.outer(f,f)+np.multiply.outer(f,f)))
    s6=sum(d[r]**2*u[r]**6 for r in reps if r!='0')
    kbar=kappa+(u2/ur)**2*eps*(1+(u2/ur)**2)+s6/ur**4+math.sqrt(2)*math.exp(beta)*math.sqrt(s6)+ur**4*gam
    out=dict(u=ur,kappa=kappa,norm_class=norm_class/ur**4,kbar=kbar,lead=kappa+(2-iota)*N**2*ur**2)
    # brute force (matrix-free), two copies
    D1=n**4
    def apply_link_op(v,Op):   # v shape (n,)*8 ordered (U1..U4,U1'..U4'); Op acts on (U_l,U_l') pairs for each link
        v=v.reshape((n,)*8)
        for l in range(4):
            v=np.moveaxis(v,(l,4+l),(0,1)); sh=v.shape
            v=(Op@v.reshape(n*n,-1)).reshape(sh)
            v=np.moveaxis(v,(0,1),(l,4+l))
        return v.reshape(D1,D1)
    def Pq(v):
        for k in range(4):
            v=Pcorner1(v,k,0); v=Pcorner1(v,k,1)
        return v
    g2=np.add.outer(fq,fq)+np.multiply.outer(fq,fq)
    def Oapply(x):
        v=Pq(x.reshape(D1,D1)); v=apply_link_op(v,Gh); v=g2*v; v=apply_link_op(v,Gh); v=Pq(v)
        return v.reshape(-1)
    if twocopy:
        op=LinearOperator((D1*D1,D1*D1),matvec=Oapply,dtype=float)
        w_=eigsh(op,k=2,which='LM',tol=1e-10,return_eigenvectors=False)
        out['norm_brute']=np.max(np.abs(w_))/ur**4
    return out
for name,G in (('Z3',group_Z3()),('S3',group_S3())):
    for beta in (0.05,0.01):
        r=run(G,beta,twocopy=True)
        print(f"{name} beta={beta}: u={r['u']:.4e} kappa={r['kappa']:.5e}  ||PO_capP||/u^4: brute={r.get('norm_brute',float('nan')):.6e} class={r['norm_class']:.6e}  bound kbar={r['kbar']:.6e}  kappa+(2-i)N^2u^2={r['lead']:.6e}")

# Exact check for Z_3 gauge theory (complex representation, iota=0, tau0=1/2), d=3, spatial torus 3x3,
# infinite time extent. Physical space: characters chi_n, n in Z_3^{links} divergence-free (dim 3^10).
import numpy as np, math, itertools
from scipy.sparse.linalg import LinearOperator, eigsh
from constants_upper import consts
L=3
def lid(x,y,i): return ((x%L)*L+(y%L))*2+i
nl=2*L*L
# oriented boundary of plaquette at (x,y): +l(x,y,0) +l(x+1,y,1) -l(x,y+1,0) -l(x,y,1)
def bd(x,y):
    v=np.zeros(nl,dtype=np.int64)
    v[lid(x,y,0)]+=1; v[lid(x+1,y,1)]+=1; v[lid(x,y+1,0)]-=1; v[lid(x,y,1)]-=1
    return v
P=[bd(x,y) for x in range(L) for y in range(L)]
assert np.all(sum(P)%3==0)
W1=np.zeros(nl,dtype=np.int64); W2=np.zeros(nl,dtype=np.int64)
for x in range(L): W1[lid(x,0,0)]=1
for y in range(L): W2[lid(0,y,1)]=1
gens=P[:8]+[W1,W2]
shape=(3,)*10
# label configuration for every basis state
grid=np.indices(shape).reshape(10,-1).T          # (59049,10)
labels=(grid@np.array(gens))%3                   # (59049,18)
assert len({tuple(r) for r in labels[:5000]})==5000
p0links=[lid(0,0,0),lid(1,0,1),lid(0,1,0),lid(0,0,1)]  # links of plaquette P[0]
om=np.exp(2j*np.pi/3)
def coeffs(fun):
    return np.array([np.real(sum(np.conj(om**(k*m))*fun(om**m) for m in range(3))/3) for k in range(3)])
def run(beta):
    c=coeffs(lambda U: math.exp(beta*U.real))          # c_k
    b=coeffs(lambda U: math.exp(beta*U.real/2))        # coefficients of w^{1/2}
    Kd=np.prod(c[labels],axis=1)                        # K_1 P on physical characters
    def shift(v,p,k):  # multiply by chi_k(U_p): n -> n + k bd(p)
        if p<8: return np.roll(v,k,axis=p)
        return np.roll(v,(-k,)*8,axis=tuple(range(8)))  # bd(P[8]) = -sum_{i<8} bd(P[i]) mod 3
    def Vop(v):
        v=v.reshape(shape)
        for p in range(9):
            v=b[0]*v+b[1]*shift(v,p,1)+b[2]*shift(v,p,-1)
        return v.reshape(-1)
    def Top(v): return Vop(Kd*Vop(v))
    n=3**10
    op=LinearOperator((n,n),matvec=Top,dtype=float)
    w,U=eigsh(op,k=2,which='LA',tol=1e-13)
    o=np.argsort(w)[::-1]; w=w[o]; U=U[:,o]
    lam0,lam1=w; Om=U[:,0]; Om*=np.sign(Om[0])
    def Aop(v):
        v=v.reshape(shape)
        return (0.5*(shift(v,0,1)+shift(v,0,-1))).reshape(-1)
    a=Om@Aop(Om); psi=Aop(Om)-a*Om
    var=psi@psi
    Vpsi=Vop(psi)
    Fmask=np.zeros(len(labels),dtype=np.int64)
    for j,l in enumerate(p0links): Fmask|=((labels[:,l]!=0).astype(np.int64)<<j)
    terms={F:np.sum(Kd[Fmask==F]*Vpsi[Fmask==F]**2)/lam0 for F in range(16)}
    cov=sum(terms.values())
    u=c[1]/c[0]
    return dict(lam0=lam0,lam1=lam1,a=a,var=var,cov=cov,G=terms[15],terms=terms,u=u,gap=math.log(lam0/lam1))
print("beta     tau0 u^4     Cov/(tau0u^4)  <G>/(tau0u^4)  min term    Var/tau0    log(Var/Cov)-4log(1/u)  gap/log(1/u)")
for beta in [0.3,0.1,0.03,0.01,4e-4,1e-4]:
    r=run(beta); u=r['u']; t=0.5*u**4
    print(f"{beta:<8g} {t:.4e}  {r['cov']/t:.6f}      {r['G']/t:.6f}       {min(r['terms'].values()):+.2e}  {r['var']/0.5:.6f}   {math.log(r['var']/r['cov'])-4*math.log(1/u):+.3e}        {r['gap']/math.log(1/u):.4f}")
    if beta<4.1e-4:
        cc=consts(beta,3,1,0)
        print(f"   bounds (N=1, iota=0, d=3): <G>/(tau0u^4) >= L={float(cc['L']):.6f};  Var/tau0 <= V={float(cc['V']):.6f};  D={float(cc['D']):.3e};  |a|={abs(r['a']):.3e}")

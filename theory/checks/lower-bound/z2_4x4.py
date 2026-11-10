# Exact finite-volume gap of Z_2 gauge theory, d=3, spatial torus 4x4 (girth 4), compared with the bounds
#   4 log(1/q_u) <= m <= 4 log(1/u) + D   (q_u = eK e^{2(K+1)theta} u, u = tanh beta, K = 12).
# Physical space: Z_2 cycles of the 4x4 torus graph, coordinates in Z_2^17 (15 plaquettes + 2 windings).
import numpy as np, math, sys
import os; sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','upper-bound'))
from scipy.sparse.linalg import LinearOperator, eigsh
from constants_upper import consts
L=4; nl=2*L*L
def lid(x,y,i): return ((x%L)*L+(y%L))*2+i
def bd(x,y):
    v=np.zeros(nl,dtype=np.int64)
    for l in (lid(x,y,0),lid(x+1,y,1),lid(x,y+1,0),lid(x,y,1)): v[l]^=1
    return v
P=[bd(x,y) for x in range(L) for y in range(L)]
W1=np.zeros(nl,dtype=np.int64); W2=np.zeros(nl,dtype=np.int64)
for x in range(L): W1[lid(x,0,0)]=1
for y in range(L): W2[lid(0,y,1)]=1
gens=P[:15]+[W1,W2]; ng=len(gens); n=2**ng
bits=((np.arange(n)[:,None]>>np.arange(ng)[None,:])&1)
labels=(bits@np.array(gens))%2                 # (n, 32)
popc=labels.sum(1)
idx=np.arange(n)
flip=[idx^(1<<k) for k in range(15)]+[idx^((1<<15)-1)]   # plaquette 16 = sum of the first 15
def run(beta):
    ch,sh=math.cosh(beta),math.sinh(beta); cb,sb=math.cosh(beta/2),math.sinh(beta/2)
    Kd=ch**(nl-popc)*sh**popc
    def Vop(v):
        for fl in flip: v=cb*v+sb*v[fl]
        return v
    op=LinearOperator((n,n),matvec=lambda v: Vop(Kd*Vop(v)),dtype=float)
    w,_=eigsh(op,k=3,which='LA',tol=1e-12)
    w=np.sort(w)[::-1]
    return math.log(w[0]/w[1])
K=12; e=math.e
print("beta      u          m_exact   4log(1/u)  m-4log(1/u)  lower 4log(1/q_u)  upper 4log(1/u)+D")
for beta in [3e-1,1e-1,3e-2,1e-2,1e-3,4e-4,1e-4]:
    u=math.tanh(beta); m=run(beta)
    eps=math.expm1(2*beta); x=e*e*K*eps
    if x<1:
        th=x/(1-x); qu=e*K*math.exp(2*(K+1)*th)*u; low=4*math.log(1/qu)
    else: low=float('nan')
    c=consts(beta,3,1,1) if beta<4.03e-4 else None
    up=4*math.log(1/u)+(float(c['D']) if c else float('nan'))
    print(f"{beta:<8g} {u:.4e}  {m:9.4f}  {4*math.log(1/u):9.4f}  {m-4*math.log(1/u):+.3e}   {low:9.4f}         {up:9.4f}")

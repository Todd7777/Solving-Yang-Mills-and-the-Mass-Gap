# High-precision check: Z_2, d=3, spatial torus 3x3. The transfer matrix preserves the homology (winding) class, so the
# torelon is the top state of the winding-(1,0) sector. Compare its energy log(lam0/lam_tor) with 3log(1/u) and with the
# lower bound 4 log(1/q_u) of Theorem mg:thm:sharp (which assumes spatial sides >= 4).
from mpmath import mp, mpf, cosh, sinh, tanh, log, exp, e, expm1
mp.dps=60
L=3; nl=2*L*L
def lid(x,y,i): return ((x%L)*L+(y%L))*2+i
pl=[]
for x in range(L):
    for y in range(L):
        pl.append(sum(1<<l for l in (lid(x,y,0),lid(x+1,y,1),lid(x,y+1,0),lid(x,y,1))))
W1=sum(1<<lid(x,0,0) for x in range(L)); W2=sum(1<<lid(0,y,1) for y in range(L))
gens=pl[:8]+[W1,W2]
def mask(i):
    m=0
    for k in range(10):
        if i>>k&1: m^=gens[k]
    return m
pop=[bin(mask(i)).count('1') for i in range(1024)]
flips=[1<<k for k in range(8)]+[0xFF]
def top(beta,sector):
    idx=[i for i in range(1024) if (i>>8)==sector]
    pos={i:j for j,i in enumerate(idx)}
    ch,sh=cosh(beta),sinh(beta); cb,sb=cosh(beta/2),sinh(beta/2)
    Kd=[ch**(nl-pop[i])*sh**pop[i] for i in idx]
    def V(v):
        for f in flips:
            v=[cb*v[j]+sb*v[pos[idx[j]^f]] for j in range(len(idx))]
        return v
    v=[mpf(1)]*len(idx); lam=mpf(0)
    for it in range(12):
        w=V([Kd[j]*x for j,x in enumerate(V(v))])
        nrm=max(abs(x) for x in w); lam_new=sum(a*b for a,b in zip(v,w))/sum(a*a for a in v)
        v=[x/nrm for x in w]
        if abs(lam_new-lam)<abs(lam_new)*mpf(10)**-50: break
        lam=lam_new
    return lam_new
K=12
for beta in [mpf('2e-7'),mpf('5e-7'),mpf('1e-6')]:
    u=tanh(beta)
    l0=top(beta,0); lt=top(beta,1)
    m_tor=log(l0/lt)
    eps=expm1(2*beta); x=e**2*K*eps; th=x/(1-x); qu=e*K*exp(2*(K+1)*th)*u
    print(f"beta={float(beta):.1e}: torelon energy {float(m_tor):.4f}, 3log(1/u)={float(3*log(1/u)):.4f}, 4log(1/q_u)={float(4*log(1/qu)):.4f}, 4log(1/u)={float(4*log(1/u)):.4f}")

# New lower bound m >= 4 log(1/q_u), q_u = e K e^{2(K+1)theta} u_*, compared with the old bound 4 log(1/q),
# q = e^3 K (e^{4 beta}-1), Part II's log(1/(6(d-1) tanh beta)), and the upper bound 4 log(1/u) + D.
import sys
import os; sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','upper-bound'))
from mpmath import mp, mpf, e, exp, log, besseli, expm1, quad, pi, cos, tanh, sqrt
from constants_upper import consts, beta0
mp.dps=30
def u_su2(b): return besseli(2,b)/besseli(1,b)          # c_{1/2}/(2 c_0)
def u_u1(b):  return besseli(1,b)/besseli(0,b)
def su3_c(b,rep):
    # Weyl integration over the maximal torus of SU(3): characters of fund (3) and trivial; w = exp((b/3) Re Tr U)
    def integrand(t1,t2):
        z=[mp.expj(t1),mp.expj(t2),mp.expj(-t1-t2)]
        trU=sum(z); w=exp(b/3*mp.re(trU))
        vd=abs((z[0]-z[1])*(z[0]-z[2])*(z[1]-z[2]))**2/6
        chi=1 if rep=='0' else mp.conj(trU)
        return mp.re(w*chi)*vd
    return quad(integrand,[0,2*pi],[0,2*pi])/(4*pi**2)
def u_su3(b):
    c0=su3_c(b,'0'); c3=su3_c(b,'3'); return c3/(3*c0)
groups={'SU(2)':(u_su2,2,1),'U(1)':(u_u1,1,0)}
for d in (3,4):
    K=4*(2*d-3); b0=beta0(d)
    print(f"d={d}: beta0={float(b0):.4e}, 4log(eK)={float(4*log(e*K)):.3f}")
    for name,(uf,N,iota) in list(groups.items())+[('SU(3)',(u_su3,3,0))]:
        for frac in (1e-3,0.5,0.999999):
            b=b0*frac
            eps=expm1(2*b); x=e**2*K*eps; th=x/(1-x)
            u=uf(b)
            qu=e*K*exp(2*(K+1)*th)*u
            q=e**3*K*expm1(4*b)
            D=consts(b,d,N,iota)['D']
            new=4*log(1/qu); old=4*log(1/q); dob=log(1/(6*(d-1)*tanh(b))); up=4*log(1/u)+D
            print(f"  {name:5s} beta={float(b):.3e}: new lower {float(new):7.3f}  old lower {float(old):7.3f}  PartII {float(dob):6.3f}  upper {float(up):7.3f}  window new {float(up-new):6.3f} old {float(up-old):6.3f}")

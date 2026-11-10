# Constants of Theorem mg:thm:fronts: delta = 4 D_e kappa_C(lambda_C) nu(lambda_T), with lambda_T >= eta_C(lambda_C),
# lambda_C >= 4 eta_T(lambda_T).  Lower bound m >= 4log(1/u) - 8(K+1)theta - log(1+delta), compared with Theorem sharp.
import sys, os
import os; sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','upper-bound'))
from mpmath import mp, mpf, e, exp, log, besseli, expm1
from constants_upper import consts, beta0
mp.dps=30
def Nk(k,d):
    if k==4: return 2*d-3
    if k==5: return mpf(0) if False else e*(e*(4*d-6))**(k-1)   # crude bound (N_5 can be nonzero via torelons)
    return e*(e*(4*d-6))**(k-1)
def constants(beta,d,u):
    K=4*(2*d-3); De=4*d-5; KS=4*(2*d-5)
    eps=expm1(2*beta); x=e**2*K*eps; th=x/(1-x); lam=exp(2*(K+1)*th)
    a=lam*u; b=lam*expm1(4*beta)
    def etaC(lC):
        r=e*KS*b*exp(lC); return 4*(d-2)*De*e*b*exp(lC)/(1-r) if r<1 else mp.inf
    def kappaC(lC):
        r=e*KS*b*exp(lC); return 16*(d-2)*De*e*b*exp(lC)/(1-r)**2 if r<1 else mp.inf
    def etaT(lT): return 2*De*sum(Nk(k,d)*exp(lT*k)*a**k/(1-a**k) for k in range(4,60))
    def nu(lT): return sum(k*Nk(k,d)*exp(lT*k)*a**(k-4) for k in range(4,60))
    lC=mpf(1); lT=etaC(lC)
    for _ in range(50):
        lC=max(4*etaT(lT),mpf(10)**-30); lT=etaC(lC)
    assert etaC(lC)<=lT*(1+1e-20) and 4*etaT(lT)<=lC*(1+1e-20)
    delta=4*De*kappaC(lC)*nu(lT)
    return dict(K=K,th=th,a=a,b=b,lT=lT,lC=lC,delta=delta)
def u_su2(b): return besseli(2,b)/besseli(1,b)
for d in (3,4):
    b0=beta0(d); K=4*(2*d-3)
    print(f"d={d}, beta0={float(b0):.3e}, 4log(eK)={float(4*log(e*K)):.2f}")
    for beta in [mpf('1e-8'),mpf('1e-7'),mpf('1e-6'),mpf('1e-5'),b0/10,b0/2,b0*0.999999]:
        u=u_su2(beta); c=constants(beta,d,u)
        loss_fronts=8*(K+1)*c['th']+log(1+c['delta']); loss_sharp=4*log(e*K)+8*(K+1)*c['th']
        D=consts(beta,d,2,1)['D']
        print(f"  beta={float(beta):.2e}: delta={float(c['delta']):.3e} (delta/beta={float(c['delta']/beta):.3e}), lambda_T={float(c['lT']):.2e}, "
              f"lower-bound loss: fronts {float(loss_fronts):.3f}, Thm sharp {float(loss_sharp):.3f}; upper D={float(D):.2e}")

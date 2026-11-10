# Constants of the improved upper bound (reflection-positivity route).
import numpy as np
from mpmath import mp, mpf, e, exp, log, sqrt, expm1
mp.dps = 40

def consts(beta, d, N, iota):
    beta = mpf(beta)
    K = 4*(2*d-3)
    eps = expm1(2*beta)
    x = e**2*K*eps
    th = x/(1-x)
    q1 = e*K*exp((K+1)*th)*eps
    tau0 = mpf(1+iota)/2
    TQ = (1+eps)**6*exp((24*d-42)*th)*(1+mpf(1)/K)**(24*(d-2))*q1**5/(1-q1)
    Tp = (1+eps)*exp((8*d-11)*th)*(1+mpf(1)/K)**K*q1**5/(1-q1)
    ubar = ((1+iota)*beta+beta**2*exp(beta))/(2*N**2)
    astar = N*ubar*(1+Tp)+sqrt(tau0)*Tp
    eta = eps/sqrt(tau0)+astar*N*ubar/tau0
    L = (1-eta)**2*(1-TQ)/(1+eps**6) - (1+astar**2/tau0)*exp(16*beta)*(mpf(N)**4/tau0**2)*TQ
    V = (1+eps)*(1+Tp)+Tp
    D = log(V/L) if L > 0 else None
    return dict(K=K, eps=eps, x=x, th=th, q1=q1, TQ=TQ, Tp=Tp, eta=eta, L=L, V=V, D=D)

def beta0(d):
    K = 4*(2*d-3)
    return log(1+1/(e**2*K*(K+2)))/2

cases = [("SU(2)",2,1),("SU(3)",3,0),("U(1)",1,0)]
if __name__ == "__main__":
 for d in (3,4,5):
     b0 = beta0(d)
     print(f"d={d}: beta0={float(b0):.6e}")
     for name,N,iota in cases:
         for frac in (1e-3, 0.1, 0.5, 0.9, 0.999999):
             c = consts(b0*frac, d, N, iota)
             print(f"  {name:6s} beta={float(b0*frac):.3e}: TQ={float(c['TQ']):.3e} Tp={float(c['Tp']):.3e} eta={float(c['eta']):.3e} L={float(c['L']):.6f} V={float(c['V']):.6f} D={float(c['D']):.3e}  D/beta={float(c['D']/(b0*frac)):.2f}")
     # largest N^4/tau0^2 allowed at beta -> beta0 (L>0): solve linear in that factor
     c = consts(b0*(1-mpf(10)**-12), d, 2, 1)
     first = (1-c['eta'])**2*(1-c['TQ'])/(1+c['eps']**6)
     print(f"  at beta0: TQ={float(c['TQ']):.4e}; L>0 iff (1+a*^2/tau0) e^(16b) N^4/tau0^2 TQ < {float(first):.6f}  => N^4/tau0^2 < ~{float(first/(exp(16*b0)*c['TQ'])):.0f}")

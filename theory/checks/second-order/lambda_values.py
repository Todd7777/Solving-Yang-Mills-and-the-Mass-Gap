# Section 10 of Part IV: Lambda(beta) = max over C-sectors of Var_{w/c0}(phi), phi a real unit plaquette function in span(chi, chi-bar)
# (the largest eigenvalue of T = P_rho (M_{1+f} - |f><f|) P_rho). Weyl integration on the maximal torus (trapezoidal rule).
# Prints -log Lambda_C and its expansion coefficients for U(1), SU(2), ..., SU(5) (fundamental representation).
import numpy as np, math
def torus(N, M):
    th=np.linspace(0,2*np.pi,M,endpoint=False)
    g=np.meshgrid(*([th]*(N-1)),indexing='ij'); ang=[x.ravel() for x in g]; ang.append(-sum(ang))
    z=np.exp(1j*np.array(ang)); vd=np.ones(z.shape[1])
    for i in range(N):
        for j in range(i+1,N): vd*=np.abs(z[i]-z[j])**2
    return z, vd/(math.factorial(N)*M**(N-1))
GRID={1:4096,2:256,3:64,4:24,5:14}
def lam(N,beta):
    if N==1:
        th=np.linspace(0,2*np.pi,GRID[1],endpoint=False); X=np.exp(1j*th); w=np.ones_like(th)/len(th)
    else:
        z,w=torus(N,GRID[N]); X=z.sum(0)
    W=np.exp(beta/N*np.real(X)); p=w*W/np.sum(W*w)
    u=np.sum(p*np.conj(X)).real/N
    phis={'+':X.real} if N==2 else {'+':math.sqrt(2)*X.real,'-':math.sqrt(2)*X.imag}
    return u,{k:np.sum(p*ph**2)-np.sum(p*ph)**2 for k,ph in phis.items()}
if __name__=="__main__":
    print("check of the quadrature: int 1 and int |chi|^2 for SU(N):",
          [round(float(np.sum(torus(N,GRID[N])[1])),12) for N in (2,3,4,5)],
          [round(float(np.sum(np.abs(torus(N,GRID[N])[0].sum(0))**2*torus(N,GRID[N])[1])),12) for N in (2,3,4,5)])
    for N in (1,2,3,4,5):
        name='U(1)' if N==1 else f'SU({N})'
        for beta in (4e-3,2e-3,1e-3):
            u,L=lam(N,beta)
            parts=[]
            for k,v in L.items():
                g=-math.log(v)
                lin={3:-1/6}.get(N,0.0)*(1 if k=='+' else -1)   # first-order term -beta/(N tau0) int(Re chi)^3 (+/- for C-even/odd)
                parts.append(f"C{k}: -log Lam = {g:+.6e}, (-log Lam - ({lin:+.4f})beta)/beta^2 = {(g-lin*beta)/beta**2:+.6f}")
            print(f"{name:6s} beta={beta:.0e}: u={u:.6e}  "+'   '.join(parts))

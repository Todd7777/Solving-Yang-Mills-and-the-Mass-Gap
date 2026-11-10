# Enumerate the local counting constants used in the planned Section 7, on spatial tori of dimension d-1 = 2, 3
# with all sides >= 4:  D_e (edges equal to or sharing a vertex with a given edge), N_4 (fronts with 4 edges containing a
# given edge: connected, every vertex degree >= 2), N_k for k = 5,6,7 versus (eD)^{k-1} with D = 4d-6, the number of
# plaquettes of the slice containing a link (2(d-2)), and K_S (plaquettes of the slice sharing a link with a given one).
import itertools, math
from collections import defaultdict
def torus(shape):
    n=len(shape); sites=list(itertools.product(*[range(L) for L in shape]))
    def add(x,i,s=1):
        y=list(x); y[i]=(y[i]+s)%shape[i]; return tuple(y)
    edges=[(x,i) for x in sites for i in range(n)]
    def ends(e): x,i=e; return (x,add(x,i))
    plaqs=[frozenset([(x,i),(add(x,i),j),(add(x,j),i),(x,j)]) for x in sites for i in range(n) for j in range(i+1,n)]
    return sites,edges,ends,plaqs
for dsp,shape in [(2,(4,4)),(2,(5,6)),(3,(4,4,4)),(3,(4,5,5))]:
    d=dsp+1
    sites,edges,ends,plaqs=torus(shape)
    vert_edges=defaultdict(set)
    for e in edges:
        for v in ends(e): vert_edges[v].add(e)
    adj={e:set().union(*[vert_edges[v] for v in ends(e)])-{e} for e in edges}
    De=max(len(adj[e])+1 for e in edges); D=max(len(adj[e]) for e in edges)
    e0=edges[0]
    # connected edge sets containing e0 of size k (lattice animals of edges), then filter min-degree>=2
    level={frozenset([e0])}; animals={1:level}
    for k in range(2,8):
        nxt=set()
        for A in level:
            for e in set().union(*[adj[a] for a in A])-A: nxt.add(A|{e})
        animals[k]=nxt; level=nxt
    def mindeg2(A):
        deg=defaultdict(int)
        for e in A:
            for v in ends(e): deg[v]+=1
        return all(c>=2 for c in deg.values())
    Nk={k:sum(1 for A in animals[k] if mindeg2(A)) for k in range(4,8)}
    pl_link=defaultdict(int)
    for P in plaqs:
        for e in P: pl_link[e]+=1
    KS=max(len({Q for Q in plaqs if Q!=P and P&Q}) for P in plaqs[:5])
    print(f"d={d} slice {shape}: D_e={De} (4d-5={4*d-5}), D={D} (4d-6={4*d-6}), N_4={Nk[4]} (2d-3={2*d-3}), "
          f"N_5={Nk[5]}, N_6={Nk[6]}, N_7={Nk[7]}; animals_k=[{', '.join(str(len(animals[k])) for k in range(4,8))}] vs (eD)^(k-1)=[{', '.join(f'{(math.e*D)**(k-1):.0f}' for k in range(4,8))}]; "
          f"plaquettes per link={max(pl_link.values())} (2(d-2)={2*(d-2)}), K_S={KS} (4(2d-5)={4*(2*d-5)})")

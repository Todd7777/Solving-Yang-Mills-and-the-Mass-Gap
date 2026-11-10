# Check of the structural lemmas for the front/segment decomposition (planned Section 7):
#  (1) fronts: vertices of T_s have degree 0 or >=2;  (2) column lemma: T_s ^ T_{s-1} subset E(S_s) (+ support links);
#  (3) segments: every front component lies in exactly one segment; each segment is touched at both ends;
#  (4) adjacency in X implies connectivity in the incidence graph;  (5) a path S_A -> S_B covers all inner layers;
#  (6) |X| = sum_segments |phi| l + sum_clusters |C|.
import random, itertools, sys
from collections import defaultdict, deque
L0,L1,L2=8,5,5; T=3            # time extent, spatial sides; B at slice T
def nb(site,mu,s=1):
    x=list(site); x[mu]=(x[mu]+s)%(L0,L1,L2)[mu]; return tuple(x)
def plaq_links(site,mu,nu):
    return frozenset([(site,mu),(nb(site,mu),nu),(nb(site,nu),mu),(site,nu)])
sites=list(itertools.product(range(L0),range(L1),range(L2)))
PLAQ={(x,mu,nu):plaq_links(x,mu,nu) for x in sites for mu in range(3) for nu in range(mu+1,3)}
def cube(x):  # 6 faces of the unit cube based at x
    f=[]
    for mu,nu in ((0,1),(0,2),(1,2)):
        lam=3-mu-nu
        f+= [(x,mu,nu),(nb(x,lam),mu,nu)]
    return set(f)
def tube(xs,t0,t1):  # straight tube over spatial plaquette at xs (spatial coords), layers t0..t1-1
    out=set()
    for t in range(t0,t1):
        x=(t%L0,)+xs
        for (y,i) in [(x,1),(nb(x,1),2),(nb(x,2),1),(x,2)]:
            out.add((y,0,i))
    return out
p0=(0,0,0)
SA=set(PLAQ[(p0,1,2)]); SB=set(PLAQ[((T,0,0),1,2)]); S=SA|SB
def free_ok(X):
    deg=defaultdict(int)
    for p in X:
        for l in PLAQ[p]: deg[l]+=1
    return all(v!=1 for l,v in deg.items() if l not in S)
def random_X(rng):
    X=set(tube((0,0),0,T))
    # optional shifted pieces: replace a layer of the tube by a detour via the neighbouring plaquette
    for _ in range(rng.randint(0,4)):
        cnt=defaultdict(int)
        for _ in range(rng.randint(1,3)):
            x=(rng.randrange(-2,T+2)%L0,rng.randrange(-1,2)%L1,rng.randrange(-1,2)%L2)
            for f in cube(x): cnt[f]^=1
        piece={f for f,v in cnt.items() if v}
        if piece and rng.random()<0.6: X|=piece          # union: branching allowed
        elif piece: X^=piece                              # mod-2 sum: deforms the tube
    if rng.random()<0.3:
        X|=tube((rng.randrange(-1,2)%L1,rng.randrange(-1,2)%L2),rng.randrange(0,T),T)|cube((T,0,0))
    return X
def analyse(X):
    temporal=[p for p in X if p[1]==0]; spatial=[p for p in X if p[1]!=0]
    Ts=defaultdict(set)                       # layer -> set of spatial links (base links of temporal plaquettes)
    for (x,mu,nu) in temporal: Ts[x[0]].add(((x[1],x[2]),nu))
    Ss=defaultdict(set)
    for p in spatial: Ss[p[0][0]].add(p)
    def ends(e):
        (a,b),i=e; v=(a,b); w=((a+(i==1))%L1,(b+(i==2))%L2); return v,w
    # (1) degrees
    for s,E in Ts.items():
        deg=defaultdict(int)
        for e in E:
            for v in ends(e): deg[v]+=1
        assert all(c>=2 for c in deg.values()),"front degree"
    # clusters
    clusters=[]
    for s,P in Ss.items():
        P=set(P)
        while P:
            q=P.pop(); comp={q}; st=[q]
            while st:
                a=st.pop()
                for b in list(P):
                    if PLAQ[a]&PLAQ[b]: P.remove(b); comp.add(b); st.append(b)
            clusters.append((s,frozenset(comp)))
    def clinks(cl):
        s,comp=cl; out=set()
        for (x,mu,nu) in comp:
            for (y,lam) in PLAQ[(x,mu,nu)]: out.add(((y[1],y[2]),lam))
        return out
    CL={cl:clinks(cl) for cl in clusters}
    # (2) column lemma
    for s in range(L0):
        diff=Ts[s]^Ts[(s-1)%L0]
        allowed=set().union(*[CL[cl] for cl in clusters if cl[0]==s]) if clusters else set()
        if s==0: allowed|={((y[1],y[2]),l) for (y,l) in SA}
        if s==T: allowed|={((y[1],y[2]),l) for (y,l) in SB}
        assert diff<=allowed,"column lemma"
    # front components per layer
    def comps(E):
        E=set(E); out=[]
        while E:
            e=E.pop(); comp={e}; st=[e]
            while st:
                a=st.pop(); va=set(ends(a))
                for b in list(E):
                    if va & set(ends(b)): E.remove(b); comp.add(b); st.append(b)
            out.append(frozenset(comp))
        return out
    FC={s:comps(Ts[s]) for s in range(L0)}
    def touches(linkset,phi):
        V=set(v for e in phi for v in ends(e))
        return any(l in phi or (set(ends(l))&V) for l in linkset)
    supp={0:{((y[1],y[2]),l) for (y,l) in SA},T:{((y[1],y[2]),l) for (y,l) in SB}}
    def touched_at(slice_,phi):
        ts=[('C',cl) for cl in clusters if cl[0]==slice_ and touches(CL[cl],phi)]
        if slice_ in supp and touches(supp[slice_],phi): ts.append(('S',slice_))
        return ts
    # (3) segments: chain identical untouched components
    nodes={}; segs=[]
    seen=set()
    for s in range(L0):
        for phi in FC[s]:
            if (s,phi) in seen: continue
            # go down to the start
            s1=s
            while True:
                prev=(s1-1)%L0
                if phi in FC[prev] and not touched_at(s1,phi) and (prev,phi) not in seen and s1!=s+1-L0:
                    s1=prev
                    if s1==s: break
                else: break
            s2=s1; layers=[]
            while True:
                layers.append(s2); seen.add((s2,phi))
                nxt=(s2+1)%L0
                if phi in FC[nxt] and not touched_at(nxt,phi) and (nxt,phi) not in seen: s2=nxt
                else: break
            bottom=touched_at(layers[0],phi); top=touched_at((layers[-1]+1)%L0,phi)
            assert bottom and top, ("segment end untouched",layers,bottom,top)
            segs.append((phi,tuple(layers),bottom,top))
    # (6) weights
    assert sum(len(phi)*len(ly) for phi,ly,_,_ in segs)+sum(len(c[1]) for c in clusters)==len(X),"weight identity"
    # incidence graph
    G=defaultdict(set)
    for k,(phi,ly,bot,top) in enumerate(segs):
        for t in bot+top: G[('seg',k)].add(t); G[t].add(('seg',k))
    for cl in clusters:
        for sl in supp:
            if cl[0]==sl and CL[cl]&supp[sl]: G[('C',cl)].add(('S',sl)); G[('S',sl)].add(('C',cl))
    # map plaquettes to nodes
    node_of={}
    for k,(phi,ly,_,_) in enumerate(segs):
        for s in ly:
            for ((a,b),i) in phi: node_of[((s,a,b),0,i)]=('seg',k)
    for cl in clusters:
        for p in cl[1]: node_of[p]=('C',cl)
    assert set(node_of)==set(X)
    def comp_of(start):
        seen2={start}; dq=deque([start])
        while dq:
            a=dq.popleft()
            for b in G[a]:
                if b not in seen2: seen2.add(b); dq.append(b)
        return seen2
    # (4) adjacency in X -> same I-component
    for p in X:
        cp=comp_of(node_of[p])
        for q in X:
            if q!=p and PLAQ[p]&PLAQ[q]: assert node_of[q] in cp,"adjacency not reflected"
        # plaquettes containing support links are linked to the support node
        if PLAQ[p]&SA: assert ('S',0) in cp,"support A"
        if PLAQ[p]&SB: assert ('S',T) in cp,"support B"
    # (5) path covering
    start=('S',0); goal=('S',T)
    prev={start:None}; dq=deque([start])
    while dq:
        a=dq.popleft()
        for b in G[a]:
            if b not in prev: prev[b]=a; dq.append(b)
    if goal in prev:
        path=[]; a=goal
        while a: path.append(a); a=prev[a]
        covered=set()
        for nd in path:
            if nd[0]=='seg': covered|=set(segs[nd[1]][1])
        assert set(range(T))<=covered, ("path does not cover inner layers",covered)
        return True,len(segs),len(clusters)
    return False,len(segs),len(clusters)
rng=random.Random(7); n_ok=0; n_conn=0; tried=0
while n_ok<400 and tried<20000:
    tried+=1
    X=random_X(rng)
    if not X or not free_ok(X): continue
    if any(all(l[0][0]!=s for p in X if p[1]==0 for l in []) for s in []): pass
    layers_hit=set(p[0][0] for p in X if p[1]==0)
    if len(layers_hit)==L0: continue        # skip wrapping configurations
    conn,ns,nc=analyse(X); n_ok+=1; n_conn+=conn
print(f"checked {n_ok} random admissible X (non-wrapping), {n_conn} of them connecting S_A to S_B: all structural statements hold")

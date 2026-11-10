# Exhaustive check of the five-plaquette lemma: if every plaquette of Y has at most one link in B and every
# link of E(Y)\B lies in >= 2 plaquettes of Y, then |Y| >= 5.  B = E(Q) (unit cube, directions 0,1,2) or
# B = E(p0).  Connected Y suffice (components inherit the hypothesis).
import itertools, sys, time
from collections import Counter
def torus(shape):
    d=len(shape)
    sites=list(itertools.product(*[range(L) for L in shape]))
    def add(x,mu,s=1):
        y=list(x); y[mu]=(y[mu]+s)%shape[mu]; return tuple(y)
    def link(x,mu): return (x,mu)
    plaqs=[]
    for x in sites:
        for mu in range(d):
            for nu in range(mu+1,d):
                plaqs.append(frozenset([link(x,mu),link(add(x,mu),nu),link(add(x,nu),mu),link(x,nu)]))
    return d,sites,add,link,plaqs
def check(shape,which,maxsize):
    d,sites,add,link,plaqs=torus(shape)
    o=tuple([0]*d)
    if which=="cube":
        cube_sites=[tuple(a if i<3 else 0 for i,a in enumerate(v)) for v in itertools.product((0,1),repeat=3) for _ in [0]]
        cube_sites=[tuple(list(v)+[0]*(d-3)) for v in itertools.product((0,1),repeat=3)]
        B=set()
        for x in cube_sites:
            for mu in range(3):
                if x[mu]==0: B.add(link(x,mu))
        assert len(B)==12
    else:
        B=set(plaqs[[i for i,p in enumerate(plaqs) if link(o,0) in p and link(o,1) in p][0]])
        assert len(B)==4
    allowed=[i for i,p in enumerate(plaqs) if len(p&B)<=1]
    excluded=len(plaqs)-len(allowed)
    lk2p={}
    for i in allowed:
        for l in plaqs[i]: lk2p.setdefault(l,[]).append(i)
    nbr={i:set() for i in allowed}
    for l,ps in lk2p.items():
        for i in ps:
            nbr[i].update(ps)
    for i in allowed: nbr[i].discard(i)
    def valid(Y):
        cnt=Counter(l for i in Y for l in plaqs[i])
        return all(c>=2 for l,c in cnt.items() if l not in B)
    found={}
    level={frozenset([i]) for i in allowed}
    for size in range(1,maxsize+1):
        good=[Y for Y in level if valid(Y)]
        found[size]=good
        if size==maxsize: break
        nxt=set()
        for Y in level:
            ext=set().union(*(nbr[i] for i in Y))-Y
            for j in ext: nxt.add(Y|{j})
        level=nxt
    return excluded,{s:len(v) for s,v in found.items()},found
for shape,which,maxsize in [((3,3,3),"cube",5),((4,4,4),"cube",5),((3,3,3),"plaq",5),((4,4,4),"plaq",5),
                            ((3,3,3,3),"cube",4),((3,3,3,3),"plaq",4),((4,3,3,3),"plaq",4)]:
    t=time.time()
    exc,cnt,found=check(shape,which,maxsize)
    msg=f"{shape} B={which}: excluded plaquettes={exc}, valid connected Y by size={cnt}"
    if maxsize==5 and cnt[5]:
        # describe size-5 solutions: are they cube caps (5 faces of a unit cube)?
        msg+=f"  [size-5 solutions: {cnt[5]}]"
    print(msg, f"({time.time()-t:.0f}s)"); sys.stdout.flush()

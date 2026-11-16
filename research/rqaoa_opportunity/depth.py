"""11/17 probe: does deeper QAOA (p=1,2,3) inside RQAOA fix the 'only sees neighbours' problem?

Full statevector simulation (small paint-shop instances so the exact optimum is known).
At every recursion step QAOA_p is re-optimised (INTERP warm start from p-1 plus 3 random starts),
the strongest |<Z_i Z_j>| over all pairs is eliminated, until 4 variables remain (then brute force).
Diagnostic: after each elimination, is the exact optimum of the reduced problem still the original optimum?
(i.e. was an optimal solution 'locked out'). Fresh seeds 800+, nothing tuned on them.
"""
import json,time,itertools
import numpy as np
from scipy.optimize import minimize
from sweep import OUT,make_instance,score,paint_cost,dumpcsv
from transfer import reduce
from paint_baselines import author_solve,AUTHOR_RSG

def all_spins(n):
    b=((np.arange(1<<n)[:,None]>>np.arange(n))&1).astype(np.int8);return 1-2*b   # rows: states, cols: qubit q

def costs(K,Z):return np.einsum('si,ij,sj->s',Z,K,Z,optimize=True)/2

def evolve(c,params,n):
    p=len(params)//2;v=np.full(len(c),1/np.sqrt(len(c)),complex)
    for ga,be in zip(params[:p],params[p:]):
        v=v*np.exp(-1j*ga*c);co=np.cos(be);si=-1j*np.sin(be)
        for q in range(n):
            blk=v.reshape(-1,2,1<<q);lo=blk[:,0,:].copy();hi=blk[:,1,:];blk[:,0,:]=co*lo+si*hi;blk[:,1,:]=si*lo+co*hi
    return v

def optimise(c,n,p,warm,rng):
    scale=1/max(np.abs(c).max(),1e-12)
    f=lambda x:-float(np.real(np.vdot(v:=evolve(c,x,n),c*v)))
    starts=[] if warm is None else [warm]
    starts+=[np.concatenate([rng.uniform(0,2*scale*np.pi/2,p),rng.uniform(0,np.pi/2,p)]) for _ in range(3)]
    best=None
    for x0 in starts:
        r=minimize(f,x0,method='L-BFGS-B',options=dict(maxiter=200))
        if best is None or r.fun<best.fun:best=r
    return best.x,-best.fun

def interp(x):
    p=len(x)//2;g,b=x[:p],x[p:];q=p+1
    up=lambda a:np.array([((i)/p)*(a[i-1] if i>0 else 0)+((p-i)/p)*(a[i] if i<p else 0) for i in range(q)])
    return np.concatenate([up(g),up(b)])

def qaoa_params(c,n,p,rng):
    x=None
    for d in range(1,p+1):
        x,e=optimise(c,n,d,None if x is None else interp(x),rng)
    return x,e

def rqaoa(K0,p,rng):
    n0=len(K0);K=K0.copy();active=list(range(n0));trace=[];energies=[]
    while len(K)>4:
        n=len(K);Z=all_spins(n);c=costs(K,Z);x,e=qaoa_params(c,n,p,rng);v=evolve(c,x,n);pr=np.abs(v)**2
        M=np.einsum('s,si,sj->ij',pr,Z,Z,optimize=True);i,j=np.triu_indices(n,1);k=int(np.argmax(np.abs(M[i,j])))
        u,w=int(i[k]),int(j[k]);s=1 if M[u,w]>=0 else -1
        energies.append(dict(remaining=n,expect=e,opt=float(c.max()),corr=float(abs(M[u,w]))))
        trace.append((active[u],active[w],s));K=reduce(K,u,w,s);active.pop(w)
    Z=all_spins(len(K));zb=Z[int(np.argmax(costs(K,Z)))];sol=dict(zip(active,zb))
    for keep,rem,s in reversed(trace):sol[rem]=s*sol[keep]
    return np.array([sol[i] for i in range(n0)]),trace,energies

def locked_step(K0,trace,opt):
    """first elimination step after which the original optimum is no longer reachable (None if never)."""
    n0=len(K0)
    for t in range(1,len(trace)+1):
        K=K0.copy();active=list(range(n0))
        for keep,rem,s in trace[:t]:
            u=active.index(keep);w=active.index(rem);K=reduce(K,u,w,s);active.pop(w)
        Z=all_spins(len(K));zb=Z[int(np.argmax(costs(K,Z)))];sol=dict(zip(active,zb))
        for keep,rem,s in reversed(trace[:t]):sol[rem]=s*sol[keep]
        z=np.array([sol[i] for i in range(n0)])
        if score(K0,z)<opt-1e-9:return t
    return None

def run(sizes,ps=(1,2,3)):
    rows=[]
    for n,seeds in sizes:
        for seed in seeds:
            K,meta=make_instance('paintshop',n,seed);Z=all_spins(n);c=costs(K,Z);opt=float(c.max())
            zopt=Z[int(np.argmax(c))];bestc=paint_cost(meta,zopt);rsg=paint_cost(meta,author_solve(meta,AUTHOR_RSG))
            line=[n,seed,'opt',bestc,'RSG',rsg]
            for p in ps:
                st=time.perf_counter();z,tr,en=rqaoa(K,p,np.random.default_rng(seed*7+p));t=time.perf_counter()-st
                cc=paint_cost(meta,z);ls=locked_step(K,tr,opt)
                ratio0=(en[0]['expect']-c.min())/(en[0]['opt']-c.min())
                rows.append(dict(n=n,seed=seed,p=p,rqaoa=cc,optimum=bestc,RSG=rsg,locked_step=ls if ls else '',steps=len(tr),qaoa_ratio_first=round(ratio0,4),first_corr=round(en[0]['corr'],4),seconds=round(t,2)))
                line+=[f'p{p}',cc,f'lock{ls}',round(ratio0,3)]
            print(*line,flush=True);dumpcsv(OUT/'depth.csv',rows)

if __name__=='__main__':
    run([(12,range(800,830)),(16,range(900,920))])

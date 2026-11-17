"""11/17 signal comparison: same RQAOA skeleton, only the 'which pair to lock' signal changes.

Instances: open-boundary LxW square lattice, random +-J couplings (planar, frustrated, no field).
Objective: maximise sum J_ij s_i s_j  (= minimise Ising energy). K_ij = J_ij.
Skeleton: at each step pick an EDGE of the current (contracted) graph with the largest |signal|,
lock s_v = sign * s_u, contract (planarity preserved), stop at 4 spins, brute force, expand.
Lock check: after every step, exhaustive optimum of the reduced problem is expanded and compared
with the original optimum -> first step at which every optimal solution was excluded.

Signals
  qaoa_p{1,2,3}   : <Z_u Z_v> of the optimised QAOA_p state (full statevector), radius p
  local_r{1,2,3}  : exact Boltzmann <s_u s_v> (beta=2) on the subgraph within graph distance r of the edge
                    -> classical signal with the SAME observation radius as qaoa_p{r}
  global_thermal  : exact Boltzmann <s_u s_v> (beta=2) on the whole graph (stands in for a converged
                    Monte Carlo; global, reference only)
  sdp             : Goemans-Williamson SDP Gram entry X_uv (global, reference only)
  random          : random edge, sign of its coupling
beta=2 and all settings fixed before any run; fresh seeds 2000+.
"""
import sys,time,json,csv
from collections import deque
import numpy as np
from sweep import OUT,dumpcsv
from transfer import reduce
from depth import all_spins,costs,evolve,qaoa_params

BETA=2.0

def lattice(L,W,seed,couplings='pm'):
    """couplings: 'pm' = random +-1 (highly degenerate ground states);
    'gauss' = standard normal (unique ground state up to global flip, so wrong locks are easier to see)."""
    rng=np.random.default_rng(seed);n=L*W;K=np.zeros((n,n))
    for x in range(L):
        for y in range(W):
            i=x*W+y
            for j in ([i+W] if x+1<L else [])+([i+1] if y+1<W else []):
                K[i,j]=K[j,i]=rng.choice([-1.,1.]) if couplings=='pm' else rng.normal()
    return K

def edges(K):
    i,j=np.nonzero(np.triu(K,1));return list(zip(i.tolist(),j.tolist()))

def neighbourhood(K,u,v,r):
    adj=[np.flatnonzero(K[a]) for a in range(len(K))];dist={u:0,v:0};q=deque([u,v])
    while q:
        a=q.popleft()
        if dist[a]>=r:continue
        for b in adj[a]:
            if b not in dist:dist[b]=dist[a]+1;q.append(b)
    return sorted(dist)

def boltz_corr(Ksub,a,b,beta):
    Z=all_spins(len(Ksub));e=costs(Ksub,Z);w=np.exp(beta*(e-e.max()));return float(np.sum(w*Z[:,a]*Z[:,b])/w.sum())

def signal(name,K,rng,cache):
    E=edges(K);n=len(K)
    if name=='random':
        u,v=E[rng.integers(len(E))];return {(u,v):np.sign(K[u,v])*1e-3}
    if name.startswith('qaoa_p'):
        p=int(name[-1]);Z=all_spins(n);c=costs(K,Z);x,_=qaoa_params(c,n,p,rng);pr=np.abs(evolve(c,x,n))**2
        return {(u,v):float(np.sum(pr*Z[:,u]*Z[:,v])) for u,v in E}
    if name.startswith('local_r'):
        r=int(name[-1]);out={}
        for u,v in E:
            nb=neighbourhood(K,u,v,r);sub=K[np.ix_(nb,nb)];out[(u,v)]=boltz_corr(sub,nb.index(u),nb.index(v),BETA)
        return out
    if name=='global_thermal':
        Z=all_spins(n);e=costs(K,Z);w=np.exp(BETA*(e-e.max()));w/=w.sum()
        return {(u,v):float(np.sum(w*Z[:,u]*Z[:,v])) for u,v in E}
    if name=='sdp':
        import cvxpy as cp
        X=cp.Variable((n,n),PSD=True);cp.Problem(cp.Maximize(cp.sum(cp.multiply(K,X))/2),[cp.diag(X)==1]).solve(solver=cp.SCS,eps=1e-6)
        Xv=X.value;return {(u,v):float(Xv[u,v]) for u,v in E}
    raise ValueError(name)

def best_expand(K,active,trace,n0):
    Z=all_spins(len(K));zb=Z[int(np.argmax(costs(K,Z)))];sol=dict(zip(active,zb))
    for keep,rem,s in reversed(trace):sol[rem]=s*sol[keep]
    return np.array([sol[i] for i in range(n0)])

def run_rqaoa(K0,name,seed):
    n0=len(K0);K=K0.copy();active=list(range(n0));trace=[];rng=np.random.default_rng(seed)
    Z0=all_spins(n0);opt=float(costs(K0,Z0).max());first_lock=None;st=time.perf_counter();sig_time=0.
    while len(K)>4 and len(edges(K))>0:
        t=time.perf_counter();S=signal(name,K,rng,None);sig_time+=time.perf_counter()-t
        mx=max(abs(v) for v in S.values());cands=[e for e,v in S.items() if abs(v)>=mx-1e-9]
        u,v=cands[rng.integers(len(cands))];s=1 if S[(u,v)]>0 else (-1 if S[(u,v)]<0 else int(np.sign(K[u,v])) or 1)
        trace.append((active[u],active[v],s));K=reduce(K,u,v,s);active.pop(v)
        z=best_expand(K,active,trace,n0);zz=z@K0@z/2
        if first_lock is None and zz<opt-1e-9:first_lock=len(trace)
    z=best_expand(K,active,trace,n0);final=float(z@K0@z/2)
    return dict(locked=int(first_lock is not None),first_lock=first_lock or '',steps=len(trace),gap=round(opt-final,6),
                signal_seconds=round(sig_time,3),total_seconds=round(time.perf_counter()-st,3))

def run(sizes,couplings='pm',outname=None):
    """Resumable: rows already in the output CSV are skipped."""
    path=OUT/(outname or f'signals_{couplings}.csv');rows=[]
    if path.exists():rows=list(csv.DictReader(open(path,encoding='utf-8-sig')))
    done={(int(r['n']),int(r['seed']),r['signal']) for r in rows}
    for (L,W),seeds,sig_names in sizes:
        for seed in seeds:
            K=lattice(L,W,seed,couplings)
            for name in sig_names:
                if (L*W,seed,name) in done:continue
                r=run_rqaoa(K,name,seed*31+len(name));r.update(L=L,W=W,n=L*W,seed=seed,signal=name,couplings=couplings);rows.append(r)
                dumpcsv(path,rows)
            print(couplings,L,W,seed,{r['signal']:(r['locked'],r['gap']) for r in rows if int(r['seed'])==seed and int(r['n'])==L*W},flush=True)

BASE=['random','local_r1','local_r2','local_r3','global_thermal','sdp','qaoa_p1','qaoa_p2']
PLANS={
    # quick check (~10 min): does everything run on this machine?
    'quick':[((4,4),range(2000,2003),BASE+['qaoa_p3'])],
    # main study: 4x4 with p<=3 (30 lattices), 4x5 with p<=2 (20 lattices)
    'full':[((4,4),range(2000,2030),BASE+['qaoa_p3']),((4,5),range(3000,3020),BASE)],
}

if __name__=='__main__':
    plan=sys.argv[1] if len(sys.argv)>1 else 'full'
    couplings=sys.argv[2] if len(sys.argv)>2 else 'gauss'
    run(PLANS[plan],couplings,f'signals_{couplings}_{plan}.csv')

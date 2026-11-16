"""11/16 probe: why p=1 RQAOA loses on binary paint shop, and whether published-style fixes help.

Fresh seeds (500+), frozen transfer parameters from discovery, no tuning on these instances.
Fixes tested (all classical CPU simulation of p=1 correlations):
  rq_local   : RQAOA + 1-flip descent
  rq_multi   : 8 RQAOA paths, each step picks randomly among the 3 strongest correlations
               (QIRO-style diversification / cheap backtracking), best after descent
  qrr        : quantum relax-and-round (Dupont & Sundar 2024 style): full-problem correlation
               matrix as Gram matrix, 64 random-hyperplane roundings + descent
  rq_tabu    : tabu warm-started from rq_local, extra budget B
  tabu       : cold tabu with the SAME total wall time as rq_tabu (3 reps)
Also: RSG (Hancl et al. code), long tabu (best known).
"""
import sys,time,json,statistics
from pathlib import Path
import numpy as np
from sweep import OUT,make_instance,score,paint_cost,dumpcsv
from benchmark import correlations,exhaustive,improve
from transfer import solve,reduce
from paint_baselines import author_solve,AUTHOR_RSG

def tabu(K,seed,budget,z0=None):
    rng=np.random.default_rng(seed);n=len(K);z=rng.choice([-1,1],n) if z0 is None else z0.copy()
    current=score(K,z);best=current;bestz=z.copy();gain=-2*z*(K@z);until=np.zeros(n,int);step=0;since=0;st=time.perf_counter()
    while time.perf_counter()-st<budget or step==0:
        allowed=(until<=step)|(current+gain>best+1e-9)
        k=int(np.argmax(np.where(allowed,gain,-np.inf)));old=z[k];g=float(gain[k]);current+=g;z[k]*=-1
        gain+=4*z*K[:,k]*old;gain[k]=-g;until[k]=step+5+int(rng.integers(5));step+=1;since+=1
        if current>best+1e-8:best=current;bestz=z.copy();since=0
        if step%128==0:current=score(K,z);gain=-2*z*(K@z)
        if since>10*n:z=rng.choice([-1,1],n);current=score(K,z);gain=-2*z*(K@z);until[:]=0;since=0
    return bestz

def rq_path(K,param,rng,topk):
    active=list(range(len(K)));orig=len(K);trace=[];K=K.copy()
    while len(K)>4 and np.abs(K).max()>1e-12:
        w=-2*K;rms=np.sqrt(np.mean(np.sum(w*w,axis=1)));gamma=param['alpha']/max(rms,1e-15)
        M=correlations(w,gamma,param['beta']);i,j=np.triu_indices(len(K),1);sco=np.abs(M[i,j])
        order=np.argsort(-sco)[:topk];k=int(order[rng.integers(len(order))]) if topk>1 else int(order[0])
        u,v=int(i[k]),int(j[k]);s=1 if M[u,v]>=0 else -1
        trace.append((active[u],active[v],s,float(sco[k])));K=reduce(K,u,v,s);active.pop(v)
    z,_=exhaustive(K);sol=dict(zip(active,z))
    for keep,rem,s,_ in reversed(trace):sol[rem]=s*sol[keep]
    return np.array([sol[i] for i in range(orig)]),trace

def qrr(K,param,rng,samples=64):
    w=-2*K;rms=np.sqrt(np.mean(np.sum(w*w,axis=1)));M=correlations(w,param['alpha']/rms,param['beta'])
    ev,V=np.linalg.eigh((M+M.T)/2);ev=np.clip(ev,0,None);X=V*np.sqrt(ev)
    best=None;bv=-np.inf
    for _ in range(samples):
        z=np.sign(X@rng.normal(size=X.shape[1]));z[z==0]=1;z,_=improve(K,z);v=score(K,z)
        if v>bv:bv=v;best=z
    return best

def run(sizes):
    param=json.loads((OUT/'transfer_model.json').read_text(encoding='utf-8'))['parameters']['paintshop']
    rows=[];errs=[]
    for n,seeds in sizes:
        for seed in seeds:
            K,meta=make_instance('paintshop',n,seed);res={}
            def add(m,z,t,rep=-1):
                c=paint_cost(meta,z);rows.append(dict(n=n,seed=seed,method=m,rep=rep,paint_changes=c,seconds=round(t,4)));return c
            st=time.perf_counter();z,tr=rq_path(K,param,None,1);trq=time.perf_counter()-st;add('rq',z,trq)
            # error diagnosis: how often a greedy elimination is wrong w.r.t. best-known solution
            zl,_=improve(K,z);add('rq_local',zl,trq)
            rng=np.random.default_rng(seed);st=time.perf_counter();best=None;bc=10**9
            for r in range(8):
                zz,_=rq_path(K,param,rng,3);zz,_=improve(K,zz);c=paint_cost(meta,zz)
                if c<bc:bc=c;best=zz
            add('rq_multi8',best,time.perf_counter()-st)
            st=time.perf_counter();zq=qrr(K,param,np.random.default_rng(seed+7));add('qrr64',zq,time.perf_counter()-st)
            st=time.perf_counter();zr=author_solve(meta,AUTHOR_RSG);add('RSG',zr,time.perf_counter()-st)
            for B in (0.25*trq,1.0*trq):
                tag='x0.25' if B<trq else 'x1'
                for rep in range(3):
                    st=time.perf_counter();zt=tabu(K,seed*10+rep,B,zl);add(f'rq_tabu_{tag}',zt,trq+time.perf_counter()-st,rep)
                    st=time.perf_counter();zc=tabu(K,seed*10+rep+5,trq+B);add(f'tabu_{tag}',zc,time.perf_counter()-st,rep)
            st=time.perf_counter();zb=min((tabu(K,999+r,max(5.0,20*trq)) for r in range(2)),key=lambda q:paint_cost(meta,q));bk=add('best_known',zb,time.perf_counter()-st)
            # diagnosis: after imposing the first f of RQAOA's eliminations, what is the best cost still reachable?
            # (long tabu on the reduced problem, expanded back). Loss that appears early = early wrong locks.
            L=len(tr);curve={}
            for f in (0.25,0.5,0.75,1.0):
                k=int(round(f*L));Kr=K.copy();active=list(range(n))
                for keep,rem,s,_ in tr[:k]:
                    u=active.index(keep);v=active.index(rem);Kr=reduce(Kr,u,v,s);active.pop(v)
                zr_=tabu(Kr,seed+int(f*100),max(1.0,4*trq));sol=dict(zip(active,zr_))
                for keep,rem,s,_ in reversed(tr[:k]):sol[rem]=s*sol[keep]
                curve[str(f)]=paint_cost(meta,np.array([sol[i] for i in range(n)]))
            corr=[t[3] for t in tr]
            errs.append(dict(n=n,seed=seed,steps=L,rq=paint_cost(meta,z),best_known=bk,reachable_after_fraction=curve,corr_first_quarter=float(np.mean(corr[:L//4])),corr_last_quarter=float(np.mean(corr[3*L//4:]))))
            print(n,seed,'rq',paint_cost(meta,z),'best',bk,'reachable',curve,'t',round(trq,2),flush=True)
            dumpcsv(OUT/'fixes.csv',rows);(OUT/'fixes_errors.json').write_text(json.dumps(errs,indent=1),encoding='utf-8')

if __name__=='__main__':
    run([(64,range(500,510)),(128,range(600,606))])

"""11/16 probe 2: "both quality and speed" variants on the same 16 fresh paint-shop instances as fixes.py.

backtrack7 : greedy RQAOA path, then 7 extra paths that resume from the greedy path's 7 LEAST-confident
             steps with the 2nd-best pair at that step, continuing greedily. Shared prefix, so extra paths
             only recompute the tail. Best after 1-flip descent.
handoff_f  : RQAOA greedy for the first fraction f of eliminations (f=0.25, 0.5; fixed in advance, not tuned),
             then tabu on the reduced problem for budget B; cold tabu gets the same total wall time.
Budgets: B = 0.25x and 1x the greedy RQAOA time. Frozen transfer parameters. No tuning on these seeds.
"""
import time,json
import numpy as np
from sweep import OUT,make_instance,score,paint_cost,dumpcsv
from benchmark import correlations,exhaustive,improve
from transfer import reduce
from fixes import tabu

def step_choice(K,param,rank):
    w=-2*K;rms=np.sqrt(np.mean(np.sum(w*w,axis=1)));M=correlations(w,param['alpha']/max(rms,1e-15),param['beta'])
    i,j=np.triu_indices(len(K),1);sco=np.abs(M[i,j]);order=np.argsort(-sco);k=int(order[rank])
    u,v=int(i[k]),int(j[k]);return u,v,(1 if M[u,v]>=0 else -1),float(sco[order[0]]),float(sco[order[0]]-sco[order[1]])

def finish(K,active,trace,n):
    z,_=exhaustive(K);sol=dict(zip(active,z))
    for keep,rem,s in reversed(trace):sol[rem]=s*sol[keep]
    return np.array([sol[i] for i in range(n)])

def greedy_from(K,active,trace,param,stop_len=4,record=False):
    K=K.copy();active=list(active);trace=list(trace);states=[]
    while len(K)>stop_len and np.abs(K).max()>1e-12:
        if record:states.append((K.copy(),list(active),list(trace)))
        u,v,s,top,margin=step_choice(K,param,0)
        if record:states[-1]=states[-1]+(margin,)
        trace.append((active[u],active[v],s));K=reduce(K,u,v,s);active.pop(v)
    return K,active,trace,states

def backtrack(K0,param,extra=7):
    n=len(K0);K,a,tr,states=greedy_from(K0,range(n),[],param,record=True)
    cands=[finish(K,a,tr,n)]
    weak=sorted(range(len(states)),key=lambda t:states[t][3])[:extra]   # smallest gap between 1st and 2nd pair
    for t in weak:
        Ks,as_,trs,_=states[t];u,v,s,_,_=step_choice(Ks,param,1)
        Ks2=reduce(Ks,u,v,s);as2=list(as_);trs2=trs+[(as2[u],as2[v],s)];as2.pop(v)
        Kf,af,trf,_=greedy_from(Ks2,as2,trs2,param);cands.append(finish(Kf,af,trf,n))
    cands=[improve(K0,z)[0] for z in cands]
    return max(cands,key=lambda z:score(K0,z))

def handoff(K0,param,frac,budget,seed):
    n=len(K0);K=K0.copy();active=list(range(n));trace=[];steps=int(frac*(n-4))
    for _ in range(steps):
        u,v,s,_,_=step_choice(K,param,0);trace.append((active[u],active[v],s));K=reduce(K,u,v,s);active.pop(v)
    zr=tabu(K,seed,budget);sol=dict(zip(active,zr))
    for keep,rem,s in reversed(trace):sol[rem]=s*sol[keep]
    return np.array([sol[i] for i in range(n)])

def run(sizes):
    param=json.loads((OUT/'transfer_model.json').read_text(encoding='utf-8'))['parameters']['paintshop'];rows=[]
    best={(int(r['n']),int(r['seed'])):int(r['paint_changes']) for r in __import__('csv').DictReader(open(OUT/'fixes.csv',encoding='utf-8-sig')) if r['method']=='best_known'}
    for n,seeds in sizes:
        for seed in seeds:
            K,meta=make_instance('paintshop',n,seed)
            def add(m,z,t,rep=-1):rows.append(dict(n=n,seed=seed,method=m,rep=rep,paint_changes=paint_cost(meta,z),seconds=round(t,4),best_known=best[(n,seed)]))
            st=time.perf_counter();Kf,a,tr,_=greedy_from(K,range(n),[],param);z=finish(Kf,a,tr,n);trq=time.perf_counter()-st;add('rq',z,trq)
            st=time.perf_counter();zb=backtrack(K,param);tb=time.perf_counter()-st;add('backtrack7',zb,tb)
            for rep in range(3):
                st=time.perf_counter();zc=tabu(K,seed*10+rep+50,tb);add('tabu_same_as_backtrack',zc,time.perf_counter()-st,rep)
            for frac in (0.25,0.5):
                for mult in (0.25,1.0):
                    for rep in range(3):
                        st=time.perf_counter();zh=handoff(K,param,frac,mult*trq,seed*10+rep);th=time.perf_counter()-st;add(f'handoff{frac}_B{mult}',zh,th,rep)
                        st=time.perf_counter();zc=tabu(K,seed*10+rep+70,th);add(f'tabu_same_as_handoff{frac}_B{mult}',zc,time.perf_counter()-st,rep)
            print(n,seed,'rq',paint_cost(meta,z),'bt',paint_cost(meta,zb),'best',best[(n,seed)],'t_rq',round(trq,2),'t_bt',round(tb,2),flush=True)
            dumpcsv(OUT/'adaptive.csv',rows)

if __name__=='__main__':
    run([(64,range(500,510)),(128,range(600,606))])

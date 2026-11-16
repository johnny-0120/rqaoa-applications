"""11/17 probe 3: confidence-triggered handoff ("loss-prediction stopping") vs fixed 25% handoff.

Rule: run greedy RQAOA; hand the reduced problem to tabu at the first step whose strongest |<ZZ>|
falls below c x (mean of the first 5 steps' strongest |<ZZ>|). Also hand off at 50% at the latest.
c is chosen on SEPARATE training seeds (n64: 1000-1009, n128: 1100-1105) from {0.95,0.9,0.85,0.8},
then frozen and evaluated on the same 16 test instances as fixes.py/adaptive.py.
Tabu budget after handoff B = 0.25 x full greedy-RQAOA time (same as adaptive.py's best setting);
cold tabu gets the identical total wall time.
"""
import sys,time,json,csv,statistics as S
import numpy as np
from sweep import OUT,make_instance,paint_cost,dumpcsv
from transfer import reduce
from fixes import tabu
from adaptive import step_choice,greedy_from,finish

def rq_with_handoff(K0,param,rule,budget,seed,maxfrac=0.5):
    n=len(K0);K=K0.copy();active=list(range(n));trace=[];tops=[];limit=int(maxfrac*(n-4))
    while len(K)>4:
        u,v,s,top,_=step_choice(K,param,0)
        if rule(tops,top) or len(trace)>=limit:break
        tops.append(top);trace.append((active[u],active[v],s));K=reduce(K,u,v,s);active.pop(v)
    zr=tabu(K,seed,budget);sol=dict(zip(active,zr))
    for keep,rem,s in reversed(trace):sol[rem]=s*sol[keep]
    return np.array([sol[i] for i in range(n)]),len(trace)/(n-4)

def make_rule(c):
    if c=='fixed25':return None
    return lambda tops,top:len(tops)>=5 and top<c*S.mean(tops[:5])

def evaluate(cases,cvals,param,best,tag,rows):
    for n,seed in cases:
        K,meta=make_instance('paintshop',n,seed)
        st=time.perf_counter();Kf,a,tr,_=greedy_from(K,range(n),[],param);trq=time.perf_counter()-st;B=0.25*trq
        for c in cvals:
            for rep in range(3):
                st=time.perf_counter()
                if c=='fixed25':
                    z,fr=rq_with_handoff(K,param,lambda t,x:False,B,seed*10+rep,maxfrac=0.25)
                else:
                    z,fr=rq_with_handoff(K,param,make_rule(c),B,seed*10+rep)
                th=time.perf_counter()-st
                rows.append(dict(set=tag,n=n,seed=seed,method=f'handoff_{c}',rep=rep,paint_changes=paint_cost(meta,z),fraction=round(fr,3),seconds=round(th,4),best_known=best[(n,seed)]))
                zc=tabu(K,seed*10+rep+90,th)
                rows.append(dict(set=tag,n=n,seed=seed,method=f'tabu_same_as_{c}',rep=rep,paint_changes=paint_cost(meta,zc),fraction='',seconds=round(th,4),best_known=best[(n,seed)]))
        print(tag,n,seed,flush=True);dumpcsv(OUT/'autohandoff.csv',rows)

def run():
    param=json.loads((OUT/'transfer_model.json').read_text(encoding='utf-8'))['parameters']['paintshop']
    cvals=[0.95,0.9,0.85,0.8];rows=[]
    train=[(64,s) for s in range(1000,1010)]+[(128,s) for s in range(1100,1106)]
    best={}
    for n,seed in train:
        K,meta=make_instance('paintshop',n,seed);best[(n,seed)]=min(paint_cost(meta,tabu(K,7+r,8.0)) for r in range(2))
    evaluate(train,cvals,param,best,'train',rows)
    gap=lambda c:S.mean(r['paint_changes']-r['best_known'] for r in rows if r['set']=='train' and r['method']==f'handoff_{c}')
    cstar=min(cvals,key=gap);print('train gaps',{c:round(gap(c),3) for c in cvals},'chosen',cstar,flush=True)
    for r in csv.DictReader(open(OUT/'fixes.csv',encoding='utf-8-sig')):
        if r['method']=='best_known':best[(int(r['n']),int(r['seed']))]=int(r['paint_changes'])
    test=[(64,s) for s in range(500,510)]+[(128,s) for s in range(600,606)]
    evaluate(test,['fixed25',cstar],param,best,'test',rows)
    (OUT/'autohandoff_choice.json').write_text(json.dumps(dict(candidates=cvals,train_gaps={str(c):gap(c) for c in cvals},chosen=cstar),indent=1),encoding='utf-8')

if __name__=='__main__':run()

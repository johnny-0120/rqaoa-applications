"""Declared follow-ups: fixed subsets, no tuning on winning instances."""
import csv,json,time,sys
from pathlib import Path
import numpy as np
from benchmark import OUT,ROOT,cloud,wnnc_dense,contract,score,error,improve,writecsv,exhaustive

def main():
    details=json.loads((OUT/'details.json').read_text(encoding='utf-8'))
    rows=list(csv.DictReader((OUT/'main.csv').open(encoding='utf-8-sig')))
    bykey={(int(r['case']),r['method']):r for r in rows}
    declared=dict(timing='local search uses per-case measured RQAOA CPU time; allow one final restart to finish, record actual time',shots='all n=12 control cases with corruption .1 or .25, 256/2048 decision shots, 3 seeds each; exact training unchanged',refinement='all n=16 seed=0 control cases with corruption .1 or .25; 257 grid, gamma range [0,4pi]',wnnc='all geometry cases, published sparse-point l5 widths [.05,.2], 40 iterations; report alongside default l0, not selected per instance',failure_audit='all main-run nonoptimal RQAOA cases; enumerate constrained optimum after every reduction')
    (OUT/'followup_protocol.json').write_text(json.dumps(declared,indent=2),encoding='utf-8')
    budget=[];wnrows=[];shotrows=[];refinements=[];failures=[]
    for d in details:
        case=d['case'];K=np.array(d['K']);truth=np.array(d['truth']);opt=d['exact_info']['optimum'];denom=np.abs(np.triu(K,1)).sum()
        rqrow=bykey[(case,'rqaoa_p1')];limit=float(rqrow['seconds']);rng=np.random.default_rng(123000+case)
        st=time.perf_counter();best=-np.inf;zbest=None;starts=0
        while time.perf_counter()-st<limit or starts==0:
            z,_=improve(K,rng.choice([-1,1],len(K)));v=score(K,z);starts+=1
            if v>best+1e-10:best=v;zbest=z.copy()
        budget.append(dict(case=case,panel=d['panel'],shape=d['shape'],n=d['n'],seed=d['seed'],corruption=d['corruption'],method='local_cpu_matched',score=best,patch_error=error(zbest,truth),is_optimal=int(opt-best<1e-7),seconds=time.perf_counter()-st,target_seconds=limit,starts=starts,rqaoa_score=float(rqrow['score'])))
        if d['panel']=='geometry':
            data=cloud(d['shape'],d['n'],d['seed']);st=time.perf_counter();z,info=wnnc_dense(data,width_min=.05,width_max=.2)
            wnrows.append(dict(case=case,panel=d['panel'],shape=d['shape'],n=d['n'],seed=d['seed'],method='wnnc_dense_numpy_l5',score=score(K,z),patch_error=error(z,truth),seconds=time.perf_counter()-st,**info))
        if d['panel']=='control' and d['n']==12 and d['corruption']>0:
            for shots in [256,2048]:
                for rep in range(3):
                    st=time.perf_counter();z,info=contract(K,shots=shots,seed=43000+case*10+rep)
                    shotrows.append(dict(case=case,shape=d['shape'],n=d['n'],seed=d['seed'],corruption=d['corruption'],shots=shots,rep=rep,score=score(K,z),patch_error=error(z,truth),is_optimal=int(opt-score(K,z)<1e-7),seconds=time.perf_counter()-st,decision_shots=info['decision_shots'],exact_score=float(rqrow['score'])))
        if d['panel']=='control' and d['n']==16 and d['seed']==0 and d['corruption']>0:
            st=time.perf_counter();z,info=contract(K,grid_count=257,gamma_max=4*np.pi)
            refinements.append(dict(case=case,shape=d['shape'],corruption=d['corruption'],original_score=float(rqrow['score']),refined_score=score(K,z),patch_error=error(z,truth),is_optimal=int(opt-score(K,z)<1e-7),seconds=time.perf_counter()-st,evaluations=info['evaluations']))
        if int(rqrow['is_optimal'])==0:
            residual=K.copy();active=list(range(len(K)));offset=0.;previous=opt;trace=[]
            for step,t in enumerate(d['outputs']['rqaoa_p1']['info']['trace'],1):
                u=active.index(t['keep']);v=active.index(t['remove']);s=t['sign'];offset+=s*residual[u,v]
                for k in range(len(residual)):
                    if k not in [u,v]:residual[u,k]=residual[k,u]=residual[u,k]+s*residual[v,k]
                residual=np.delete(np.delete(residual,v,axis=0),v,axis=1);active.pop(v)
                _,info=exhaustive(residual);current=offset+info['optimum']
                assert current<=previous+1e-7
                trace.append(dict(step=step,restricted_optimum=current,loss=max(0.,previous-current),**t));previous=current
            assert abs(previous-float(rqrow['score']))<1e-7
            failures.append(dict(case=case,shape=d['shape'],n=d['n'],seed=d['seed'],corruption=d['corruption'],original_optimum=opt,final_score=float(rqrow['score']),trace=trace))
        if case%12==0:print('followup',case,'/',len(details),flush=True)
    writecsv(OUT/'cpu_matched.csv',budget);writecsv(OUT/'wnnc_sparse_width.csv',wnrows);writecsv(OUT/'shots.csv',shotrows);writecsv(OUT/'parameter_refinement.csv',refinements)
    (OUT/'failure_audit.json').write_text(json.dumps(failures,indent=2),encoding='utf-8')
    print(json.dumps(dict(cpu_matched=len(budget),wnnc_l5=len(wnrows),shot_runs=len(shotrows),refinements=len(refinements),failures=len(failures))),flush=True)

if __name__=='__main__':main()

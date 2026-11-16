"""Predeclared follow-up scope: shot sensitivity, cutoff, scale, failure audit.
Readout-only sensitivity is deliberately not a calibrated device noise model.
"""
import csv,json,time
from pathlib import Path
import numpy as np
import networkx as nx
from rqaoa_experiment import graph,rqaoa,exact,spins,cuts,optimize

OUT=Path('results')

def write(name,rows):
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    shotrows=[]
    # Fixed 3 graphs, 10 independent measurement seeds each. Angles ideal-trained.
    # This isolates errors in recursive pair selection, not optimizer shot noise.
    for family,p,gs in [('bipartite',.8,0),('random',.5,4),('near_bipartite',.15,0)]:
        g,w=graph(family,12,gs,p);opt,_,_=exact(w)
        for shots in (32,256,4096):
            for err in (0.,.03):
                for seed in range(10):
                    r=rqaoa(w,shots=shots,seed=seed,readout_error=err)
                    shotrows.append(dict(family=family,n=12,graph_seed=gs,p=p,shots_per_stage=shots,readout_error=err,measurement_seed=seed,ratio=r['cut']/opt,seconds=r['seconds']))
            write('shot_sensitivity.csv',shotrows)
            print('shots',family,shots,flush=True)
    cutoff=[]
    for family,p in [('random',.5),('random',.8),('near_bipartite',.3)]:
        for seed in range(5):
            g,w=graph(family,16,seed,p);opt,_,_=exact(w)
            for nc in (2,4,8):
                r=rqaoa(w,threshold=nc)
                cutoff.append(dict(family=family,n=16,p=p,seed=seed,cutoff=nc,ratio=r['cut']/opt,seconds=r['seconds'],evaluations=r['evaluations']))
    write('cutoff_sensitivity.csv',cutoff)
    # Bipartite opt is certified by a linear-time coloring; no exponential oracle.
    scale=[]
    for n in (30,50):
        for p in (.2,.8):
            g,w=graph('bipartite',n,0,p);opt=g.number_of_edges()
            st=time.perf_counter();colors=nx.bipartite.color(g);z=np.array([1. if colors[i] else -1. for i in range(n)])
            classical=float(cuts(w,z[None,:])[0]);ct=time.perf_counter()-st
            fit=optimize(w);r=rqaoa(w)
            scale.append(dict(family='bipartite',n=n,p=p,seed=0,edges=opt,qaoa_expected_ratio=fit['expectation']/opt,rqaoa_ratio=r['cut']/opt,rqaoa_seconds=r['seconds'],rqaoa_evaluations=r['evaluations'],coloring_ratio=classical/opt,coloring_seconds=ct))
            write('scale.csv',scale);print('scale',n,p,r['seconds'],flush=True)
    traces=json.loads((OUT/'traces.json').read_text(encoding='utf-8'));fails=[]
    rows=list(csv.DictReader((OUT/'benchmark.csv').open(encoding='utf-8-sig')))
    for row in rows:
        if float(row['rqaoa_ratio'])>=1.-1e-9: continue
        d=traces[int(row['id'])]; n=d['n'];w=np.zeros((n,n))
        for i,j in d['edges']:w[i,j]=w[j,i]=1.
        opt,_,v=exact(w);optimal=spins(n)[v>=opt-1e-9];remaining=len(optimal)
        for stage,t in enumerate(d['rqaoa']['trace'],1):
            keep=optimal[:,t['remove']]==t['sign']*optimal[:,t['keep']]
            optimal=optimal[keep]
            if not len(optimal):
                fails.append(dict(id=int(row['id']),family=d['family'],n=n,seed=d['seed'],p=d['p'],stage=stage,optimal_assignments_before=remaining,first_incompatible_constraint=t,result_ratio=float(row['rqaoa_ratio'])))
                break
            remaining=len(optimal)
    (OUT/'failure_audit.json').write_text(json.dumps(fails,indent=2),encoding='utf-8')
    print('DONE',len(shotrows),'shot trials',len(cutoff),'cutoff trials',len(scale),'scale trials',len(fails),'failures',flush=True)

if __name__=='__main__':main()

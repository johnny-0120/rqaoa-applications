"""Canonical recursive greedy from arXiv:2509.15294v2 appendix B.2.

The earlier first-occurrence pair-insertion baseline is retained under its own
descriptive name. It must not be attributed to the paper's RG algorithm.
"""
import csv,json,time,importlib.util,hashlib
import numpy as np
from sweep import OUT,make_instance,paint_cost,score,dumpcsv

def author_module(name):
    path=OUT.parent/'sources'/name
    spec=importlib.util.spec_from_file_location(name[:-3],path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

AUTHOR_RG=author_module('BPSP_Recursive.py').recursive_greedy
AUTHOR_RSG=author_module('BPSP_Star.py').recursive_star_greedy

def author_solve(meta,fn):
    seq=meta['sequence'];positions={a:tuple(i for i,b in enumerate(seq) if a==b) for a in set(seq)}
    colors,cost=fn(seq,positions)
    for a,(i,j) in positions.items():assert colors[i]+colors[j]==1
    z=np.array([1-2*colors[positions[a][0]] for a in range(len(positions))])
    assert paint_cost(meta,z)==cost
    return z

def recursive_greedy(meta):
    def reconstruct(seq):
        if len(seq)==2:return [1,-1]
        car=seq[-1];m=seq.index(car);sub=seq[:m]+seq[m+1:-1];colors=reconstruct(sub)
        if m==0:first=colors[0]
        elif m==len(colors):first=colors[-1]
        elif colors[m-1]==colors[m]:first=colors[m-1]
        else:first=-colors[-1]
        return colors[:m]+[first]+colors[m:]+[-first]
    seq=meta['sequence'];colors=reconstruct(seq)
    assignment={}
    for a,c in zip(seq,colors):
        if a not in assignment:assignment[a]=c
        else:assert assignment[a]==-c
    return np.array([assignment[a] for a in range(len(assignment))])

def validate():
    from sweep import paint_greedy
    seq=[4,0,0,2,1,1,4,3,2,5,5,3];seen=set();parity=[]
    for a in seq:parity.append(1 if a not in seen else -1);seen.add(a)
    meta=dict(sequence=seq,occurrence_sign=parity)
    assert paint_cost(meta,paint_greedy(meta))==6
    assert paint_cost(meta,recursive_greedy(meta))==5
    assert paint_cost(meta,author_solve(meta,AUTHOR_RG))==5
    assert paint_cost(meta,author_solve(meta,AUTHOR_RSG))==4
    for n in [6,16,64]:
        for seed in range(8):
            K,meta=make_instance('paintshop',n,seed);z=recursive_greedy(meta)
            assert paint_cost(meta,z)==round(meta['paint_cost_constant']-score(K,z))
            assert paint_cost(meta,z)==paint_cost(meta,author_solve(meta,AUTHOR_RG))
            author_solve(meta,AUTHOR_RSG)
    return dict(published_six_car_example_greedy=6,published_six_car_example_RG=5,published_six_car_example_RSG=4,random_feasibility_checks=24,author_RG_agrees=True)

def run():
    checks=validate();rows=[]
    for phase in ['discovery','confirmation']:
        details=json.loads((OUT/(phase+'_details.json')).read_text(encoding='utf-8'))
        for d in details:
            if d['family']!='paintshop':continue
            meta=d['metadata'];K=np.array(d['K'])
            for name,fn in [('paint_RG_canonical',recursive_greedy),('paint_RG_author',lambda m:author_solve(m,AUTHOR_RG)),('paint_RSG_author',lambda m:author_solve(m,AUTHOR_RSG))]:
                st=time.perf_counter();z=fn(meta);elapsed=time.perf_counter()-st
                rows.append(dict(phase=phase,n=d['n'],seed=d['seed'],method=name,score=score(K,z),paint_changes=paint_cost(meta,z),seconds=elapsed,z=json.dumps(z.tolist()),earlier_baseline_name='first_occurrence_pair_insertion_not_paper_RG'))
    dumpcsv(OUT/'paint_canonical_RG.csv',rows)
    (OUT/'paint_baseline_validation.json').write_text(json.dumps(checks,indent=2),encoding='utf-8');print(checks,flush=True)
    print('Canonical RG rows',len(rows),flush=True)

if __name__=='__main__':run()

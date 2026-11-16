"""Follow-up: put dipole geometric interactions in the SAME Ising objective.
Addresses observed misspecification of the simple centroid/normal graph.
This model is a research construction, not claimed novel or the Dipole paper's
global optimization objective. Never uses ground-truth signs to construct K.
"""
import json,time
import numpy as np
from benchmark import OUT,cloud,dipole_field,exhaustive,mst,spectral,sdp,local,anneal,contract,score,error,writecsv

def dipole_graph(data):
    points=data['points'];normals=data['normals'];n=len(points);K=np.zeros((n,n))
    for j in range(n):
        field=dipole_field(points[j],normals[j],points.reshape(-1,3)).reshape(points.shape)
        K[:,j]=np.sum(field*normals,axis=(1,2))
    K=(K+K.T)/2;np.fill_diagonal(K,0)
    return K/max(np.abs(K).max(),1e-15)

def run():
    base=json.loads((OUT/'details.json').read_text(encoding='utf-8'));base=[d for d in base if d['panel']=='geometry']
    (OUT/'geometry_graph_protocol.json').write_text(json.dumps(dict(reason='Naive local dot-product graph failed on torus/slab; exploratory stronger geometric model',cases='all 36 geometry cases unchanged',kernel='sum over points in patch i of n_i dot dipole-field from patch j, symmetrized; no truth used',scope='exploratory follow-up on development cases, not held-out validation',methods=['exact','mst','spectral','sdp256','local32','anneal8','edge_contract','rqaoa_p1','local_cpu_matched']),indent=2),encoding='utf-8')
    rows=[];details=[]
    for d in base:
        data=cloud(d['shape'],d['n'],d['seed']);st=time.perf_counter();K=dipole_graph(data);graph_time=time.perf_counter()-st;truth=data['truth']
        st=time.perf_counter();ze,ex=exhaustive(K,truth);et=time.perf_counter()-st
        methods=[('exact',lambda:(ze,ex)),('mst',lambda:mst(K)),('spectral',lambda:spectral(K)),('sdp256',lambda:sdp(K,d['seed'])),('local32',lambda:local(K,d['seed'],32)),('anneal8',lambda:anneal(K,d['seed'])),('edge_contract',lambda:contract(K,'edge')),('rqaoa_p1',lambda:contract(K))]
        outputs={}
        for name,fn in methods:
            st=time.perf_counter();z,info=fn();elapsed=et if name=='exact' else time.perf_counter()-st;v=score(K,z)
            assert v<=ex['optimum']+1e-7
            rows.append(dict(case=d['case'],shape=d['shape'],n=d['n'],seed=d['seed'],method=name,patch_error=error(z,truth),score=v,optimality_gap=max(0.,ex['optimum']-v),is_optimal=int(ex['optimum']-v<1e-7),seconds=elapsed,graph_seconds=graph_time,truth_error_of_exact=error(ze,truth),optimum_truth_error_min=ex['optimum_truth_error_min'],optimum_truth_error_max=ex['optimum_truth_error_max']))
            outputs[name]=dict(z=np.asarray(z).astype(int).tolist(),info=info)
        limit=elapsed;rng=np.random.default_rng(990000+d['case']);best=-np.inf;bestz=None;starts=0;st=time.perf_counter()
        while time.perf_counter()-st<limit or starts==0:
            z,_=local(K,int(rng.integers(2**30)),1);v=score(K,z);starts+=1
            if v>best+1e-10:best=v;bestz=z.copy()
        rows.append(dict(case=d['case'],shape=d['shape'],n=d['n'],seed=d['seed'],method='local_cpu_matched',patch_error=error(bestz,truth),score=best,optimality_gap=max(0.,ex['optimum']-best),is_optimal=int(ex['optimum']-best<1e-7),seconds=time.perf_counter()-st,graph_seconds=graph_time,truth_error_of_exact=error(ze,truth),optimum_truth_error_min=ex['optimum_truth_error_min'],optimum_truth_error_max=ex['optimum_truth_error_max']))
        details.append(dict(case=d['case'],shape=d['shape'],n=d['n'],seed=d['seed'],K=K.tolist(),truth=truth.tolist(),exact_info=ex,outputs=outputs))
        writecsv(OUT/'geometry_graph.csv',rows);(OUT/'geometry_graph_details.json').write_text(json.dumps(details),encoding='utf-8')
        rq=next(r for r in rows if r['case']==d['case'] and r['method']=='rqaoa_p1');print('geometry_graph',d['case'],d['shape'],d['n'],'RQerr',rq['patch_error'],'exacterr',error(ze,truth),flush=True)

if __name__=='__main__':run()

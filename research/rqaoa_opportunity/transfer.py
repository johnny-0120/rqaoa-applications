"""Discovery-only median parameter transfer, confirmation-only evaluation.

Existing idea (optimization-free/parameter-transfer RQAOA), not a novelty claim.
Fit scale-free gamma * RMS weighted degree and beta on discovery residual traces.
No test-instance parameter optimization, model search or threshold tuning.
"""
import json,time
import numpy as np
from sweep import OUT,score,dumpcsv,timed_local
from benchmark import correlations,exhaustive,improve

def reduce(K,u,v,s):
    K=K.copy()
    for k in range(len(K)):
        if k not in (u,v):K[u,k]=K[k,u]=K[u,k]+s*K[v,k]
    return np.delete(np.delete(K,v,0),v,1)

def learn():
    ds=json.loads((OUT/'discovery_details.json').read_text(encoding='utf-8'));values={}
    for d in ds:
        K=np.array(d['K']);active=list(range(len(K)))
        for t in d['outputs']['rqaoa_p1']['info']['trace']:
            w=-2*K;scale=max(np.abs(w).max(),1e-15);rms=np.sqrt(np.mean(np.sum(w*w,axis=1)))
            values.setdefault(d['family'],[]).append([t['gamma']/scale*rms,t['beta']])
            u=active.index(t['keep']);v=active.index(t['remove']);K=reduce(K,u,v,t['sign']);active.pop(v)
    params={family:dict(alpha=float(np.median(np.array(a)[:,0])),beta=float(np.median(np.array(a)[:,1])),training_residuals=len(a)) for family,a in values.items()}
    (OUT/'transfer_model.json').write_text(json.dumps(dict(rule='per-family medians from discovery only; gamma=alpha/RMS weighted degree',training='n16/20 seeds0..3',parameters=params),indent=2),encoding='utf-8')
    return params

def solve(K,param):
    active=list(range(len(K)));orig=K.copy();trace=[]
    while len(K)>4 and np.abs(K).max()>1e-12:
        w=-2*K;rms=np.sqrt(np.mean(np.sum(w*w,axis=1)));gamma=param['alpha']/max(rms,1e-15)
        M=correlations(w,gamma,param['beta']);i,j=np.triu_indices(len(K),1);sco=np.abs(M[i,j]);k=int(np.flatnonzero(sco>=sco.max()-1e-9)[0]);u,v=int(i[k]),int(j[k]);s=1 if M[u,v]>=0 else -1
        trace.append(dict(keep=active[u],remove=active[v],sign=s,gamma=gamma,beta=param['beta'],remaining=len(K)));K=reduce(K,u,v,s);active.pop(v)
    z,_=exhaustive(K);sol=dict(zip(active,z))
    for t in reversed(trace):sol[t['remove']]=t['sign']*sol[t['keep']]
    z=np.array([sol[i] for i in range(len(orig))]);return z,trace

def run():
    params=learn();ds=json.loads((OUT/'confirmation_details.json').read_text(encoding='utf-8'));rows=[];details=[]
    (OUT/'transfer_protocol.json').write_text(json.dumps(dict(selection='All confirmation cases, no exclusions',parameter_training='Discovery only, per-family median with RMS scaling',test_optimization='none',classical_comparison='same CPU wall budget local search; also local32, SDP1024 and full-RQ from main experiment',postprocessing='single local improvement also reported separately, and SDP/local baselines retained'),indent=2),encoding='utf-8')
    for d in ds:
        K=np.array(d['K']);st=time.perf_counter();z,trace=solve(K,params[d['family']]);elapsed=time.perf_counter()-st
        st=time.perf_counter();zp,_=improve(K,z);polish_seconds=time.perf_counter()-st
        st=time.perf_counter();zc,inf=timed_local(K,990000+d['seed'],elapsed);cseconds=time.perf_counter()-st
        for name,zz,seconds in [('rqaoa_transfer',z,elapsed),('rqaoa_transfer_local',zp,elapsed+polish_seconds),('local_transfer_cpu_matched',zc,cseconds)]:
            v=score(K,zz)
            rows.append(dict(family=d['family'],n=d['n'],seed=d['seed'],method=name,score=v,seconds=seconds,normalized_score=v/d['weight_sum'],paint_changes=d['metadata']['paint_cost_constant']-v if d['family']=='paintshop' else None,exact_optimum=d['exact_optimum'],original_reference=d['best_observed']))
        details.append(dict(family=d['family'],n=d['n'],seed=d['seed'],transfer_z=z.astype(int).tolist(),polished_z=zp.astype(int).tolist(),matched_local_z=zc.astype(int).tolist(),trace=trace,local_budget=inf))
    dumpcsv(OUT/'transfer.csv',rows);(OUT/'transfer_details.json').write_text(json.dumps(details),encoding='utf-8')
    print(json.dumps(params),flush=True);print('Transfer rows',len(rows),flush=True)

if __name__=='__main__':run()

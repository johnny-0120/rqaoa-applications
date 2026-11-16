"""Bounded 128-variable literature-inspired dense-Ising probe.
First near-zero QAOA peak, edge-only reduction, cutoff8, own validated formulas.
Not an exact reproduction or certified parameter optimum of the cited paper.
"""
import json,time
import numpy as np
from scipy.optimize import minimize_scalar
from sweep import OUT,score,sdp_polish,timed_local,dumpcsv
from benchmark import components,correlations,exhaustive,local
from fast_confirmation import timed_tabu
from transfer import reduce

def first_peak(K):
    w=-2*K;scale=max(np.abs(w).max(),1e-15);w=w/scale
    rms=max(np.sqrt(np.mean(np.sum(w*w,axis=1))),1e-15);calls=0
    def at(alpha):
        nonlocal calls
        calls+=1;gamma=alpha/rms;i,j,a,b=components(w,gamma);weights=w[i,j];A=-.5*np.dot(weights,a);B=.25*np.dot(weights,b)
        beta=(np.arctan2(A,B)%(2*np.pi))/4;v=float(.5*weights.sum()-B+np.hypot(A,B));return v,beta
    xs=[0.,.1];ys=[at(0.)[0],at(.1)[0]]
    while xs[-1]<4:
        xs.append(xs[-1]+.1);ys.append(at(xs[-1])[0])
        if ys[-1]<ys[-2]-1e-12:break
    k=int(np.argmax(ys));a,b=xs[max(0,k-1)],xs[min(len(xs)-1,k+1)]
    if b>a:
        r=minimize_scalar(lambda alpha:-at(alpha)[0],bounds=(a,b),method='bounded',options={'xatol':1e-6});alpha=float(r.x)
    else:alpha=xs[k]
    _,beta=at(alpha);gamma=alpha/rms
    return correlations(w,gamma,beta),dict(alpha=alpha,gamma=gamma,beta=beta,evaluations=calls)

def rq(K):
    active=list(range(len(K)));trace=[];calls=0
    while len(K)>8:
        M,fit=first_peak(K);calls+=fit['evaluations'];i,j=np.triu_indices(len(K),1);valid=np.abs(K[i,j])>1e-12;i,j=i[valid],j[valid]
        if not len(i):break
        scores=np.abs(M[i,j]);idx=int(np.flatnonzero(scores>=scores.max()-1e-9)[0]);u,v=int(i[idx]),int(j[idx]);s=1 if M[u,v]>=0 else -1
        trace.append(dict(remaining=len(K),keep=active[u],remove=active[v],sign=s,**fit));K=reduce(K,u,v,s);active.pop(v)
        if len(K)%24==0:print(' dense residual',len(K),flush=True)
    z,_=exhaustive(K);solution=dict(zip(active,z))
    for t in reversed(trace):solution[t['remove']]=t['sign']*solution[t['keep']]
    return np.array([solution[i] for i in range(len(solution))]),dict(trace=trace,evaluations=calls)

def run():
    configs=[(.3,400),(.3,401),(.7,400),(.7,401)]
    (OUT/'dense128_protocol.json').write_text(json.dumps(dict(n=128,configs=configs,weights='rounded Gaussian mean50 variance30',rqaoa='near-zero first peak, edge-only, cutoff8; own implementation inspired by arXiv2501.16419v4',comparators=['SDP1024','SDP1024 + local','local32','CPU-matched local','CPU-matched tabu'],scope='four-instance feasibility probe, no global optimum or significance claim',advantage_metric='raw score and weighted MaxCut, avoid signed energy ratio'),indent=2),encoding='utf-8')
    rows=[];details=[]
    for p,seed in configs:
        rng=np.random.default_rng(seed);mask=np.triu(rng.random((128,128))<p,1);weights=np.round(rng.normal(50,np.sqrt(30),(128,128)));K=np.where(mask,-weights/2,0.);K+=K.T;constant=-np.triu(K,1).sum();outputs={}
        def add(name,z,seconds,info):
            v=score(K,z);rows.append(dict(n=128,density=p,seed=seed,method=name,score=v,maxcut=v+constant,seconds=seconds));outputs[name]=dict(z=z.astype(int).tolist(),info=info)
        st=time.perf_counter();z,info=rq(K);budget=time.perf_counter()-st;add('rqaoa_first_peak',z,budget,info)
        st=time.perf_counter();raw,polished,info=sdp_polish(K,seed);elapsed=time.perf_counter()-st;add('sdp1024',raw,info['raw_seconds'],info);add('sdp1024_local',polished,elapsed,info)
        st=time.perf_counter();z,info=local(K,seed,32);add('local32',z,time.perf_counter()-st,info)
        st=time.perf_counter();z,info=timed_local(K,seed+12,budget);add('local_cpu_matched',z,time.perf_counter()-st,info)
        st=time.perf_counter();z,info=timed_tabu(K,seed+55,budget);add('tabu_cpu_matched',z,time.perf_counter()-st,info)
        details.append(dict(n=128,density=p,seed=seed,K=K.tolist(),outputs=outputs));dumpcsv(OUT/'dense128.csv',rows);(OUT/'dense128_details.json').write_text(json.dumps(details),encoding='utf-8')
        print('dense128 done',p,seed,[(r['method'],r['score']) for r in rows if r['density']==p and r['seed']==seed],flush=True)

if __name__=='__main__':run()

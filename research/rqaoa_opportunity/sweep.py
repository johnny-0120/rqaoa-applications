"""Discovery/confirmation of RQAOA opportunities. No advantage assumed.

New outputs isolated here. Classical CPU simulation, p=1, exact correlations.
Small instances have exhaustive references and exact standard-QAOA sample metrics.
Larger ones have only best-known incumbents and numerical SDP relaxations.
"""
from pathlib import Path
import argparse,csv,hashlib,importlib.metadata,json,platform,sys,time
import numpy as np
import networkx as nx

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'research/normal_orientation'))
from benchmark import fit_qaoa,contract,score,local,anneal,sdp,improve,sign
sys.path.insert(0,str(ROOT/'research'))
from q9_depth_noise import state
OUT=Path(__file__).resolve().parent/'results';OUT.mkdir(parents=True,exist_ok=True)
FAMILIES=['paintshop','dense_positive','sk_gaussian','signed_regular3']

def dumpcsv(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def make_instance(family,n,seed):
    rng=np.random.default_rng(seed+100000*FAMILIES.index(family)+n*1000)
    meta={}
    if family=='paintshop':
        seq=rng.permutation(np.repeat(np.arange(n),2));seen=np.zeros(n,int);parity=[]
        for a in seq:parity.append(1 if seen[a]==0 else -1);seen[a]+=1
        parity=np.array(parity);K=np.zeros((n,n));constant=(2*n-1)/2
        for i in range(2*n-1):
            a,b=int(seq[i]),int(seq[i+1]);coupling=parity[i]*parity[i+1]/2
            if a==b:constant-=coupling
            else:K[a,b]+=coupling;K[b,a]+=coupling
        meta=dict(sequence=seq.tolist(),occurrence_sign=parity.tolist(),paint_cost_constant=constant)
    elif family=='dense_positive':
        mask=np.triu(rng.random((n,n))<.7,1)
        weights=np.round(rng.normal(50,np.sqrt(30),(n,n)))
        K=np.where(mask,-weights/2,0.);K+=K.T
        meta=dict(edge_probability=.7,weight_mean=50,weight_variance=30,scaling='K=-weight/2; score plus constant is MaxCut')
    elif family=='sk_gaussian':
        K=np.triu(rng.normal(0,1/np.sqrt(n),(n,n)),1);K+=K.T
        meta=dict(gaussian_mean=0,gaussian_std=1/np.sqrt(n))
    elif family=='signed_regular3':
        G=nx.random_regular_graph(3,n,seed=seed+1000*n);K=np.zeros((n,n))
        for i,j in sorted(G.edges()):K[i,j]=K[j,i]=rng.choice([-1.,1.])
        meta=dict(degree=3,edge_weights='independent equiprobable +-1')
    else:raise ValueError(family)
    return K,meta

def paint_cost(meta,z):
    colors=np.array(meta['occurrence_sign'])*z[np.array(meta['sequence'])]
    return int(np.count_nonzero(colors[1:]!=colors[:-1]))

def paint_greedy(meta):
    seq=meta['sequence'];assignment={};seen=set();previous=1
    for a in seq:
        if a not in seen:assignment[a]=previous;seen.add(a);color=previous
        else:color=-assignment[a]
        previous=color
    return np.array([assignment[i] for i in range(len(assignment))])

def paint_first_occurrence_insertion(meta):
    # Own first-occurrence insertion heuristic; NOT paper's canonical RG.
    seq=meta['sequence'];order=list(dict.fromkeys(seq));assigned={}
    for car in order:
        selected=set(assigned)|{car};sub=[a for a in seq if a in selected]
        def count(s):
            vals=[];seen=set()
            for a in sub:
                base=s if a==car else assigned[a]
                vals.append(base if a not in seen else -base);seen.add(a)
            return sum(a!=b for a,b in zip(vals,vals[1:]))
        assigned[car]=1 if count(1)<=count(-1) else -1
    return np.array([assigned[i] for i in range(len(assigned))])

def all_scores(K):
    n=len(K);values=np.empty(1<<n)
    for start in range(0,len(values),32768):
        nums=np.arange(start,min(start+32768,len(values)),dtype=np.uint64)
        z=1-2*((nums[:,None]>>np.arange(n,dtype=np.uint64))&1).astype(float)
        values[start:start+len(nums)]=np.sum((z@K)*z,axis=1)/2
    return values

def best_expectation(values,prob,shots):
    order=np.argsort(values,kind='stable');v=values[order];cdf=np.cumsum(prob[order]);cdf=np.clip(cdf,0,1);cdf[-1]=1
    return float(np.dot(v,np.diff(np.r_[0.,cdf**shots])))

def qaoa_distribution(K,fit,values):
    # cut(w) = sum(w)/4 + score(K)/scale for w=-2K/scale.
    cost=fit['normalized_w'].sum()/4+values/fit['scale']
    prob=np.abs(state(cost,np.array([fit['gamma'],fit['beta']])))**2
    return prob/prob.sum()

def sdp_polish(K,seed):
    # Same 1024 random hyperplanes as the cited dense-Ising comparison, plus
    # a separately labelled classical single-spin polish of ALL roundings.
    import cvxpy as cp
    started=time.perf_counter();n=len(K);X=cp.Variable((n,n),symmetric=True)
    problem=cp.Problem(cp.Maximize(cp.sum(cp.multiply(K,X))/2),[cp.diag(X)==1,X>>0])
    problem.solve(solver='SCS',eps=1e-6,max_iters=30000,verbose=False)
    if X.value is None:raise RuntimeError(problem.status)
    eigen,V=np.linalg.eigh((X.value+X.value.T)/2);factor=V@np.diag(np.sqrt(np.maximum(eigen,0)))
    rng=np.random.default_rng(770000+seed);zs=sign(rng.normal(size=(1024,n))@factor.T)
    vals=np.sum((zs@K)*zs,axis=1)/2;k=int(np.argmax(vals));bestz=zs[k].copy();rawz=bestz.copy();best=float(vals[k])
    raw_seconds=time.perf_counter()-started
    for z in zs:
        improved,_=improve(K,z);v=score(K,improved)
        if v>best+1e-9:best=v;bestz=improved.copy()
    return rawz,bestz,dict(raw_seconds=raw_seconds,relaxation=float(problem.value),status=problem.status,rounds=1024,min_eigenvalue=float(eigen.min()),diag_error=float(np.max(np.abs(np.diag(X.value)-1))))

def timed_local(K,seed,limit):
    rng=np.random.default_rng(seed);start=time.perf_counter();best=-np.inf;bestz=None;starts=0
    while time.perf_counter()-start<limit or starts==0:
        z,_=improve(K,rng.choice([-1,1],len(K)));v=score(K,z);starts+=1
        if v>best+1e-9:best=v;bestz=z.copy()
    return bestz,dict(starts=starts,limit=limit,actual_seconds=time.perf_counter()-start)

def instance(family,n,seed,phase):
    K,meta=make_instance(family,n,seed);denom=float(np.triu(np.abs(K),1).sum());outputs={};rows=[]
    def add(name,z,info,seconds):
        value=score(K,np.asarray(z));outputs[name]=dict(z=np.asarray(z).astype(int).tolist(),score=value,info=info)
        rows.append(dict(phase=phase,family=family,n=n,seed=seed,method=name,score=value,seconds=seconds,normalized_score=value/denom,exact_gap=None,normalized_gap=None,reference_kind=None,evaluations=info.get('evaluations',0),starts=info.get('starts',0),paint_changes=paint_cost(meta,np.asarray(z)) if family=='paintshop' else None))
    values=None;exact=None
    if n<=20:
        st=time.perf_counter();values=all_scores(K);exact=float(values.max());idx=int(np.argmax(values));z=1-2*((idx>>np.arange(n))&1)
        add('exact',z,dict(certified='enumeration',states=len(values)),time.perf_counter()-st)
    st=time.perf_counter();M,fit=fit_qaoa(K);fitseconds=time.perf_counter()-st
    expectation=float(np.sum(K*M)/2)
    if n<=20:
        st=time.perf_counter();prob=qaoa_distribution(K,fit,values);distseconds=time.perf_counter()-st
        optprob=float(prob[values>=exact-1e-7].sum());shots=256*(n-4)
        qmetrics=dict(expectation=expectation,best256=best_expectation(values,prob,256),best_decision_matched=best_expectation(values,prob,shots),matched_shots=shots,one_shot_opt_probability=optprob,prob_success256=float(-np.expm1(256*np.log1p(-min(optprob,1.-1e-15)))),fit_seconds=fitseconds,statevector_seconds=distseconds,fit_evaluations=fit['evaluations'])
        rng=np.random.default_rng(990000+seed);samples=rng.choice(len(prob),32,p=prob);best=-np.inf;bestz=None;st=time.perf_counter()
        for x in samples:
            z=1-2*((int(x)>>np.arange(n))&1);z,_=improve(K,z);v=score(K,z)
            if v>best:best=v;bestz=z.copy()
        add('qaoa_sample32_local',bestz,dict(sample_count=32,fit_evaluations=fit['evaluations']),fitseconds+distseconds+time.perf_counter()-st)
    else:qmetrics=dict(expectation=expectation,fit_seconds=fitseconds,fit_evaluations=fit['evaluations'])
    st=time.perf_counter();z,info=contract(K);rqseconds=time.perf_counter()-st;add('rqaoa_p1',z,info,rqseconds)
    st=time.perf_counter();z,info=contract(K,'edge');add('edge_contract',z,info,time.perf_counter()-st)
    for starts in [1,32]:
        st=time.perf_counter();z,info=local(K,seed,starts);add('local'+str(starts),z,info,time.perf_counter()-st)
    st=time.perf_counter();z,info=anneal(K,seed,starts=16,sweeps=160);add('anneal16',z,info,time.perf_counter()-st)
    st=time.perf_counter();raw,polished,info=sdp_polish(K,seed);elapsed=time.perf_counter()-st
    add('sdp1024',raw,info,info['raw_seconds']);add('sdp1024_local',polished,info,elapsed)
    st=time.perf_counter();z,info=timed_local(K,330000+seed,rqseconds);add('local_cpu_matched',z,info,time.perf_counter()-st)
    if family=='paintshop':
        for name,fn in [('paint_greedy',paint_greedy),('paint_first_occurrence_insertion',paint_first_occurrence_insertion)]:
            st=time.perf_counter();z=fn(meta);add(name,z,{},time.perf_counter()-st)
    reference=exact if exact is not None else max(r['score'] for r in rows)
    for r in rows:
        r['exact_gap']=max(0.,reference-r['score']) if exact is not None else None
        r['normalized_gap']=max(0.,reference-r['score'])/(2*denom)
        r['reference_kind']='exact' if exact is not None else 'best_observed_no_certificate'
    data=dict(phase=phase,family=family,n=n,seed=seed,K=K.tolist(),metadata=meta,weight_sum=denom,exact_optimum=exact,best_observed=reference,standard_qaoa=qmetrics,outputs=outputs)
    print(f'{phase} {family} n={n} seed={seed}: RQ={outputs["rqaoa_p1"]["score"]:.5f}, LS32={outputs["local32"]["score"]:.5f}, SDP={outputs["sdp1024"]["score"]:.5f}, RQtime={rqseconds:.2f}s',flush=True)
    return rows,data

def validate():
    from qiskit.quantum_info import Statevector
    from rqaoa_experiment import circuit
    checks=[]
    for family in FAMILIES:
        K,meta=make_instance(family,6,123);values=all_scores(K);M,fit=fit_qaoa(K)
        p=qaoa_distribution(K,fit,values);q=Statevector.from_instruction(circuit(fit['normalized_w'],fit['gamma'],fit['beta'])).probabilities()
        err=float(np.max(np.abs(p-q)));assert err<1e-10
        assert abs(np.dot(values,p)-np.sum(K*M)/2)<1e-9
        if family=='paintshop':
            for x in range(64):
                z=1-2*((x>>np.arange(6))&1);assert abs(paint_cost(meta,z)-(meta['paint_cost_constant']-score(K,z)))<1e-9
            for fn in [paint_greedy,paint_first_occurrence_insertion]:
                z=fn(meta);assert set(z)<=set([-1,1])
        # Expected best-of-s must stay between expectation and exhaustive best.
        eb=best_expectation(values,p,256);assert np.dot(values,p)-1e-9<=eb<=values.max()+1e-9
        checks.append(dict(family=family,statevector_error=err))
    (OUT/'validation.json').write_text(json.dumps(checks,indent=2),encoding='utf-8');print(checks)

def run(phase):
    if phase=='pilot':configs=[(f,16,0) for f in FAMILIES]
    elif phase=='discovery':configs=[(f,n,s) for f in FAMILIES for n in [16,20] for s in range(4)]
    else:
        discovery=json.loads((OUT/'discovery_details.json').read_text(encoding='utf-8'));rank=[]
        # Locked screening criterion: mean improvement over same decision-shot
        # standard QAOA, divided by 2 sum|K|. Does not select on isolated wins.
        for f in FAMILIES:
            ds=[d for d in discovery if d['family']==f]
            gain=np.mean([(d['outputs']['rqaoa_p1']['score']-d['standard_qaoa']['best_decision_matched'])/(2*d['weight_sum']) for d in ds]);rank.append((float(gain),f))
        rank.sort(reverse=True);selected=[f for _,f in rank[:2]]
        (OUT/'selection.json').write_text(json.dumps(dict(criterion='mean normalized RQ minus standard-QAOA matched decision-shot expected best',ranking=rank,selected=selected,interpretation='selection measures benefit versus standard QAOA, not classical advantage'),indent=2),encoding='utf-8')
        configs=[(f,n,s) for f in selected for n,seeds in [(20,range(10,20)),(32,range(20,26)),(64,range(30,34))] for s in seeds]
    protocol=dict(phase=phase,configs=configs,p=1,cutoff=4,source='self-generated graphs following stated literature families, not author datasets',ideal_correlations=True,quantum_hardware=False,screening='exploratory; pilot seed0 is part of discovery; confirmation uses new seeds',standard_qaoa='exact distribution for n<=20, expected best256 and best256*(n-4); aligns decision sampling only, not training shots',classical='local1/local32/16start160sweepSA/SDP1024 with and without local polish/edge contraction/CPU-time-matched local',exact_reference='enumerate n<=20; larger best observed and numerical SDP only',gamma='129-point [0,2pi] search after scaling max coupling, beta analytic; no global certificate')
    (OUT/(phase+'_protocol.json')).write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    rows=[];details=[];started=time.perf_counter()
    for i,(family,n,seed) in enumerate(configs,1):
        rr,d=instance(family,n,seed,phase);rows+=rr;details.append(d)
        dumpcsv(OUT/(phase+'.csv'),rows);(OUT/(phase+'_details.json')).write_text(json.dumps(details),encoding='utf-8')
        print('progress',i,'/',len(configs),'seconds',round(time.perf_counter()-started,1),flush=True)
    (OUT/(phase+'_receipt.json')).write_text(json.dumps(dict(cases=len(configs),rows=len(rows),seconds=time.perf_counter()-started,python=sys.version,platform=platform.platform(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2),encoding='utf-8')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('task',choices=['validate','pilot','discovery','confirmation']);a=p.parse_args()
    validate() if a.task=='validate' else run(a.task)

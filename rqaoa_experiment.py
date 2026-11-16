"""Reproducible p=1 MaxCut / Recursive QAOA benchmark (CPU, ideal).

Independent implementation of p=1 light-cone expectation identities.
References: arXiv:1910.08980, arXiv:2012.03421, arXiv:2408.13207.
Circuit convention: exp(-i beta sum X) exp(-i gamma C)|+>,
C=sum_{i<j} w_ij (1-Z_i Z_j)/2. Qiskit uses rzz(-gamma*w).
No quantum API, credentials, or QPU used.
"""
import argparse, csv, json, time, platform, importlib.metadata
from pathlib import Path
import numpy as np
import networkx as nx
from scipy.optimize import minimize_scalar
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator

OUT = Path('results')

def components(w, gamma):
    n = len(w)
    i, j = np.triu_indices(n, 1)
    mask = (np.arange(n)[None, :] != i[:, None]) & (np.arange(n)[None, :] != j[:, None])
    wi, wj = w[i], w[j]
    pi = np.prod(np.where(mask, np.cos(gamma*wi), 1.), axis=1)
    pj = np.prod(np.where(mask, np.cos(gamma*wj), 1.), axis=1)
    a = -.5*np.sin(gamma*w[i,j])*(pi+pj)
    minus = np.prod(np.where(mask, np.cos(gamma*(wi-wj)), 1.), axis=1)
    plus = np.prod(np.where(mask, np.cos(gamma*(wi+wj)), 1.), axis=1)
    b = .5*(minus-plus)
    return i, j, a, b

def correlations(w, gamma, beta):
    i,j,a,b = components(w,gamma)
    m = np.eye(len(w))
    m[i,j] = m[j,i] = a*np.sin(4*beta)+b*np.sin(2*beta)**2
    return m

def optimize(w):
    """Global beta optimum at each gamma; deterministic gamma grid + refinement.
    Finite gamma search is a heuristic, not a global-optimum certificate.
    All benchmark weights and reduced weights are integers.
    """
    start=time.perf_counter(); calls=0
    if np.max(np.abs(w)) < 1e-12:
        return dict(gamma=0., beta=0., expectation=0., evaluations=0, seconds=0.)
    def at(gamma):
        nonlocal calls
        calls+=1
        i,j,a,b = components(w,gamma); weights=w[i,j]
        A=-.5*np.dot(weights,a); B=.25*np.dot(weights,b)
        base=.5*weights.sum()-B
        beta=(np.arctan2(A,B) % (2*np.pi))/4
        return base+np.hypot(A,B), beta
    count=max(129, 32*int(np.max(np.abs(w)))+1)
    grid=np.linspace(0,2*np.pi,count)
    values=np.array([at(g)[0] for g in grid])
    peaks=[k for k in range(1,count-1) if values[k]>=values[k-1] and values[k]>=values[k+1]]
    peaks=sorted(peaks,key=lambda k: -values[k])[:4]
    best=int(np.argmax(values)); gamma=grid[best]; value,beta=at(gamma)
    for k in peaks:
        r=minimize_scalar(lambda g:-at(g)[0],bounds=(grid[k-1],grid[k+1]),method='bounded',options={'xatol':1e-9})
        val,bet=at(r.x)
        if val>value+1e-10: gamma,value,beta=float(r.x),float(val),float(bet)
    return dict(gamma=float(gamma), beta=float(beta), expectation=float(value),evaluations=calls,seconds=time.perf_counter()-start)

def circuit(w,gamma,beta):
    c=QuantumCircuit(len(w)); c.h(range(len(w)))
    for i,j in zip(*np.triu_indices(len(w),1)):
        if abs(w[i,j])>1e-12: c.rzz(-gamma*float(w[i,j]),int(i),int(j))
    c.rx(2*beta,range(len(w)))
    return c

def spins(n):
    return 1.-2.*((np.arange(1<<n,dtype=np.uint64)[:,None] >> np.arange(n,dtype=np.uint64))&1).astype(float)

def cuts(w,z):
    return .25*(w.sum()-np.einsum('bi,ij,bj->b',z,w,z,optimize=True))

def exact(w):
    z=spins(len(w)); v=cuts(w,z); k=int(np.argmax(v))
    return float(v[k]), z[k], v

def rqaoa(w,threshold=4,mode='quantum',seed=0,shots=0,readout_error=0.):
    start=time.perf_counter(); orig=w.copy(); active=list(range(len(w))); trace=[]; calls=0
    rng=np.random.default_rng(seed)
    while len(active)>threshold:
        n=len(w)
        if mode=='quantum':
            fit=optimize(w); calls+=fit['evaluations']
            if shots:
                # All ZZ values estimated jointly from the same Z-basis shots.
                prob=Statevector.from_instruction(circuit(w,fit['gamma'],fit['beta'])).probabilities()
                samples=rng.choice(len(prob),size=shots,p=prob)
                zz=1.-2.*((samples[:,None] >> np.arange(n))&1)
                if readout_error:
                    zz *= np.where(rng.random(zz.shape)<readout_error,-1.,1.)
                corr=zz.T@zz/shots
            else: corr=correlations(w,fit['gamma'],fit['beta'])
            scores=np.abs(corr); np.fill_diagonal(scores,-1)
            # Standard all-pair RQAOA; lexicographic tie-break within 1e-9.
            candidates=np.argwhere(np.triu(scores>=scores.max()-1e-9,1))
            i,j=map(int,candidates[0]); sign=1 if corr[i,j]>=0 else -1
            strength=float(corr[i,j]); gap=float(scores.max()-np.partition(scores[np.triu_indices(n,1)],-2)[-2]) if n>2 else 0.
        elif mode=='random':
            i,j=sorted(rng.choice(n,2,replace=False).tolist()); sign=int(rng.choice([-1,1])); strength=0.; gap=0.;fit={}
        elif mode=='edge':
            scores=np.abs(w); np.fill_diagonal(scores,-1); i,j=map(int,np.argwhere(np.triu(scores>=scores.max()-1e-9,1))[0])
            sign=-1 if w[i,j]>=0 else 1;strength=0.;gap=0.;fit={}
        else: raise ValueError(mode)
        trace.append(dict(remaining=n,keep=active[i],remove=active[j],sign=sign,correlation=strength,gap=gap,nonzero_edges=int(np.count_nonzero(np.triu(w,1))),max_weight=float(np.abs(w).max()),fit=fit))
        # z_j=sign*z_i. Constants do not affect subsequent optimizations.
        # Values are always reported on the original graph after reconstruction.
        w=w.copy()
        for k in range(n):
            if k not in (i,j): w[i,k]=w[k,i]=w[i,k]+sign*w[j,k]
        w=np.delete(np.delete(w,j,axis=0),j,axis=1);active.pop(j)
    _,z,_=exact(w); solution=dict(zip(active,z))
    for t in reversed(trace): solution[t['remove']]=t['sign']*solution[t['keep']]
    z=np.array([solution[i] for i in range(len(orig))])
    return dict(cut=float(cuts(orig,z[None,:])[0]),spins=z.astype(int).tolist(),trace=trace,evaluations=calls,seconds=time.perf_counter()-start)

def local_search(w,seed,starts=32):
    rng=np.random.default_rng(seed); best=-np.inf; start=time.perf_counter()
    for _ in range(starts):
        z=rng.choice([-1.,1.],len(w))
        for __ in range(10000):
            gain=z*(w@z); k=int(np.argmax(gain))
            if gain[k]<1e-10: break
            z[k]*=-1
        best=max(best,cuts(w,z[None,:])[0])
    return float(best),time.perf_counter()-start

def graph(family,n,seed,p):
    rng=np.random.default_rng(seed)
    if family in ('bipartite','near_bipartite'):
        g=nx.bipartite.random_graph(n//2,n-n//2,p if family=='bipartite' else .6,seed=seed)
        if family=='near_bipartite':
            for i in range(n):
                for j in range(i+1,n):
                    if (i<n//2)==(j<n//2) and rng.random()<p:g.add_edge(i,j)
    elif family=='random': g=nx.gnp_random_graph(n,p,seed=seed)
    elif family=='path': g=nx.path_graph(n)
    elif family=='cycle':g=nx.cycle_graph(n)
    elif family=='star':g=nx.star_graph(n-1)
    elif family=='complete':g=nx.complete_graph(n)
    else:raise ValueError(family)
    return g,nx.to_numpy_array(g,nodelist=range(n))

def verify():
    OUT.mkdir(exist_ok=True);rng=np.random.default_rng(1107);errors=[];suberrors=[]
    for n in (2,3,5,7):
        for trial in range(4):
            w=rng.integers(-3,4,size=(n,n)).astype(float);w=np.triu(w,1);w=w+w.T
            ga,be=rng.uniform(-2,2,size=2); z=spins(n)
            p=Statevector.from_instruction(circuit(w,ga,be)).probabilities()
            actual=z.T@(p[:,None]*z)
            errors.append(float(np.max(np.abs(actual-correlations(w,ga,be)))))
            # Explicitly verify full original objective under a substitution,
            # including the constant absent from the reduced cut representation.
            i,j,s=0,n-1,-1; reduced=np.delete(np.delete(w,j,0),j,1)
            for k in range(n-1):
                if k!=i:reduced[i,k]=reduced[k,i]=w[i,k]+s*w[j,k]
            zr=spins(n-1);zo=np.insert(zr,j,s*zr[:,i],axis=1)
            diff=cuts(w,zo)-cuts(reduced,zr)
            suberrors.append(float(np.ptp(diff)))
    g,w=graph('cycle',6,0,0); f=optimize(w)
    c=circuit(w,f['gamma'],f['beta']);c.measure_all()
    counts=AerSimulator().run(c,shots=4096,seed_simulator=1107).result().get_counts()
    val=sum(k*cuts(w,np.array([[1 if b=='0' else -1 for b in key[::-1]]]))[0] for key,k in counts.items())/4096
    result=dict(correlation_cases=len(errors),max_correlation_error=max(errors),max_substitution_error=max(suberrors),aer_cycle6_expected=f['expectation'],aer_cycle6_shots=4096,aer_cycle6_empirical_cut=val)
    assert max(errors)<1e-10 and max(suberrors)<1e-10
    assert abs(val-f['expectation'])<.15
    (OUT/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(result,flush=True)

def run():
    OUT.mkdir(exist_ok=True); rows=[]; detail=[]; started=time.perf_counter()
    configs=[(f,p) for f in ('bipartite','random') for p in (.2,.5,.8)]+[('near_bipartite',p) for p in (.05,.15,.3)]
    cases=[(f,n,s,p) for f,p in configs for n in (12,16) for s in range(5)]
    cases += [(f,n,0,0.) for f in ('path','cycle','star','complete') for n in (12,16)]
    for idx,(family,n,seed,p) in enumerate(cases):
        t=time.perf_counter();g,w=graph(family,n,seed,p);opt,_,values=exact(w)
        f=optimize(w);prob=Statevector.from_instruction(circuit(w,f['gamma'],f['beta'])).probabilities()
        rng=np.random.default_rng(10000+seed)
        batches=rng.choice(len(prob),size=(20,256),p=prob);best=values[batches].max(axis=1)
        rq=rqaoa(w,seed=seed);edge=rqaoa(w,mode='edge',seed=seed);rand=rqaoa(w,mode='random',seed=seed)
        greedy,gt=local_search(w,seed)
        # Same post-training final sampling allowance as 256 per RQAOA stage.
        matched=values[rng.choice(len(prob),size=(20,256*(n-4)),p=prob)].max(axis=1)
        row=dict(id=idx,family=family,n=n,seed=seed,p=p,edges=g.number_of_edges(),density=nx.density(g),triangles=sum(nx.triangles(g).values())//3,optimum=opt,qaoa_expected_ratio=f['expectation']/opt,qaoa_success_probability=float(prob[values>=opt-1e-9].sum()),qaoa_best256_mean_ratio=float(best.mean()/opt),qaoa_best256_opt_frequency=float(np.mean(best>=opt-1e-9)),qaoa_matched_finalshots_mean_ratio=float(matched.mean()/opt),rqaoa_ratio=rq['cut']/opt,edge_contract_ratio=edge['cut']/opt,random_contract_ratio=rand['cut']/opt,local32_ratio=greedy/opt,qaoa_evaluations=f['evaluations'],rqaoa_evaluations=rq['evaluations'],qaoa_fit_seconds=f['seconds'],rqaoa_seconds=rq['seconds'],local32_seconds=gt,case_seconds=time.perf_counter()-t)
        rows.append(row);detail.append(dict(id=idx,family=family,n=n,seed=seed,p=p,edges=list(g.edges()),qaoa_fit=f,rqaoa=rq))
        with (OUT/'benchmark.csv').open('w',newline='',encoding='utf-8-sig') as h:
            writer=csv.DictWriter(h,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        (OUT/'traces.json').write_text(json.dumps(detail,indent=2),encoding='utf-8')
        print(f'{idx+1}/{len(cases)} {family} n={n} p={p} seed={seed}: E={row["qaoa_expected_ratio"]:.3f} best256={row["qaoa_best256_mean_ratio"]:.3f} RQ={row["rqaoa_ratio"]:.3f} LS={row["local32_ratio"]:.3f} {row["case_seconds"]:.1f}s',flush=True)
    env=dict(python=platform.python_version(),platform=platform.platform(),packages={p:importlib.metadata.version(p) for p in ('numpy','scipy','networkx','qiskit','qiskit-aer','matplotlib')},cases=len(rows),total_seconds=time.perf_counter()-started,parameters=dict(p=1,n_cutoff=4,graph_seeds=list(range(5)),sample_batches=20,shots_per_batch=256,local_starts=32,gamma_grid_min=129,gamma_peaks_refined=4,tie_tolerance=1e-9))
    (OUT/'environment.json').write_text(json.dumps(env,indent=2),encoding='utf-8')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('task',choices=['verify','run']);args=parser.parse_args()
    verify() if args.task=='verify' else run()

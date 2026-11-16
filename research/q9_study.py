"""Question 9 extension: graph structure, paired regular rewiring, p=1 RQAOA.
New results are isolated under research/q9_results/; old pilot is untouched.
"""
import csv,json,sys,time,hashlib,math
from pathlib import Path
import numpy as np
import networkx as nx
from scipy.optimize import minimize_scalar
from qiskit.quantum_info import Statevector

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from rqaoa_experiment import components,correlations,optimize,circuit,spins,cuts,local_search
OUT=ROOT/'research/q9_results';OUT.mkdir(exist_ok=True,parents=True)

def dumpcsv(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as h:
        writer=csv.DictWriter(h,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def restricted_fit(w):
    """Optimize within Bae et al. gamma/beta region, with finite grid refinement.
    Best beta on its interval is exact for each gamma. No global certificate.
    """
    calls=0;st=time.perf_counter(); bound=np.pi/(2*np.abs(w).max());epsilon=1e-8
    def evaluate(gamma):
        nonlocal calls
        calls+=1;i,j,a,b=components(w,gamma);weights=w[i,j]
        A=-.5*np.dot(weights,a);B=.25*np.dot(weights,b);base=.5*weights.sum()-B
        angle=np.arctan2(A,B)%(2*np.pi)
        candidates=[4*epsilon,np.pi-4*epsilon]
        if 4*epsilon<angle<np.pi-4*epsilon:candidates.append(angle)
        values=[base+A*np.sin(t)+B*np.cos(t) for t in candidates]
        k=int(np.argmax(values));return float(values[k]),float(candidates[k]/4)
    grid=np.linspace(-bound,bound,129);vals=np.array([evaluate(g)[0] for g in grid])
    k=int(np.argmax(vals));ga=float(grid[k]);value,be=evaluate(ga)
    peaks=sorted([i for i in range(1,128) if vals[i]>=vals[i-1] and vals[i]>=vals[i+1]],key=lambda i:-vals[i])[:4]
    for k in peaks:
        r=minimize_scalar(lambda g:-evaluate(g)[0],bounds=(grid[k-1],grid[k+1]),method='bounded',options={'xatol':1e-10})
        v,beta=evaluate(r.x)
        if v>value+1e-10:ga,value,be=float(r.x),v,beta
    return dict(gamma=ga,beta=be,expectation=value,evaluations=calls,seconds=time.perf_counter()-st)

def recursive(w,variant,initial_fit=None):
    start=time.perf_counter();original=w.copy();active=list(range(len(w)));trace=[];evaluations=0
    while len(active)>4 and np.max(np.abs(w))>1e-12:
        fit=initial_fit if not trace and initial_fit is not None and variant!='restricted' else (restricted_fit(w) if variant=='restricted' else optimize(w))
        evaluations+=fit['evaluations'];corr=correlations(w,fit['gamma'],fit['beta'])
        i,j=np.triu_indices(len(w),1)
        if variant!='all':valid=np.abs(w[i,j])>1e-12;i,j=i[valid],j[valid]
        scores=np.abs(corr[i,j]);k=int(np.flatnonzero(scores>=scores.max()-1e-9)[0]);u,v=int(i[k]),int(j[k]);s=1 if corr[u,v]>=0 else -1
        trace.append(dict(remaining=len(w),keep=active[u],remove=active[v],sign=s,correlation=float(corr[u,v]),candidate_gap=float(scores.max()-np.partition(scores,-2)[-2]) if len(scores)>1 else 0.,edge_count=int(np.count_nonzero(np.triu(w,1))),max_weight=float(np.abs(w).max()),fit=fit))
        w=w.copy()
        for a in range(len(w)):
            if a not in (u,v):w[u,a]=w[a,u]=w[u,a]+s*w[v,a]
        w=np.delete(np.delete(w,v,0),v,1);active.pop(v)
    zs=spins(len(w));values=cuts(w,zs);z=zs[int(np.argmax(values))];result=dict(zip(active,z))
    for t in reversed(trace):result[t['remove']]=t['sign']*result[t['keep']]
    z=np.array([result[i] for i in range(len(original))]);value=float(cuts(original,z[None,:])[0])
    return dict(cut=value,spins=z.astype(int).tolist(),trace=trace,evaluations=evaluations,seconds=time.perf_counter()-start)

def gw(w,seed):
    import cvxpy as cp
    start=time.perf_counter();n=len(w);X=cp.Variable((n,n),symmetric=True)
    problem=cp.Problem(cp.Maximize(.25*(w.sum()-cp.sum(cp.multiply(w,X)))),[cp.diag(X)==1,X>>0])
    problem.solve(solver='SCS',eps=1e-6,max_iters=20000,verbose=False)
    assert problem.status in ['optimal','optimal_inaccurate'],problem.status
    eigen,U=np.linalg.eigh((X.value+X.value.T)/2);V=U@np.diag(np.sqrt(np.maximum(eigen,0)))
    rng=np.random.default_rng(70000+seed);z=np.where(rng.standard_normal((256,n))@V.T>=0,1.,-1.)
    values=cuts(w,z);return dict(cut=float(values.max()),seconds=time.perf_counter()-start,sdp_value=float(problem.value),status=problem.status,diagonal_error=float(np.abs(np.diag(X.value)-1).max()),minimum_eigenvalue=float(eigen.min()))

def family_graph(family,n,seed):
    if family=='path':g=nx.path_graph(n)
    elif family=='star':g=nx.star_graph(n-1)
    elif family=='cycle':g=nx.cycle_graph(n)
    elif family=='complete':g=nx.complete_graph(n)
    elif family=='tree':g=nx.random_labeled_tree(n,seed=seed)
    elif family=='bipartite':g=nx.bipartite.random_graph(n//2,n-n//2,.5,seed=seed)
    elif family=='random':g=nx.gnp_random_graph(n,.5,seed=seed)
    elif family=='community':g=nx.stochastic_block_model([n//2,n-n//2],[[.7,.15],[.15,.7]],seed=seed)
    else:raise ValueError(family)
    return nx.convert_node_labels_to_integers(g)

def rewired_chain(seed):
    """Connected regular bipartite base; cumulative degree-preserving switches.
    Conditioning on the ability to perform all 12 prescribed switches is explicit.
    """
    for attempt in range(200):
        rng=np.random.default_rng(90000+1000*seed+attempt);g=nx.Graph();g.add_nodes_from(range(16))
        for _ in range(4):
            for trial in range(10000):
                perm=rng.permutation(8)
                if all(not g.has_edge(i,8+int(perm[i])) for i in range(8)):
                    g.add_edges_from((i,8+int(perm[i])) for i in range(8));break
            else:raise RuntimeError('Could not create bipartite regular graph')
        if not nx.is_connected(g):continue
        snapshots={0:g.copy()};success=True
        for step in range(1,13):
            cross=[(min(a,b),max(a,b)) for a,b in g.edges() if (a<8)!=(b<8)]
            candidates=[]
            for x in range(len(cross)):
                for y in range(x+1,len(cross)):
                    a,b=cross[x];c,d=cross[y]
                    if len({a,b,c,d})==4 and not g.has_edge(a,c) and not g.has_edge(b,d):candidates.append((a,b,c,d))
            rng.shuffle(candidates)
            for a,b,c,d in candidates:
                h=g.copy();h.remove_edges_from([(a,b),(c,d)]);h.add_edges_from([(a,c),(b,d)])
                if nx.is_connected(h):g=h;break
            else:success=False;break
            if step in (1,2,4,8,12):snapshots[step]=g.copy()
        if success:return snapshots,attempt
    raise RuntimeError(f'No complete rewiring chain for seed {seed}')

def best_of(prob,values,shots=256):
    unique,inv=np.unique(values,return_inverse=True);mass=np.bincount(inv,weights=prob);cdf=np.clip(np.cumsum(mass),0,1);cdf[-1]=1
    return float(np.dot(unique,np.diff(np.r_[0.,cdf**shots])))

def run():
    protocol=json.loads((ROOT/'research/q9_protocol.json').read_text(encoding='utf-8'));(OUT/'protocol.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding='utf-8')
    cases=[]
    for family in protocol['families']:
        for n in (12,16):
            for seed in (range(10) if family in ('tree','bipartite','community','random') else [0]):
                cases.append(dict(panel='families',family=family,n=n,seed=seed,switches=-1,construction_attempt=0,g=family_graph(family,n,seed)))
    for seed in range(20):
        chain,attempt=rewired_chain(seed)
        for k,g in chain.items():cases.append(dict(panel='rewiring',family='regular4',n=16,seed=seed,switches=k,construction_attempt=attempt,g=g))
    rows=[];details=[];errors=[];start=time.perf_counter()
    for idx,case in enumerate(cases):
        g=case['g'];n=case['n'];seed=case['seed'];w=nx.to_numpy_array(g,nodelist=range(n));z=spins(n);values=cuts(w,z);opt=float(values.max());m=g.number_of_edges()
        fit=optimize(w);prob=Statevector.from_instruction(circuit(w,fit['gamma'],fit['beta'])).probabilities();err=abs(float(prob@values)-fit['expectation']);errors.append(err);assert err<1e-8
        e=recursive(w,'edge',fit);a=recursive(w,'all',fit);r=recursive(w,'restricted');ls,lstime=local_search(w,seed,32);sd=gw(w,seed)
        for res in (e,a,r):
            zz=res['spins'];assert sum(zz[i]!=zz[j] for i,j in g.edges())==res['cut'];assert res['cut']<=opt+1e-9
        rng=np.random.default_rng(310000+idx);batches=values[rng.choice(len(prob),size=(30,256),p=prob)].max(axis=1)
        row={k:v for k,v in case.items() if k!='g'};row['id']=idx
        row.update(edges=m,density=nx.density(g),degree_variance=float(np.var([d for _,d in g.degree()])),clustering=nx.average_clustering(g),triangles=sum(nx.triangles(g).values())//3,components=nx.number_connected_components(g),is_bipartite=int(nx.is_bipartite(g)),core_degeneracy=max(nx.core_number(g).values()),optimum=opt,opt_cut_degeneracy=int(np.count_nonzero(values>=opt-1e-9)//2),frustration_fraction=1-opt/m,qaoa_expected_ratio=fit['expectation']/opt,qaoa_success_probability=float(prob[values>=opt-1e-9].sum()),qaoa_best256_exact_ratio=best_of(prob,values)/opt,qaoa_best256_batch_mean=float(batches.mean()/opt),qaoa_best256_batch_sd=float(batches.std(ddof=1)/opt),edge_rq_ratio=e['cut']/opt,all_rq_ratio=a['cut']/opt,restricted_rq_ratio=r['cut']/opt,local32_ratio=ls/opt,gw256_ratio=sd['cut']/opt,qaoa_fit_seconds=fit['seconds'],edge_rq_seconds=e['seconds'],restricted_rq_seconds=r['seconds'],local32_seconds=lstime,gw_seconds=sd['seconds'],edge_rq_evaluations=e['evaluations'],restricted_rq_evaluations=r['evaluations'],gw_sdp_value=sd['sdp_value'],gw_min_eigenvalue=sd['minimum_eigenvalue'],gw_diagonal_error=sd['diagonal_error'])
        if case['panel']=='rewiring':assert m==32 and all(d==4 for _,d in g.degree()) and nx.is_connected(g)
        rows.append(row);details.append(dict(id=idx,metadata={k:v for k,v in case.items() if k!='g'},edges=sorted([sorted(edge) for edge in g.edges()]),qaoa_fit=fit,edge_rq=e,all_rq=a,restricted_rq=r))
        dumpcsv(OUT/'main.csv',rows)
        (OUT/'graphs_and_traces.json').write_text(json.dumps(details,indent=2),encoding='utf-8')
        if idx%10==0 or idx==len(cases)-1:print(f'{idx+1}/{len(cases)} {case["panel"]} {case["family"]} k={case["switches"]} E={row["qaoa_expected_ratio"]:.3f} best={row["qaoa_best256_exact_ratio"]:.3f} edge={row["edge_rq_ratio"]:.3f} restricted={row["restricted_rq_ratio"]:.3f}',flush=True)
    hashes=[hashlib.sha256(json.dumps(d['edges']).encode()).hexdigest() for d in details]
    # Exact isomorphism check across distinct paired seed bases is reported separately if needed;
    # labeled duplicate detection does not prove non-isomorphism.
    receipt=dict(cases=len(rows),family_cases=sum(r['panel']=='families' for r in rows),rewiring_cases=sum(r['panel']=='rewiring' for r in rows),unique_labeled_edge_lists=len(set(hashes)),seconds=time.perf_counter()-start,max_qiskit_expectation_error=max(errors),all_decoded_scores_checked=True,all_regular_invariants_checked=True,gw_note='Numerical SCS relaxation plus eigenvalue clipping and 256 hyperplane roundings; not rigorous SDP certificates')
    (OUT/'main_audit.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(receipt,flush=True)

if __name__=='__main__':run()

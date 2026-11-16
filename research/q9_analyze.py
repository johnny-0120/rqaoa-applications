"""Aggregate Question9 results, paired bootstrap and independent diagnostics."""
import csv,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'research'))
from rqaoa_experiment import cuts,spins,components
from q9_study import OUT
from scipy.optimize import minimize_scalar

def read(name):return list(csv.DictReader((OUT/name).open(encoding='utf-8-sig')))
def avg(rows,key):return float(np.mean([float(r[key]) for r in rows]))
def bootstrap(values,seed=1107):
    values=np.array(values,float);rng=np.random.default_rng(seed);sample=rng.choice(values,size=(20000,len(values)),replace=True).mean(axis=1)
    return dict(mean=float(values.mean()),lower=float(np.quantile(sample,.025)),upper=float(np.quantile(sample,.975)),unit='graph-seed paired difference',replicates=20000)

def dense_fit(w):
    i,j=np.triu_indices(len(w),1);weights=w[i,j]
    def f(ga):
        _,_,a,b=components(w,ga);A=-.5*np.dot(weights,a);B=.25*np.dot(weights,b)
        return float(.5*weights.sum()-B+np.hypot(A,B))
    grid=np.linspace(0,2*np.pi,4097);v=np.array([f(g) for g in grid]);best=float(v.max())
    peaks=sorted([k for k in range(1,len(v)-1) if v[k]>=v[k-1] and v[k]>=v[k+1]],key=lambda k:-v[k])[:12]
    for k in peaks:
        r=minimize_scalar(lambda ga:-f(ga),bounds=(grid[k-1],grid[k+1]),method='bounded',options={'xatol':1e-12});best=max(best,float(-r.fun))
    return best

def main():
    rows=read('main.csv');details=json.loads((OUT/'graphs_and_traces.json').read_text(encoding='utf-8'))
    keys=['qaoa_expected_ratio','qaoa_best256_exact_ratio','edge_rq_ratio','all_rq_ratio','restricted_rq_ratio','local32_ratio','gw256_ratio']
    fam={f:dict(n=len(a),**{k:avg(a,k) for k in keys},edge_success=sum(float(x['edge_rq_ratio'])>1-1e-9 for x in a)) for f in ['path','star','cycle','tree','complete','bipartite','community','random'] if(a:=[x for x in rows if x['family']==f])}
    controls={str(k):dict(n=len(a),**{t:avg(a,t) for t in keys+['frustration_fraction','triangles','opt_cut_degeneracy']},qaoa_cut_fraction=float(np.mean([float(x['qaoa_expected_ratio'])*float(x['optimum'])/32 for x in a])),edge_success=sum(float(x['edge_rq_ratio'])>1-1e-9 for x in a)) for k in [0,1,2,4,8,12] if(a:=[x for x in rows if x['panel']=='rewiring' and int(x['switches'])==k])}
    zero={int(x['seed']):x for x in rows if x['panel']=='rewiring' and int(x['switches'])==0};last={int(x['seed']):x for x in rows if x['panel']=='rewiring' and int(x['switches'])==12}
    pairs={k:bootstrap([float(last[s][k])-float(zero[s][k]) for s in sorted(zero)]) for k in keys}
    pairs['cut_fraction_change']=bootstrap([(float(last[s]['qaoa_expected_ratio'])*float(last[s]['optimum'])-float(zero[s]['qaoa_expected_ratio'])*32)/32 for s in sorted(zero)])
    for k in ('qaoa_expected_ratio','qaoa_best256_exact_ratio','edge_rq_ratio'):
        matrix=np.array([[float(next(x for x in rows if x['panel']=='rewiring' and int(x['seed'])==seed and int(x['switches'])==step)[k]) for step in [0,1,2,4,8,12]] for seed in range(20)])
        assert matrix.shape==(20,6)
    positive=[];sampling_counts=[];st=time.perf_counter()
    for row,d in zip(rows,details):
        # Independent original-edge enumeration checks for all decoded RQAOA cuts.
        for variant in ('edge_rq','all_rq','restricted_rq'):
            z=d[variant]['spins'];actual=sum(z[i]!=z[j] for i,j in d['edges']);assert actual==d[variant]['cut']
        if float(row['edge_rq_ratio'])<1-1e-9:
            n=int(row['n']);w=np.zeros((n,n))
            for i,j in d['edges']:w[i,j]=w[j,i]=1
            z=spins(n);score=cuts(w,z);alive=np.ones(len(z),bool);total=0.
            for step,t in enumerate(d['edge_rq']['trace'],1):
                before=float(score[alive].max());relation=z[:,t['remove']]==t['sign']*z[:,t['keep']];alive&=relation;after=float(score[alive].max());loss=before-after;total+=loss
                if loss>1e-9:positive.append(dict(id=int(row['id']),panel=row['panel'],family=row['family'],seed=int(row['seed']),switches=int(row['switches']),step=step,loss=loss,before=before,after=after,correlation=t['correlation'],candidate_gap=t['candidate_gap']))
            assert abs(total-(float(row['optimum'])-d['edge_rq']['cut']))<1e-9
        # Probability-of-at-least-one optimum is exact for independent draws.
        success=float(row['qaoa_success_probability']);sampling_counts.append(1-(1-success)**256)
    # Optimization-resolution sensitivity: all failed roots, small independent check.
    resolution=[]
    for row,d in zip(rows,details):
        if float(row['edge_rq_ratio'])>=1-1e-9:continue
        n=int(row['n']);w=np.zeros((n,n))
        for i,j in d['edges']:w[i,j]=w[j,i]=1
        original=d['qaoa_fit']['expectation'];dense=dense_fit(w)
        resolution.append(dict(id=int(row['id']),original_expected=original,dense_expected=dense,gain=dense-original))
    summary=dict(scope='Question9 structure study, no loss predictor trained',families=fam,rewiring=controls,paired_endpoint_changes=pairs,total_success={k:sum(float(x[k])>=1-1e-9 for x in rows) for k in ['edge_rq_ratio','restricted_rq_ratio','local32_ratio','gw256_ratio']},variant_cut_differences=dict(edge_vs_all=sum(abs(float(x['edge_rq_ratio'])-float(x['all_rq_ratio']))>1e-9 for x in rows),edge_vs_restricted=sum(abs(float(x['edge_rq_ratio'])-float(x['restricted_rq_ratio']))>1e-9 for x in rows)),mean_timings={k:avg(rows,k) for k in ['qaoa_fit_seconds','edge_rq_seconds','restricted_rq_seconds','local32_seconds','gw_seconds']},failure_loss_events=positive,optimization_resolution_checks=resolution,max_dense_root_gain=max(x['gain'] for x in resolution),caution='Dense scan only checks original failed-graph roots, not every reduced instance or global optimality certificate.',diagnostic_seconds=time.perf_counter()-st)
    (OUT/'analysis.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(dict(success=summary['total_success'],paired=summary['paired_endpoint_changes'],variants=summary['variant_cut_differences'],dense_gain=summary['max_dense_root_gain'],loss_events=len(positive)),indent=2))

if __name__=='__main__':main()

"""Paired held-out summaries, independent rescoring and audit receipts."""
import csv,json,hashlib,statistics
import numpy as np
from sweep import OUT,paint_cost

def read(name):return list(csv.DictReader((OUT/name).open(encoding='utf-8-sig')))
def mean(x):return float(statistics.mean(x))

def fresh_summary():
    rows=read('fast_confirmation.csv');groups={}
    for r in rows:groups.setdefault((int(r['n']),r['method'],int(r['seed'])),[]).append(r)
    case={key:dict(swaps=mean(float(r['paint_changes']) for r in rr),seconds=float(statistics.median(float(r['seconds']) for r in rr))) for key,rr in groups.items()}
    result=[];comparisons=[]
    for n in [32,64]:
        methods=sorted({key[1] for key in case if key[0]==n});seeds=sorted({key[2] for key in case if key[0]==n})
        for m in methods:
            a=[case[(n,m,s)] for s in seeds]
            result.append(dict(n=n,method=m,cases=len(a),mean_swaps=mean(x['swaps'] for x in a),median_seconds=float(statistics.median(x['seconds'] for x in a))))
            if m=='rqaoa_transfer':continue
            diffs=np.array([case[(n,'rqaoa_transfer',s)]['swaps']-case[(n,m,s)]['swaps'] for s in seeds])
            rng=np.random.default_rng(11212);boot=diffs[rng.integers(0,len(diffs),(5000,len(diffs)))].mean(1)
            comparisons.append(dict(n=n,comparator=m,wins=int(np.sum(diffs< -1e-7)),ties=int(np.sum(np.abs(diffs)<=1e-7)),losses=int(np.sum(diffs>1e-7)),mean_RQ_minus_comparator_swaps=float(diffs.mean()),paired_bootstrap95=np.quantile(boot,[.025,.975]).tolist()))
    return result,comparisons

def analyze():
    result={};errors=[];checks=0;bound_errors=[];discover=[];confirm=[]
    for phase in ['discovery','confirmation']:
        ds=json.loads((OUT/(phase+'_details.json')).read_text(encoding='utf-8'));table=read(phase+'.csv');lookup={(r['family'],int(r['n']),int(r['seed']),r['method']):r for r in table}
        if phase=='discovery':discover=ds
        else:confirm=ds
        for d in ds:
            K=np.array(d['K']);assert np.max(np.abs(K-K.T))<1e-12 and np.max(np.abs(np.diag(K)))<1e-12
            for name,o in d['outputs'].items():
                z=np.array(o['z']);assert set(z)<=set([-1,1]) and len(z)==len(K)
                value=sum(K[i,j]*z[i]*z[j] for i in range(len(K)) for j in range(i+1,len(K)))
                row=lookup[(d['family'],d['n'],d['seed'],name)];errors += [abs(value-o['score']),abs(value-float(row['score']))];checks+=1
                if d['exact_optimum'] is not None:assert value<=d['exact_optimum']+1e-6
                if d['family']=='paintshop':assert abs(paint_cost(d['metadata'],z)-(d['metadata']['paint_cost_constant']-value))<1e-7
                if name=='sdp1024' and d['exact_optimum'] is not None:bound_errors.append(d['exact_optimum']-o['info']['relaxation'])
            if d['exact_optimum'] is not None:
                q=d['standard_qaoa'];assert q['expectation']<=q['best256']+1e-6 and q['best256']<=q['best_decision_matched']+1e-6 and q['best_decision_matched']<=d['exact_optimum']+1e-6
        summaries=[]
        for family in sorted({d['family'] for d in ds}):
            for n in sorted({d['n'] for d in ds if d['family']==family}):
                a=[d for d in ds if d['family']==family and d['n']==n];rq=lambda d:d['outputs']['rqaoa_p1']['score']
                entry=dict(family=family,n=n,cases=len(a))
                if n<=20:
                    entry.update(rq_exact=sum(abs(rq(d)-d['exact_optimum'])<1e-7 for d in a),rq_minus_qaoa256_normalized=mean((rq(d)-d['standard_qaoa']['best256'])/(2*d['weight_sum']) for d in a),rq_minus_qaoa_decision_matched_normalized=mean((rq(d)-d['standard_qaoa']['best_decision_matched'])/(2*d['weight_sum']) for d in a))
                entry['comparators']={}
                for m in ['local32','sdp1024','sdp1024_local','local_cpu_matched']:
                    diffs=[rq(d)-d['outputs'][m]['score'] for d in a]
                    entry['comparators'][m]=dict(wins=sum(x>1e-7 for x in diffs),ties=sum(abs(x)<=1e-7 for x in diffs),losses=sum(x< -1e-7 for x in diffs))
                summaries.append(entry)
        result[phase]=summaries
    fs,fc=fresh_summary();result['fresh_paint_summary']=fs;result['fresh_paint_comparisons']=fc
    # Rescore frozen-transfer and fresh-challenge stored solutions independently.
    td=json.loads((OUT/'transfer_details.json').read_text());lookup={(d['family'],d['n'],d['seed']):d for d in confirm};tr=read('transfer.csv')
    for d in td:
        parent=lookup[(d['family'],d['n'],d['seed'])];K=np.array(parent['K'])
        for field,m in [('transfer_z','rqaoa_transfer'),('polished_z','rqaoa_transfer_local'),('matched_local_z','local_transfer_cpu_matched')]:
            z=np.array(d[field]);value=sum(K[i,j]*z[i]*z[j] for i in range(len(K)) for j in range(i+1,len(K)));row=next(r for r in tr if r['family']==d['family'] and int(r['n'])==d['n'] and int(r['seed'])==d['seed'] and r['method']==m);errors.append(abs(value-float(row['score'])));checks+=1
    fresh=json.loads((OUT/'fast_confirmation_details.json').read_text())
    for d in fresh:
        K=np.array(d['K'])
        for o in d['outputs'].values():
            z=np.array(o['z']);v=sum(K[i,j]*z[i]*z[j] for i in range(len(K)) for j in range(i+1,len(K)));errors.append(abs(v-o['score']));assert abs(paint_cost(d['metadata'],z)-(d['metadata']['paint_cost_constant']-v))<1e-7;checks+=1
    if (OUT/'dense128_details.json').exists():
        dense=json.loads((OUT/'dense128_details.json').read_text());dr=read('dense128.csv');result['dense128']=[]
        for d in dense:
            K=np.array(d['K'])
            for m,o in d['outputs'].items():
                z=np.array(o['z']);v=sum(K[i,j]*z[i]*z[j] for i in range(len(K)) for j in range(i+1,len(K)));r=next(r for r in dr if float(r['density'])==d['density'] and int(r['seed'])==d['seed'] and r['method']==m);errors.append(abs(v-float(r['score'])));checks+=1
            rqscore=float(next(r['score'] for r in dr if float(r['density'])==d['density'] and int(r['seed'])==d['seed'] and r['method']=='rqaoa_first_peak'))
            result['dense128'].append(dict(density=d['density'],seed=d['seed'],RQscore=rqscore,comparisons={m:float(next(r['score'] for r in dr if float(r['density'])==d['density'] and int(r['seed'])==d['seed'] and r['method']==m))-rqscore for m in ['sdp1024','sdp1024_local','local32','local_cpu_matched','tabu_cpu_matched']}))
    result['audit']=dict(rescored_outputs=checks,max_rescore_error=max(errors),max_small_SDP_bound_violation=max(bound_errors),note='Numerical SDP not a formal certificate; large references only feasible incumbents')
    assert max(errors)<1e-7
    result['file_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*.csv')}
    (OUT/'analysis.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result['audit']));print('Fresh summaries',json.dumps(fs))

if __name__=='__main__':analyze()

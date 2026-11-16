"""Fresh held-out paintshop challenge after identifying fast-transfer signal.
Frozen transfer trained on initial discovery only; 40 new instances, no tuning.
"""
import json,time,statistics
import numpy as np
from sweep import OUT,make_instance,score,local,contract,paint_cost,sdp_polish,timed_local,dumpcsv
from transfer import solve
from paint_baselines import author_solve,AUTHOR_RSG,recursive_greedy

def timed_tabu(K,seed,budget):
    rng=np.random.default_rng(seed);n=len(K);z=rng.choice([-1,1],n);current=score(K,z);best=current;bestz=z.copy();gain=-2*z*(K@z);until=np.zeros(n,int);step=0;since=0;restarts=0;st=time.perf_counter()
    while time.perf_counter()-st<budget or step==0:
        allowed=(until<=step)|(current+gain>best+1e-9)
        k=int(np.argmax(np.where(allowed,gain,-np.inf)))
        old=z[k];g=float(gain[k]);current+=g;z[k]*=-1
        gain+=4*z*K[:,k]*old;gain[k]=-g
        until[k]=step+5+int(rng.integers(5));step+=1;since+=1
        if current>best+1e-8:best=current;bestz=z.copy();since=0
        if step%128==0:
            current=score(K,z);gain=-2*z*(K@z)
        if since>10*n:
            z=rng.choice([-1,1],n);current=score(K,z);gain=-2*z*(K@z);until[:]=0;since=0;restarts+=1
    return bestz,dict(moves=step,restarts=restarts,actual_seconds=time.perf_counter()-st)

def validate_tabu():
    rng=np.random.default_rng(1112);K=rng.normal(size=(8,8));K=np.triu(K,1);K+=K.T;z=rng.choice([-1,1],8);gain=-2*z*(K@z);errors=[]
    for step in range(30):
        k=int(rng.integers(8));old=z[k];g=gain[k];before=score(K,z);z[k]*=-1;gain+=4*z*K[:,k]*old;gain[k]=-g
        errors += [abs(score(K,z)-before-g),float(np.max(np.abs(gain+2*z*(K@z))))]
    assert max(errors)<1e-10
    return float(max(errors))

def run():
    model=json.loads((OUT/'transfer_model.json').read_text(encoding='utf-8'));param=model['parameters']['paintshop']
    protocol=dict(reason='New 40-instance challenge after fast RQAOA transfer showed signal in first confirmation; no refitting',family='paintshop',sizes_and_seeds={'32':list(range(100,120)),'64':list(range(200,220))},transfer=param,quantum_evolution_depth=1,correlations='exact CPU',repetitions='transfer median of3 times; local/tabu each3 seeds at same per-case median CPU budget',classical=['RSG author code','canonical RG','edge contraction','local32','SDP1024 plus local polish','CPU-matched local','CPU-matched tabu'],tabu='tenure5..9; aspiration; restart after10*n nonimproving moves; no target-instance tuning',optimality='No exact certificate for n32/64; report swaps directly',tabu_delta_validation=validate_tabu())
    (OUT/'fast_confirmation_protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    rows=[];details=[]
    for n,seeds in [(32,range(100,120)),(64,range(200,220))]:
      for seed in seeds:
        K,meta=make_instance('paintshop',n,seed);times=[];zs=[]
        for rep in range(3):
            st=time.perf_counter();z,trace=solve(K,param);times.append(time.perf_counter()-st);zs.append(z)
        assert all(np.array_equal(z,zs[0]) for z in zs);budget=statistics.median(times);outputs={}
        def add(name,z,seconds,rep=-1,info=None):
            value=score(K,z);changes=paint_cost(meta,z);assert abs(changes-(meta['paint_cost_constant']-value))<1e-9
            rows.append(dict(n=n,seed=seed,method=name,rep=rep,score=value,paint_changes=changes,seconds=seconds,budget_seconds=budget))
            outputs[name+str(rep)]=dict(z=z.astype(int).tolist(),score=value,info=info)
        add('rqaoa_transfer',zs[0],budget,info=dict(times=times,trace=trace))
        st=time.perf_counter();z,info=local(K,seed,32);add('local32',z,time.perf_counter()-st,info=info)
        st=time.perf_counter();z,info=contract(K,'edge');add('edge_contract',z,time.perf_counter()-st,info=info)
        st=time.perf_counter();z=recursive_greedy(meta);add('paint_RG',z,time.perf_counter()-st)
        st=time.perf_counter();z=author_solve(meta,AUTHOR_RSG);add('paint_RSG_author',z,time.perf_counter()-st)
        st=time.perf_counter();raw,polished,info=sdp_polish(K,seed);sdpseconds=time.perf_counter()-st;add('sdp1024',raw,info['raw_seconds'],info=info);add('sdp1024_local',polished,sdpseconds,info=info)
        for rep in range(3):
            st=time.perf_counter();z,info=timed_local(K,seed*100+rep,budget);add('local_cpu_matched',z,time.perf_counter()-st,rep,info)
            st=time.perf_counter();z,info=timed_tabu(K,seed*1000+rep,budget);add('tabu_cpu_matched',z,time.perf_counter()-st,rep,info)
        details.append(dict(n=n,seed=seed,K=K.tolist(),metadata=meta,outputs=outputs))
        dumpcsv(OUT/'fast_confirmation.csv',rows);(OUT/'fast_confirmation_details.json').write_text(json.dumps(details),encoding='utf-8')
        case=[r for r in rows if r['n']==n and r['seed']==seed];rq=next(r['paint_changes'] for r in case if r['method']=='rqaoa_transfer');tabu=statistics.mean(r['paint_changes'] for r in case if r['method']=='tabu_cpu_matched');print('fresh paint',n,seed,'RQ',rq,'tabu mean',round(tabu,3),'budget',round(budget,4),flush=True)
    print('Fresh cases',len(details),'rows',len(rows),flush=True)

if __name__=='__main__':run()

"""Question9 supplementary p=1..3 and small Aer gate-noise experiments."""
import json,sys,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from qiskit import QuantumCircuit,transpile
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel,depolarizing_error,ReadoutError

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'research'))
from rqaoa_experiment import optimize,cuts,spins,circuit
from q9_study import OUT,family_graph,dumpcsv,best_of
import networkx as nx

def state(cost,parameters):
    p=len(parameters)//2;n=int(np.log2(len(cost)));v=np.ones(len(cost),complex)/np.sqrt(len(cost))
    for ga,be in zip(parameters[:p],parameters[p:]):
        v*=np.exp(-1j*ga*cost);c=np.cos(be);s=-1j*np.sin(be)
        for q in range(n):
            block=v.reshape(-1,2,1<<q);lo=block[:,0,:].copy();hi=block[:,1,:].copy()
            block[:,0,:]=c*lo+s*hi;block[:,1,:]=s*lo+c*hi
    return v

def qiskit_layers(w,x):
    p=len(x)//2;c=QuantumCircuit(len(w));c.h(range(len(w)))
    for ga,be in zip(x[:p],x[p:]):
        for i,j in zip(*np.triu_indices(len(w),1)):
            if w[i,j]:c.rzz(-float(ga*w[i,j]),int(i),int(j))
        c.rx(2*float(be),range(len(w)))
    return c

def validate():
    rng=np.random.default_rng(4401);errors=[]
    for p in (1,2,3):
        g=family_graph('random',7,2);w=nx.to_numpy_array(g);cost=cuts(w,spins(7));x=rng.uniform(.1,1,size=2*p)
        errors.append(float(np.max(np.abs(np.abs(state(cost,x))**2-Statevector.from_instruction(qiskit_layers(w,x)).probabilities()))))
    assert max(errors)<1e-10
    return max(errors)

def depth():
    error=validate();rows=[];st=time.perf_counter();families=['path','star','cycle','tree','complete','bipartite','community','random']
    for family in families:
        g=family_graph(family,10,0);w=nx.to_numpy_array(g);cost=cuts(w,spins(10));opt=float(cost.max());f=optimize(w)
        for seed in range(3):
            previous=np.array([f['gamma'],f['beta']]);prevval=f['expectation']
            for p in (1,2,3):
                start=time.perf_counter();calls=0
                if p==1:x=previous.copy();value=prevval;calls=f['evaluations']
                else:
                    base=np.r_[previous[:p-1],0.,previous[p-1:],0.];x=base.copy();value=prevval
                    rng=np.random.default_rng(86000+seed*100+p)
                    grid=np.linspace(0,1,p);oldgrid=np.linspace(0,1,p-1)
                    warm=np.r_[np.interp(grid,oldgrid,previous[:p-1]),np.interp(grid,oldgrid,previous[p-1:])]*(p-1)/p
                    init=warm+rng.normal(0,.08,size=2*p)
                    def objective(params):
                        nonlocal calls
                        calls+=1;pr=np.abs(state(cost,params))**2;return -float(pr@cost)
                    result=minimize(objective,init,method='L-BFGS-B',options={'maxfun':800,'maxiter':110,'ftol':1e-11,'gtol':1e-6})
                    if -result.fun>value:x=result.x;value=float(-result.fun)
                pr=np.abs(state(cost,x))**2
                rows.append(dict(family=family,n=10,graph_seed=0,optimizer_seed=seed,p=p,optimum=opt,edges=g.number_of_edges(),expected_ratio=value/opt,success_probability=float(pr[cost>=opt-1e-9].sum()),best256_ratio=best_of(pr,cost)/opt,evaluations=calls,seconds=time.perf_counter()-start,parameters=json.dumps(x.tolist())))
                previous=x;prevval=value
        dumpcsv(OUT/'depth.csv',rows);print('depth',family,flush=True)
    (OUT/'depth_audit.json').write_text(json.dumps(dict(rows=len(rows),max_numpy_vs_qiskit_probability_error=error,seconds=time.perf_counter()-st,p1_note='Deterministic analytic optimizer repeated across optimizer-seed labels, not independent random p1 training.',higher_depth_budget='Nominal maxfun 800 per optimizer seed; scipy may finish a finite-difference batch after the cap. Actual calls recorded in CSV.',optimizer='L-BFGS-B with interpolated warm starts; post-pilot refinement, retained Powell pilot separately',rqaoa_depth='RQAOA remains p=1; depth panel is standard QAOA only'),indent=2),encoding='utf-8')

def noisy_samples(w,fit,error,seed,shots=4096):
    model=NoiseModel()
    if error:
        model.add_all_qubit_quantum_error(depolarizing_error(error,2),['cx'])
        model.add_all_qubit_readout_error(ReadoutError([[.99,.01],[.01,.99]]))
    simulator=AerSimulator(method='density_matrix',noise_model=model,max_parallel_threads=1)
    c=circuit(w,fit['gamma'],fit['beta']);c.measure_all()
    tc=transpile(c,basis_gates=['h','rx','rz','cx','measure'],optimization_level=0,seed_transpiler=1107)
    counts=simulator.run(tc,shots=shots,seed_simulator=seed).result().get_counts()
    zz=[]
    for bitstring,count in counts.items():zz.extend([[1 if b=='0' else -1 for b in bitstring[::-1]]]*count)
    zz=np.array(zz,float)
    # get_counts groups identical outcomes. Shuffle before forming best-of-256
    # batches; grouped count order is not a chronological measurement stream.
    np.random.default_rng(seed+900000).shuffle(zz,axis=0)
    return zz,dict(cx_count=int(tc.count_ops().get('cx',0)),depth=tc.depth())

def noise():
    rows=[];start=time.perf_counter()
    for family in ('cycle','bipartite','random'):
        g=family_graph(family,8,0);original=nx.to_numpy_array(g);v=cuts(original,spins(8));opt=float(v.max())
        for error in (0.,.005,.02):
            for seed in range(5):
                w=original.copy();fit=optimize(w);samples,gate=noisy_samples(w,fit,error,1000+seed)
                score=cuts(original,samples);rawmean=float(score.mean()/opt);success=float(np.mean(score>=opt-1e-9));batches=score.reshape(16,256).max(axis=1)
                active=list(range(8));trace=[];totalcx=0;totalshots=0
                while len(w)>4 and np.abs(w).max()>1e-12:
                    fit=optimize(w);samples,resources=noisy_samples(w,fit,error,2000+100*seed+len(w));corr=samples.T@samples/len(samples)
                    i,j=np.triu_indices(len(w),1);valid=np.abs(w[i,j])>1e-12;i,j=i[valid],j[valid];mag=np.abs(corr[i,j]);k=int(np.flatnonzero(mag>=mag.max()-1e-9)[0]);u,vv=int(i[k]),int(j[k]);s=1 if corr[u,vv]>=0 else -1
                    trace.append((active[u],active[vv],s));totalcx+=resources['cx_count'];totalshots+=4096
                    w=w.copy()
                    for a in range(len(w)):
                        if a not in (u,vv):w[u,a]=w[a,u]=w[u,a]+s*w[vv,a]
                    w=np.delete(np.delete(w,vv,0),vv,1);active.pop(vv)
                zr=spins(len(w));val=cuts(w,zr);solution=dict(zip(active,zr[int(np.argmax(val))]))
                for i,j,s in reversed(trace):solution[j]=s*solution[i]
                decoded=np.array([[solution[i] for i in range(8)]]);rqcut=float(cuts(original,decoded)[0])
                rows.append(dict(family=family,n=8,graph_seed=0,two_qubit_depolarizing_parameter=error,readout_error=.01 if error else 0.,measurement_seed=seed,shots_per_circuit=4096,qaoa_expected_empirical_ratio=rawmean,qaoa_success_empirical_probability=success,qaoa_best256_batch_ratio=float(batches.mean()/opt),edge_rq_noisy_ratio=rqcut/opt,qaoa_cx=gate['cx_count'],qaoa_depth=gate['depth'],rq_total_cx_per_once_each_stage=totalcx,rq_total_evaluation_shots=totalshots))
            dumpcsv(OUT/'noise.csv',rows);print('noise',family,error,flush=True)
    (OUT/'noise_audit.json').write_text(json.dumps(dict(rows=len(rows),seconds=time.perf_counter()-start,backend='Aer density_matrix',parameter_training='ideal analytic p1 at every recursive stage',noise='Depolarizing parameter on every decomposed CX; 1% readout when parameter>0; no calibrated hardware, T1/T2, routing or noisy training',rq_variant='edge-only unrestricted domain'),indent=2),encoding='utf-8')

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='noise':noise()
    else:depth()

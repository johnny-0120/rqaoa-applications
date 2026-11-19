"""11/18 check teammate's claims + test the conditions he proposed.

Part A (step0): on the ORIGINAL graph (first RQAOA step), per signal:
  copyJ   : fraction of edges where sign(signal) == sign(J)          (claim: qaoa p1,p2 = 100%)
  signacc : fraction of edges where sign(signal) == ground-state s_u s_v   (claim: p1=p2=0.870, local_r2=0.911)
  top_ok  : the edge RQAOA would lock first (max |signal|) has the correct sign
  frust_hit: among edges the ground state VIOLATES (s_u s_v != sign J), fraction the signal gets right
  gamma*,beta : optimised QAOA angles (gamma scaled by max|J|) -> 'small angle' claim
Part B (lock): full RQAOA lock rate on new graph families (triangular, 3D cubic, random 3-regular).
Usage: python claims.py step0 <family> <couplings> ; python claims.py lock <family> <couplings>
"""
import sys,os,csv,time
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('MKL_NUM_THREADS','1')
import numpy as np
from multiprocessing import Pool
from sweep import OUT,dumpcsv
from depth import all_spins,costs,evolve,qaoa_params
import signals as S

def weights(rng,m,couplings):
    return rng.choice([-1.,1.],m) if couplings=='pm' else rng.normal(size=m)

def graph(family,seed,couplings):
    rng=np.random.default_rng(seed);E=[]
    if family in('sq44','sq45','tri44'):
        L,W={'sq44':(4,4),'sq45':(4,5),'tri44':(4,4)}[family]
        for x in range(L):
            for y in range(W):
                i=x*W+y
                if x+1<L:E.append((i,i+W))
                if y+1<W:E.append((i,i+1))
                if family=='tri44' and x+1<L and y+1<W:E.append((i,i+W+1))
        n=L*W
    elif family=='cub224':                       # 2x2x4 open cubic lattice, non-planar-ish 3D
        dims=(2,2,4);n=16;idx=lambda a,b,c:(a*2+b)*4+c
        for a in range(2):
            for b in range(2):
                for c in range(4):
                    i=idx(a,b,c)
                    if a+1<2:E.append((i,idx(a+1,b,c)))
                    if b+1<2:E.append((i,idx(a,b+1,c)))
                    if c+1<4:E.append((i,idx(a,b,c+1)))
    elif family=='cub233':                       # 2x3x3 = 18 spins, genuinely 3D (non-planar)
        n=18;idx=lambda a,b,c:(a*3+b)*3+c
        for a in range(2):
            for b in range(3):
                for c in range(3):
                    i=idx(a,b,c)
                    if a+1<2:E.append((i,idx(a+1,b,c)))
                    if b+1<3:E.append((i,idx(a,b+1,c)))
                    if c+1<3:E.append((i,idx(a,b,c+1)))
    elif family=='reg3_16':
        import networkx as nx
        G=nx.random_regular_graph(3,16,seed=int(seed));n=16;E=[tuple(sorted(e)) for e in G.edges()]
    else:raise ValueError(family)
    K=np.zeros((n,n));w=weights(rng,len(E),couplings)
    for (i,j),x in zip(E,w):K[i,j]=K[j,i]=x
    return K

def ground_states(K):
    Z=all_spins(len(K));c=costs(K,Z);return Z[c>=c.max()-1e-9]

QSIG=['qaoa_p1','qaoa_p2','qaoa_p3']
CSIG=['naive','local_r1','local_r2','local_r3','sdp']

_orig_signal=S.signal
def sig(name,K,rng,cache=None):
    if name=='naive':return {(u,v):float(K[u,v]) for u,v in S.edges(K)}   # lock the strongest coupling, follow its sign
    return _orig_signal(name,K,rng,cache)
S.signal=sig   # run_rqaoa looks up S.signal, so it now also understands 'naive'

def step0_one(args):
    family,couplings,seed,names=args;K=graph(family,seed,couplings);rng=np.random.default_rng(seed*7+1)
    E=S.edges(K);G=ground_states(K)
    # edge "truth": +1/-1 if every ground state agrees, 0 if ground states disagree (degenerate, pm)
    tr={}
    for u,v in E:
        vals=set((G[:,u]*G[:,v]).tolist());tr[(u,v)]=vals.pop() if len(vals)==1 else 0
    rows=[];n=len(K);Z=all_spins(n);c=costs(K,Z)
    for name in names:
        t=time.perf_counter();extra=dict(gamma_scaled='',beta='',approx_ratio='')
        if name.startswith('qaoa_p'):
            p=int(name[-1]);x,e=qaoa_params(c,n,p,rng);pr=np.abs(evolve(c,x,n))**2
            s={(u,v):float(np.sum(pr*Z[:,u]*Z[:,v])) for u,v in E}
            sc=np.abs(K).max();extra=dict(gamma_scaled=';'.join(f'{g*sc:.3f}' for g in x[:p]),beta=';'.join(f'{b:.3f}' for b in x[p:]),
                                          approx_ratio=round((e-c.min())/(c.max()-c.min()),4))
        else:s=sig(name,K,rng)
        secs=time.perf_counter()-t
        sg={e:np.sign(v) for e,v in s.items()}
        known=[e for e in E if tr[e]!=0];viol=[e for e in known if tr[e]!=np.sign(K[e])]
        top=max(E,key=lambda e:abs(s[e]))
        rows.append(dict(family=family,couplings=couplings,seed=seed,n=n,signal=name,
            copyJ=round(np.mean([sg[e]==np.sign(K[e]) for e in E]),4),
            signacc=round(np.mean([sg[e]==tr[e] for e in known]),4) if known else '',
            n_viol=len(viol),frust_hit=round(np.mean([sg[e]==tr[e] for e in viol]),4) if viol else '',
            top_ok=int(tr[top]==0 or sg[top]==tr[top]),max_abs=round(max(abs(v) for v in s.values()),4),seconds=round(secs,2),**extra))
    return rows

def lock_one(args):
    family,couplings,seed,name=args;K=graph(family,seed,couplings)
    r=S.run_rqaoa(K,name,seed*31+len(name));r.update(family=family,couplings=couplings,seed=seed,n=len(K),signal=name);return r

def main(mode,family,couplings,nseeds,names,procs):
    path=OUT/f'claims_{mode}_{family}_{couplings}.csv';rows=[]
    if path.exists():rows=list(csv.DictReader(open(path,encoding='utf-8-sig')))
    seeds=range(5000,5000+nseeds)
    if mode=='step0':
        done={int(r['seed']) for r in rows};jobs=[(family,couplings,s,names) for s in seeds if s not in done]
        fn=step0_one
    else:
        done={(int(r['seed']),r['signal']) for r in rows};jobs=[(family,couplings,s,nm) for s in seeds for nm in names if (s,nm) not in done]
        jobs.sort(key=lambda j:-int(j[3][-1]) if j[3].startswith('qaoa') else 0)   # slowest jobs first
        fn=lock_one
    with Pool(procs) as pool:
        for out in pool.imap_unordered(fn,jobs):
            rows+=out if isinstance(out,list) else [out];dumpcsv(path,rows)
            print(mode,family,couplings,len(rows),flush=True)

if __name__=='__main__':
    mode,family,couplings=sys.argv[1:4];nseeds=int(sys.argv[4]) if len(sys.argv)>4 else 30
    names=sys.argv[5].split(',') if len(sys.argv)>5 else (['random']+CSIG+QSIG if mode=='lock' else CSIG+QSIG)
    main(mode,family,couplings,nseeds,names,int(os.environ.get('PROCS','11')))

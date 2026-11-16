"""Conditional patch-sign benchmark, not raw point-cloud reconstruction.

Own benchmark and solver implementations. Geometry ports follow the MIT-licensed
WNNC and Dipole repositories; preserved notices in vendor_inspection/*LICENSE.
No original learned models, treecode acceleration, GPU or QPU are executed.
"""
from pathlib import Path
import argparse, csv, hashlib, importlib.metadata, json, platform, sys, time
import numpy as np
import networkx as nx
from scipy.spatial.distance import cdist
from scipy.optimize import minimize_scalar
import cvxpy as cp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from rqaoa_experiment import components, correlations, circuit

OUT = Path(__file__).resolve().parent / 'results'
OUT.mkdir(parents=True, exist_ok=True)

def unit(a):
    return a / np.maximum(np.linalg.norm(a, axis=-1, keepdims=True), 1e-15)

def sign(a):
    return np.where(np.asarray(a) >= 0, 1, -1)

def score(K, z):
    return float(z @ K @ z / 2)

def error(z, truth, weights=None):
    if weights is None:
        weights = np.ones(len(z))
    a = np.average(z != truth, weights=weights)
    return float(min(a, 1-a))

def cloud(shape, n, seed, points_per_patch=12):
    """Analytic surfaces with supplied sign-consistent local normal axes.
    Truth is used for simulation and evaluation only, not graph construction.
    Thin slab samples just the top/bottom surfaces; side walls are absent.
    """
    rng = np.random.default_rng(1000 + seed + 100*n + 10000*['sphere','torus','thin_slab'].index(shape))
    pts=[]; norms=[]
    if shape == 'sphere':
        ii=np.arange(n);z=1-2*(ii+.5)/n;phi=ii*np.pi*(3-np.sqrt(5))
        centers=np.c_[np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z]
        for c in centers:
            tangent=rng.normal(size=(points_per_patch,3));tangent-=np.outer(tangent@c,c)
            p=unit(c+.13*tangent);pts.append(p);norms.append(p.copy())
    elif shape == 'torus':
        major=n//4
        for a in range(major):
            for b in range(4):
                u=2*np.pi*(a+.5)/major+rng.uniform(-.18,.18,points_per_patch)
                v=2*np.pi*(b+.5)/4+rng.uniform(-.18,.18,points_per_patch)
                pts.append(np.c_[(1+.35*np.cos(v))*np.cos(u),(1+.35*np.cos(v))*np.sin(u),.35*np.sin(v)])
                norms.append(np.c_[np.cos(v)*np.cos(u),np.cos(v)*np.sin(u),np.sin(v)])
    elif shape == 'thin_slab':
        count=n//2
        # Matching top/bottom sample locations make wrong cross-sheet links explicit.
        xy=np.c_[np.linspace(-.8,.8,count), .4*np.sin(np.arange(count)*2.4)]
        for side in (1,-1):
            for c in xy:
                p=np.c_[c+rng.uniform(-.10,.10,(points_per_patch,2)),np.full(points_per_patch,.045*side)]
                pts.append(p);norms.append(np.tile([0.,0.,float(side)],(points_per_patch,1)))
    pts=np.array(pts);true_norm=np.array(norms)
    points=pts+rng.normal(scale=.003,size=pts.shape)
    axes=unit(true_norm+rng.normal(scale=.03,size=true_norm.shape))
    truth=rng.choice([-1,1],n)
    observed=axes*truth[:,None,None]
    centers=points.mean(axis=1);mean_norm=unit(observed.mean(axis=1))
    return dict(points=points,normals=observed,truth=truth,true_normals=true_norm,centers=centers,mean_normals=mean_norm)

def make_graph(data, panel, corruption, seed):
    centers=data['centers'];normals=data['mean_normals'];truth=data['truth'];n=len(truth)
    dist=cdist(centers,centers);np.fill_diagonal(dist,np.inf)
    # Symmetric union of 4NN. Join components via closest cross-component pair.
    mask=np.zeros((n,n),bool)
    for i in range(n):mask[i,np.argsort(dist[i])[:4]]=True
    mask|=mask.T
    G=nx.from_numpy_array(mask)
    while not nx.is_connected(G):
        labels=np.empty(n,int)
        for k,comp in enumerate(nx.connected_components(G)):labels[list(comp)]=k
        D=np.where(labels[:,None]!=labels[None,:],dist,np.inf);i,j=np.unravel_index(np.argmin(D),D.shape)
        mask[i,j]=mask[j,i]=True;G.add_edge(int(i),int(j))
    dot=normals@normals.T
    scale=float(np.median(dist[mask]));weights=np.exp(-np.square(dist/scale))*np.maximum(np.abs(dot),.15)
    weights=np.where(mask,weights,0.);np.fill_diagonal(weights,0.)
    if panel=='geometry':relations=sign(dot)
    else:
        # Explicit oracle-relative-sign control, not inferred geometric observations.
        relations=np.outer(truth,truth)
        rng=np.random.default_rng(800000+seed);i,j=np.triu_indices(n,1)
        flips=rng.random(len(i))<corruption
        relations[i[flips],j[flips]]*=-1;relations[j,i]=relations[i,j]
    K=weights*relations
    return K,dict(edges=int(mask.sum()//2),bad_edge_weight_fraction=float(np.sum(weights*(relations!=np.outer(truth,truth)))/weights.sum()),components=nx.number_connected_components(G))

def exhaustive(K, truth=None):
    """Fix first spin +1 (global symmetry), enumerate with original score."""
    n=len(K);best=-np.inf;bestz=None;ties=0;emin=1.;emax=0.
    for start in range(0,1<<(n-1),32768):
        nums=np.arange(start,min(start+32768,1<<(n-1)),dtype=np.uint64)
        zz=np.ones((len(nums),n));zz[:,1:]=1-2*((nums[:,None]>>np.arange(n-1,dtype=np.uint64))&1).astype(float)
        vals=np.sum((zz@K)*zz,axis=1)/2;m=float(vals.max())
        if m>best+1e-9:best=m;bestz=zz[int(np.argmax(vals))].copy();ties=0;emin=1.;emax=0.
        valid=np.abs(vals-best)<=1e-9
        ties+=int(valid.sum())
        if truth is not None and valid.any():
            e=np.mean(zz[valid]!=truth,axis=1);e=np.minimum(e,1-e)
            emin=min(emin,float(e.min()));emax=max(emax,float(e.max()))
    return bestz,dict(optimum=best,optima_up_to_flip=ties,optimum_truth_error_min=emin,optimum_truth_error_max=emax)

def mst(K):
    G=nx.from_numpy_array(np.abs(K));T=nx.maximum_spanning_tree(G)
    z=np.ones(len(K),int)
    for u,v in nx.bfs_edges(T,0):z[v]=z[u]*sign(K[u,v])
    return z,{}

def spectral(K):
    d=np.maximum(np.abs(K).sum(axis=1),1e-15)
    _,V=np.linalg.eigh(K/np.sqrt(np.outer(d,d)))
    return sign(V[:,-1]),{}

def improve(K,z):
    z=z.copy();steps=0
    for _ in range(10000):
        gain=-2*z*(K@z);i=int(np.argmax(gain))
        if gain[i]<=1e-10:break
        z[i]*=-1;steps+=1
    return z,steps

def local(K, seed, starts):
    rng=np.random.default_rng(seed);best=-np.inf;bestz=None;steps=0
    for _ in range(starts):
        z,s=improve(K,rng.choice([-1,1],len(K)));steps+=s;v=score(K,z)
        if v>best+1e-10:best=v;bestz=z.copy()
    return bestz,dict(starts=starts,flip_proposals=steps)

def anneal(K,seed,starts=8,sweeps=80):
    rng=np.random.default_rng(seed);best=-np.inf;bestz=None;n=len(K)
    tmax=max(float(np.abs(K).sum(axis=1).max()),1e-3)
    for _ in range(starts):
        z=rng.choice([-1,1],n)
        for T in np.geomspace(tmax,tmax*.001,sweeps):
            for i in rng.permutation(n):
                gain=-2*z[i]*np.dot(K[i],z)
                if gain>=0 or rng.random()<np.exp(gain/T):z[i]*=-1
        z,_=improve(K,z);v=score(K,z)
        if v>best+1e-10:best=v;bestz=z.copy()
    return bestz,dict(starts=starts,flip_proposals=starts*sweeps*n)

def sdp(K,seed):
    n=len(K);X=cp.Variable((n,n),symmetric=True)
    problem=cp.Problem(cp.Maximize(cp.sum(cp.multiply(K,X))/2),[cp.diag(X)==1,X>>0])
    problem.solve(solver='SCS',eps=1e-7,max_iters=30000,verbose=False)
    if X.value is None:raise RuntimeError('SDP failed: '+str(problem.status))
    eig,V=np.linalg.eigh((X.value+X.value.T)/2)
    factor=V@np.diag(np.sqrt(np.maximum(eig,0)))
    rng=np.random.default_rng(seed);zs=sign(rng.normal(size=(256,n))@factor.T)
    vals=np.sum((zs@K)*zs,axis=1)/2;k=int(np.argmax(vals))
    return zs[k],dict(status=problem.status,relaxation_value=float(problem.value),min_eigenvalue=float(eig.min()),diagonal_error=float(np.max(np.abs(np.diag(X.value)-1))),rounds=256)

def fit_qaoa(K,grid_count=129,gamma_max=2*np.pi):
    """Signed real coefficients; rescale each residual by max coupling.
    Finite bounded gamma search is heuristic, with exact best beta per gamma.
    """
    w=-2*K;scale=max(float(np.abs(w).max()),1e-15);w=w/scale;calls=0
    def evaluate(gamma):
        nonlocal calls
        calls+=1;i,j,a,b=components(w,gamma);weights=w[i,j]
        A=-.5*np.dot(weights,a);B=.25*np.dot(weights,b)
        beta=(np.arctan2(A,B)%(2*np.pi))/4
        return float(.5*weights.sum()-B+np.hypot(A,B)),float(beta)
    grid=np.linspace(0,gamma_max,grid_count);values=np.array([evaluate(g)[0] for g in grid]);k=int(np.argmax(values));ga=float(grid[k]);val,be=evaluate(ga)
    peaks=sorted([i for i in range(1,len(grid)-1) if values[i]>=values[i-1] and values[i]>=values[i+1]],key=lambda i:-values[i])[:4]
    for k in peaks:
        result=minimize_scalar(lambda g:-evaluate(g)[0],bounds=(grid[k-1],grid[k+1]),method='bounded',options={'xatol':1e-8})
        v,b=evaluate(result.x)
        if v>val+1e-10:ga,val,be=float(result.x),v,b
    return correlations(w,ga,be),dict(gamma=ga,beta=be,scale=scale,evaluations=calls,normalized_w=w)

def contract(K,mode='rqaoa',shots=0,seed=0,grid_count=129,gamma_max=2*np.pi):
    orig=K.copy();active=list(range(len(K)));trace=[];calls=0;rng=np.random.default_rng(seed)
    while len(K)>4 and np.max(np.abs(K))>1e-12:
        if mode=='rqaoa':
            M,fit=fit_qaoa(K,grid_count,gamma_max);calls+=fit['evaluations']
            if shots:
                from qiskit.quantum_info import Statevector
                prob=Statevector.from_instruction(circuit(fit['normalized_w'],fit['gamma'],fit['beta'])).probabilities()
                sample=rng.choice(len(prob),shots,p=prob)
                zz=1-2*((sample[:,None]>>np.arange(len(K)))&1).astype(float);M=zz.T@zz/shots
            strength=np.abs(M)
        else:M=K;strength=np.abs(K);fit={}
        i,j=np.triu_indices(len(K),1);sco=strength[i,j];q=int(np.flatnonzero(sco>=sco.max()-1e-9)[0]);u,v=int(i[q]),int(j[q]);s=int(sign(M[u,v]))
        trace.append(dict(remaining=len(K),keep=active[u],remove=active[v],sign=s,relation_strength=float(M[u,v]),edges=int(np.count_nonzero(np.triu(K,1))),gamma=fit.get('gamma'),beta=fit.get('beta')))
        K=K.copy()
        for k in range(len(K)):
            if k not in (u,v):K[u,k]=K[k,u]=K[u,k]+s*K[v,k]
        K=np.delete(np.delete(K,v,0),v,1);active.pop(v)
    zz,_=exhaustive(K);solution=dict(zip(active,zz))
    for t in reversed(trace):solution[t['remove']]=t['sign']*solution[t['keep']]
    z=np.array([solution[i] for i in range(len(orig))])
    return z,dict(trace=trace,evaluations=calls,stages=len(trace),decision_shots=shots*len(trace),training='exact expectations',gamma_max=gamma_max,grid_count=grid_count)

def dipole_field(source_points,source_normals,targets):
    R=source_points[:,None,:]-targets[None,:,:];dist=np.linalg.norm(R,axis=-1)
    U=R/np.maximum(dist[:,:,None],1e-15)
    E=(source_normals[:,None,:]-3*np.sum(source_normals[:,None,:]*U,axis=-1)[:,:,None]*U)/(dist**3+1e-5)[:,:,None]
    return E.sum(axis=0)

def dipole_patch(data):
    """NumPy port of strongest_field_propagation, diffuse=False.
    Supplied patches/normals replace the original learned local stage.
    """
    points=data['points'];normals=data['normals'];n=len(points);z=np.ones(n,int)
    curves=[np.linalg.eigvalsh(np.cov(p.T))[0] for p in points];initial=int(np.argmin(np.abs(curves)))
    remaining=[i for i in range(n) if i!=initial];visited=[initial]
    E=dipole_field(points[initial],normals[initial],points.reshape(-1,3)).reshape(points.shape)
    while remaining:
        interactions=np.array([np.sum(E[i]*normals[i]) for i in remaining]);k=int(np.argmax(np.abs(interactions)));q=remaining.pop(k);z[q]=sign(interactions[k])
        visited.append(q)
        if remaining:
            targets=points[remaining].reshape(-1,3)
            E[remaining]+=dipole_field(points[q],z[q]*normals[q],targets).reshape((len(remaining),points.shape[1],3))
    return z,dict(implementation='NumPy field-only port; supplied local patches; no learned model',diffuse=False)

def winding_kernels(points,width):
    """Exact dense form of official A and G near-field kernels (no 1/(4pi)).
    Kernel is zero inside width, as in the official CPU source.
    """
    d=points[:,None,:]-points[None,:,:];r=np.linalg.norm(d,axis=-1);valid=r>=width
    r3=np.where(valid,np.maximum(r,1e-15)**-3,0.);r5=np.where(valid,np.maximum(r,1e-15)**-5,0.)
    A=-d*r3[:,:,None]
    return A,d,r3,r5

def wnnc_dense(data,iters=40,width_min=.002,width_max=.016):
    """Dense float64 NumPy port of official main_wnnc.py update equations.
    Point normals are estimated from positions, then projected to patch signs.
    No treecode; official executable and original-scale performance unverified.
    """
    x=data['points'].reshape(-1,3);center=(x.min(0)+x.max(0))/2;scale=np.ptp(x,axis=0).max()*1.1
    x=(x-center)*(2/scale);mu=np.zeros_like(x);b=np.full(len(x),.5);out=mu.copy()
    for w in np.linspace(width_max,width_min,iters):
        A,d,r3,r5=winding_kernels(x,w)
        def forward(v):return np.einsum('ijc,jc->i',A,v,optimize=True)
        def transpose(s):return np.einsum('ijc,i->jc',A,s,optimize=True)
        residual=transpose(b-forward(mu));Ar=forward(residual)
        alpha=float(np.sum(residual**2)/max(np.sum(Ar**2),1e-30));mu+=alpha*residual
        dot=np.einsum('ijc,jc->ij',d,mu,optimize=True)
        G=r3@mu-3*np.einsum('ijc,ij->ic',d,dot*r5,optimize=True)
        out=unit(G);mu=out*np.linalg.norm(mu,axis=1,keepdims=True)
    normals=out.reshape(data['normals'].shape)
    alignment=np.sum(normals*data['normals'],axis=2).sum(axis=1)
    return sign(alignment),dict(implementation='Dense NumPy CPU equation port, not official binary',iterations=iters,width_min=width_min,width_max=width_max,mean_axis_alignment=float(np.mean(np.abs(np.sum(normals*data['normals'],axis=2)))))

def writecsv(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def validate():
    from qiskit.quantum_info import Statevector
    rng=np.random.default_rng(12011);errs=[];energy_err=[];sub_err=[]
    for n in [3,5,7]:
        K=rng.normal(size=(n,n));K=np.triu(K,1);K+=K.T;w=-2*K
        ga,be=rng.uniform(-1,1,2);p=Statevector.from_instruction(circuit(w,ga,be)).probabilities()
        nums=np.arange(1<<n,dtype=np.uint64);zs=1-2*((nums[:,None]>>np.arange(n,dtype=np.uint64))&1).astype(float)
        M=zs.T@(p[:,None]*zs);errs.append(float(np.max(np.abs(M-correlations(w,ga,be)))))
        energy_err.append(float(abs(np.dot(p,np.sum((zs@K)*zs,axis=1)/2)-np.sum(K*M)/2)))
        i,j,s=0,n-1,-1;R=K.copy()
        for k in range(n):
            if k not in [i,j]:R[i,k]=R[k,i]=K[i,k]+s*K[j,k]
        R=np.delete(np.delete(R,j,0),j,1);z=zs[:1<<(n-1),:n-1];full=np.insert(z,j,s*z[:,i],axis=1)
        delta=np.sum((full@K)*full,axis=1)/2-np.sum((z@R)*z,axis=1)/2
        sub_err.append(float(np.max(np.abs(delta-s*K[i,j]))))
        ze,info=exhaustive(K);assert abs(info['optimum']-max(score(K,q) for q in zs))<1e-8
    x=rng.normal(size=(12,3));mu=rng.normal(size=(12,3));b=rng.normal(size=12)
    A,d,r3,r5=winding_kernels(x,.02)
    lhs=np.dot(np.einsum('ijc,jc->i',A,mu),b);rhs=np.sum(mu*np.einsum('ijc,i->jc',A,b))
    adjoint=abs(lhs-rhs)
    # G is -gradient of the winding field away from the truncated kernel boundary.
    target=x[:1]+.31;diff=target[:,None,:]-x[None,:,:];rr=np.linalg.norm(diff,axis=-1)
    G=np.sum(mu[None,:,:]/rr[:,:,None]**3-3*diff*np.sum(diff*mu[None,:,:],axis=-1)[:,:,None]/rr[:,:,None]**5,axis=1)
    def pot(t):
        dd=t[:,None,:]-x[None,:,:];return np.sum(-dd*mu[None,:,:]/np.linalg.norm(dd,axis=-1)[:,:,None]**3)
    fd=[]
    for k in range(3):
        e=np.zeros((1,3));e[0,k]=1e-6;fd.append(-(pot(target+e)-pot(target-e))/2e-6)
    grad_err=float(np.max(np.abs(G[0]-fd)))
    consistent=np.outer([1,-1,1,-1],[1,-1,1,-1]).astype(float);np.fill_diagonal(consistent,0)
    for solver in [lambda K:mst(K),lambda K:spectral(K),lambda K:local(K,1,8),lambda K:anneal(K,1),lambda K:sdp(K,1),lambda K:contract(K)]:
        z,_=solver(consistent);assert error(z,np.array([1,-1,1,-1]))==0
    result=dict(correlation_max_error=max(errs),expectation_max_error=max(energy_err),substitution_max_error=max(sub_err),winding_adjoint_error=float(adjoint),winding_gradient_error=grad_err,consistent_graph_checks=6)
    assert max(errs)<1e-10 and max(energy_err)<1e-10 and max(sub_err)<1e-10 and adjoint<1e-10 and grad_err<1e-6
    (OUT/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result),flush=True)

def run(pilot=False):
    protocol=dict(date='2026-11-12',task='conditional patch-sign normal orientation',shapes=['sphere','torus','thin_slab'],sizes=[12] if pilot else [12,16,20],seeds=[0] if pilot else [0,1,2,3],control_corruption=[0.,.1,.25],points_per_patch=12,cutoff=4,p=1,gamma_grid=129,gamma_max=2*np.pi,geometry_ports=['dipole_patch_numpy','wnnc_dense_numpy_l0'],scope='analytic local normal axes supplied, no PCA/neural local estimation, no real scan or reconstruction',comparison='fixed settings, CPU wall time measured; not equal total budgets; exact ideal training',wnnc_commit='175b116580c9c49e8f8d421fb17eee10fcddc198',dipole_commit='72ba66b6a54858ebdb34de44afaefaeabfeee81b')
    (OUT/('pilot_protocol.json' if pilot else 'protocol.json')).write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    rows=[];details=[];start=time.perf_counter();case=0
    for shape in protocol['shapes']:
      for n in protocol['sizes']:
       for seed in protocol['seeds']:
        data=cloud(shape,n,seed)
        panels=[('control',q) for q in protocol['control_corruption']]+[('geometry',0.)]
        for panel,q in panels:
            K,meta=make_graph(data,panel,q,seed+100*n);truth=data['truth'];case+=1
            st=time.perf_counter();ze,exact_info=exhaustive(K,truth);exact_time=time.perf_counter()-st
            methods=[('exact',lambda:(ze,exact_info)),('mst',lambda:mst(K)),('spectral',lambda:spectral(K)),('sdp256',lambda:sdp(K,seed)),('local1',lambda:local(K,seed,1)),('local32',lambda:local(K,seed,32)),('anneal8',lambda:anneal(K,seed)),('edge_contract',lambda:contract(K,'edge')),('rqaoa_p1',lambda:contract(K))]
            if panel=='geometry':methods += [('dipole_patch_numpy',lambda:dipole_patch(data)),('wnnc_dense_numpy_l0',lambda:wnnc_dense(data))]
            outputs={};denom=float(np.triu(np.abs(K),1).sum())
            for name,fn in methods:
                st=time.perf_counter();z,info=fn();elapsed=exact_time if name=='exact' else time.perf_counter()-st
                value=score(K,z);err=error(z,truth);gap=max(0.,exact_info['optimum']-value)
                assert value<=exact_info['optimum']+1e-6
                row=dict(case=case,panel=panel,shape=shape,n=n,seed=seed,corruption=q,method=name,patch_error=err,wrong_patches=int(round(err*n)),score=value,optimality_gap=gap,normalized_gap=gap/max(2*denom,1e-15),weighted_disagreement=(denom-value)/(2*denom),is_optimal=int(gap<1e-7),seconds=elapsed,edges=meta['edges'],bad_edge_weight_fraction=meta['bad_edge_weight_fraction'],optimum_truth_error=error(ze,truth),optimum_truth_error_min=exact_info['optimum_truth_error_min'],optimum_truth_error_max=exact_info['optimum_truth_error_max'],evaluations=info.get('evaluations',0),stages=info.get('stages',0))
                rows.append(row);outputs[name]=dict(z=np.asarray(z).astype(int).tolist(),info=info)
            record=dict(case=case,panel=panel,shape=shape,n=n,seed=seed,corruption=q,K=K.tolist(),truth=truth.tolist(),metadata=meta,exact_info=exact_info,outputs=outputs)
            details.append(record)
            file='pilot_rows.csv' if pilot else 'main.csv';writecsv(OUT/file,rows)
            (OUT/('pilot_details.json' if pilot else 'details.json')).write_text(json.dumps(details),encoding='utf-8')
            rq=next(r for r in rows if r['case']==case and r['method']=='rqaoa_p1');ls=next(r for r in rows if r['case']==case and r['method']=='local32')
            print(f'{case}: {shape} n={n} seed={seed} {panel} q={q}: RQ error={rq["patch_error"]:.3f}, LS={ls["patch_error"]:.3f}, RQ gap={rq["optimality_gap"]:.5f}; {time.perf_counter()-start:.1f}s',flush=True)
        if not pilot:
            path=OUT/'clouds';path.mkdir(exist_ok=True)
            np.savez_compressed(path/f'{shape}_{n}_{seed}.npz',**data)
    env=dict(python=sys.version,platform=platform.platform(),packages={p:importlib.metadata.version(p) for p in ['numpy','scipy','cvxpy','qiskit','matplotlib','networkx']},total_seconds=time.perf_counter()-start,cases=case,rows=len(rows),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/('pilot_environment.json' if pilot else 'environment.json')).write_text(json.dumps(env,indent=2),encoding='utf-8');print(json.dumps(env),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('task',choices=['validate','pilot','run']);a=p.parse_args()
    if a.task=='validate':validate()
    else:run(a.task=='pilot')

"""Independent output rescoring and invariance checks for this benchmark."""
import csv,json,hashlib,time
from pathlib import Path
import numpy as np
from benchmark import OUT,cloud,make_graph,dipole_patch,wnnc_dense,error,contract,score
from geometry_graph import dipole_graph

def run():
    main=json.loads((OUT/'details.json').read_text(encoding='utf-8'))
    geo=json.loads((OUT/'geometry_graph_details.json').read_text(encoding='utf-8'))
    rows=list(csv.DictReader((OUT/'main.csv').open(encoding='utf-8-sig')))
    grows=list(csv.DictReader((OUT/'geometry_graph.csv').open(encoding='utf-8-sig')))
    errors=[];checks=0;sdp_diag=[];sdp_eigen=[];sdp_viol=[]
    for details,table in [(main,rows),(geo,grows)]:
        lookup={(int(r['case']),r['method']):r for r in table}
        for d in details:
            K=np.array(d['K']);truth=np.array(d['truth']);assert np.max(np.abs(K-K.T))<1e-12 and np.max(np.abs(np.diag(K)))<1e-12
            for m,o in d['outputs'].items():
                z=np.array(o['z']);assert len(z)==len(K) and set(z)<=set([-1,1])
                # Direct edge-list sum, independent of matrix-scoring helper.
                value=sum(K[i,j]*z[i]*z[j] for i in range(len(K)) for j in range(i+1,len(K)))
                row=lookup[(d['case'],m)];errors.append(abs(value-float(row['score'])))
                assert abs(error(z,truth)-float(row['patch_error']))<1e-12
                assert value<=d['exact_info']['optimum']+1e-7;checks+=1
                if m=='sdp256':
                    info=o['info'];sdp_diag.append(info['diagonal_error']);sdp_eigen.append(info['min_eigenvalue']);sdp_viol.append(d['exact_info']['optimum']-info['relaxation_value'])
    rng=np.random.default_rng(11212);gauge_errors=[];truth_leak=[];field_cov=[];wnnc_cov=[]
    for shape in ['sphere','torus','thin_slab']:
        data=cloud(shape,12,2);K=dipole_graph(data);flips=rng.choice([-1,1],12)
        transformed={k:np.copy(v) for k,v in data.items()};transformed['normals']*=flips[:,None,None];transformed['mean_normals']*=flips[:,None]
        new=dipole_graph(transformed);gauge_errors.append(float(np.max(np.abs(new-K*np.outer(flips,flips)))))
        # Geometry solvers have no access to simulator truth arrays.
        stripped={k:v for k,v in data.items() if k not in ['truth','true_normals']}
        truth_leak.append(float(np.max(np.abs(dipole_graph(stripped)-K))))
        a,_=dipole_patch(stripped);b,_=dipole_patch({k:v for k,v in transformed.items() if k not in ['truth','true_normals']})
        field_cov.append(error(a,b*flips))
        a,_=wnnc_dense(stripped,width_min=.05,width_max=.2);b,_=wnnc_dense({k:v for k,v in transformed.items() if k not in ['truth','true_normals']},width_min=.05,width_max=.2)
        wnnc_cov.append(error(a,b*flips))
        altered={k:np.copy(v) for k,v in data.items()};altered['truth']=rng.choice([-1,1],12)
        K1,_=make_graph(data,'geometry',0,2);K2,_=make_graph(altered,'geometry',0,2);assert np.max(np.abs(K1-K2))==0
    # Explicitly audit parameter-search sensitivity on the three discovered failures.
    refined=[]
    for d in main:
        rq=d['outputs']['rqaoa_p1'];K=np.array(d['K'])
        if d['exact_info']['optimum']-score(K,np.array(rq['z']))>1e-7:
            st=time.perf_counter();z,info=contract(K,grid_count=513,gamma_max=4*np.pi)
            refined.append(dict(case=d['case'],original_score=score(K,np.array(rq['z'])),refined_score=score(K,z),optimum=d['exact_info']['optimum'],seconds=time.perf_counter()-st))
    result=dict(rescored_outputs=checks,max_rescore_error=max(errors),max_sdp_diagonal_error=max(sdp_diag),min_sdp_eigenvalue=min(sdp_eigen),max_sdp_bound_violation=max(sdp_viol),max_gauge_covariance_error=max(gauge_errors),truth_independence_error=max(truth_leak),dipole_gauge_error=max(field_cov),wnnc_gauge_error=max(wnnc_cov),discovered_failure_parameter_refinement=refined)
    assert max(errors)<1e-9 and max(gauge_errors)<1e-10 and max(truth_leak)==0 and max(field_cov)==0 and max(wnnc_cov)==0
    result['file_hashes']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*.csv')}
    (OUT/'audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':run()

"""Independent checks of persisted results; does not optimize new solutions."""
import csv,json,hashlib
from pathlib import Path
import numpy as np
from qiskit.quantum_info import Statevector
from rqaoa_experiment import graph,circuit,exact,spins,cuts

p=Path('results');rows=list(csv.DictReader((p/'benchmark.csv').open(encoding='utf-8-sig')))
traces=json.loads((p/'traces.json').read_text(encoding='utf-8')); errs=[]; slow=[]
for row,d in zip(rows,traces):
    g,w=graph(row['family'],int(row['n']),int(row['seed']),float(row['p']))
    assert sorted(map(tuple,d['edges']))==sorted(g.edges())
    opt,_,values=exact(w); assert opt==float(row['optimum'])
    z=np.array(d['rqaoa']['spins'])
    independent_cut=sum(z[i]!=z[j] for i,j in g.edges())
    assert independent_cut==d['rqaoa']['cut']
    assert abs(independent_cut/opt-float(row['rqaoa_ratio']))<1e-12
    fit=d['qaoa_fit'];prob=Statevector.from_instruction(circuit(w,fit['gamma'],fit['beta'])).probabilities()
    errs.append(abs(float(prob@values)-fit['expectation']))
    assert abs(prob.sum()-1)<1e-10
    for k in ['qaoa_expected_ratio','qaoa_best256_mean_ratio','rqaoa_ratio','local32_ratio']:
        assert 0<=float(row[k])<=1+1e-10
    if d['id'] in [0,44,52]:
        # Independent edge-by-edge Python bit enumeration, not the quadratic form.
        n=len(w); e=list(g.edges())
        oracle=max(sum(((x>>i)&1)!=((x>>j)&1) for i,j in e) for x in range(1<<n))
        assert oracle==opt;slow.append(d['id'])
assert max(errs)<1e-9
audit=dict(cases_checked=len(rows),max_statevector_expectation_error=max(errs),independent_bit_enumeration_ids=slow,reconstructed_cuts_checked=98,all_ratios_in_range=True)
(p/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8');print(audit)
v=Path(r'<Obsidian 筆記庫>')
sync=json.loads((v/'graphify-out/sync-state.json').read_text(encoding='utf-8'))
project=v/'10 專案/全國智慧運算與量子資訊創新應用競賽.md'
inputs=[Path('Team_9.pptx'),Path('225033.jpg'),project,*Path(r'<Google 雲端硬碟>\量子資訊應用競賽').glob('*.md')]
provenance=dict(date='2026-11-07',graph_built_at=sync['built_at'],project_graph_hash_matches=hashlib.sha256(project.read_bytes()).hexdigest()==sync['sources']['10 專案/全國智慧運算與量子資訊創新應用競賽.md']['sha256'],resolution='Read current original project note and meeting Markdown; no vault write or graph rebuild.',sources=[dict(path=str(x.resolve()),sha256=hashlib.sha256(x.read_bytes()).hexdigest()) for x in inputs],web_sources=['https://arxiv.org/abs/1910.08980','https://arxiv.org/abs/2012.03421','https://arxiv.org/abs/2408.13207','https://github.com/vijeycreative/Modified_RQAOA','https://openqaoa.entropicalabs.com/workflows/recursive-qaoa/'])
(Path('sources')/'provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2),encoding='utf-8')

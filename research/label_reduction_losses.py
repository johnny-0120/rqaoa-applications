"""Offline reduction-loss labels for existing RQAOA traces.

This is an analysis/oracle tool for small graphs, NOT a deployable solver.
It enumerates original assignments and enforces the recorded constraints.
Labels must never be passed into the online choice rule at evaluation time.
No benchmark results or original solver files are modified.
"""
import csv
import json
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from rqaoa_experiment import cuts,spins

def main():
    traces=json.loads((ROOT/'results/traces.json').read_text(encoding='utf-8'))
    rows=[]; graph_summary=[]
    for case in traces:
        n=case['n']; w=np.zeros((n,n))
        for i,j in case['edges']:w[i,j]=w[j,i]=1.
        z=spins(n);score=cuts(w,z);alive=np.ones(len(z),dtype=bool)
        root_best=float(score.max());total_loss=0.
        for step,record in enumerate(case['rqaoa']['trace'],1):
            i,j,s=record['keep'],record['remove'],record['sign']
            same=z[:,i]==z[:,j]
            best_same=float(score[alive & same].max())
            best_opposite=float(score[alive & ~same].max())
            before=max(best_same,best_opposite)
            chosen=best_same if s==1 else best_opposite
            loss=before-chosen;total_loss+=loss
            rows.append(dict(graph_id=case['id'],family=case['family'],original_n=n,seed=case['seed'],edge_probability=case['p'],step=step,remaining=record['remaining'],keep=i,remove=j,chosen_sign=s,abs_correlation=abs(record['correlation']),top_gap=record['gap'],current_optimum=before,best_same=best_same,best_opposite=best_opposite,irreversible_loss=loss,loses_any_score=int(loss>1e-9),loss_over_original_edge_count=loss/max(len(case['edges']),1)))
            alive &= same if s==1 else ~same
        final_optimum=float(score[alive].max())
        assert abs(root_best-final_optimum-total_loss)<1e-8
        assert abs(final_optimum-case['rqaoa']['cut'])<1e-8
        graph_summary.append(dict(graph_id=case['id'],original_optimum=root_best,final_cut=case['rqaoa']['cut'],sum_reduction_loss=total_loss,first_positive_loss_step=next((r['step'] for r in rows if r['graph_id']==case['id'] and r['irreversible_loss']>1e-9),None)))
    folder=ROOT/'research/results';folder.mkdir(parents=True,exist_ok=True)
    with (folder/'reduction_loss_labels.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary=dict(status='offline_oracle_labels_only_not_predictor_validation',graphs=len(graph_summary),labeled_steps=len(rows),positive_loss_steps=sum(r['loses_any_score'] for r in rows),all_telescope_checks_passed=True,all_terminal_scores_match_original=True,positive_loss_records=[r for r in rows if r['irreversible_loss']>1e-9],graphs_with_loss=[r for r in graph_summary if r['sum_reduction_loss']>1e-9])
    (folder/'reduction_loss_audit.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()

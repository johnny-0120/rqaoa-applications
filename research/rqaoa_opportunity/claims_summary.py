"""Summarise claims_*.csv (step0 and lock)."""
import csv,glob,collections as C,numpy as np
from sweep import OUT

def f(x):
    try:return float(x)
    except:return None

for path in sorted(glob.glob(str(OUT/'claims_step0_*.csv'))):
    rows=list(csv.DictReader(open(path,encoding='utf-8-sig')));g=C.defaultdict(list)
    for r in rows:g[r['signal']].append(r)
    print('\n==',path.split('claims_step0_')[1][:-4],f'({len({r["seed"] for r in rows})} graphs)')
    print(f'  {"signal":10s} {"copyJ":>6s} {"copy100%":>8s} {"signacc":>7s} {"frust_hit":>9s} {"top_ok":>6s} {"sec":>6s}  gamma*/beta (median)')
    for s in ['naive','local_r1','local_r2','local_r3','sdp','qaoa_p1','qaoa_p2','qaoa_p3']:
        v=g.get(s)
        if not v:continue
        m=lambda k:np.mean([f(r[k]) for r in v if f(r[k]) is not None])
        ang=''
        if s.startswith('qaoa'):
            G=np.array([[float(a) for a in r['gamma_scaled'].split(';')] for r in v]);B=np.array([[float(a) for a in r['beta'].split(';')] for r in v])
            ang=' '.join(f'{a:.2f}' for a in np.median(G,0))+' / '+' '.join(f'{a:.2f}' for a in np.median(B,0))
        print(f'  {s:10s} {m("copyJ"):6.3f} {np.mean([f(r["copyJ"])==1 for r in v]):8.2f} {m("signacc"):7.3f} {m("frust_hit"):9.3f} {m("top_ok"):6.2f} {m("seconds"):6.1f}  {ang}')

for path in sorted(glob.glob(str(OUT/'claims_lock_*.csv'))):
    rows=list(csv.DictReader(open(path,encoding='utf-8-sig')));g=C.defaultdict(list)
    for r in rows:g[r['signal']].append(r)
    print('\n== LOCK',path.split('claims_lock_')[1][:-4])
    for s in ['random','naive','local_r1','local_r2','local_r3','sdp','qaoa_p1','qaoa_p2','qaoa_p3']:
        v=g.get(s)
        if not v:continue
        print(f'  {s:10s} N={len(v):3d} locked={sum(int(r["locked"]) for r in v):3d} meangap={np.mean([float(r["gap"]) for r in v]):.3f} sig_s={np.mean([float(r["signal_seconds"]) for r in v]):7.1f}')

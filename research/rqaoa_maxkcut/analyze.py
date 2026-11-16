import csv, json, sys
from collections import defaultdict
import numpy as np
out = {}
for s in sys.argv[1:]:
    rows = list(csv.DictReader(open(f"results/{s}.csv", encoding="utf-8")))
    c = defaultdict(list); t = defaultdict(list)
    for r in rows:
        key = (int(r["n"]), r["param"], int(r["seed"]))
        c[key + (r["method"],)].append(float(r["cost"])); t[key + (r["method"],)].append(float(r["time_s"]))
    inst = sorted({k[:3] for k in c if k[:3] + ("best_known",) in c})
    methods = sorted({k[3] for k in c})
    cost = {(i, m): np.mean(c[i + (m,)]) for i in inst for m in methods if i + (m,) in c}
    groups = sorted({i[:2] for i in inst})
    res = {}
    for g in groups:
        ii = [i for i in inst if i[:2] == g]
        gr = {"instances": len(ii), "methods": {}}
        for m in methods:
            v = [cost[(i, m)] for i in ii if (i, m) in cost]
            gaps = [cost[(i, m)] - cost[(i, "best_known")] for i in ii if (i, m) in cost]
            tt = [np.mean(t[i + (m,)]) for i in ii if i + (m,) in t]
            gr["methods"][m] = dict(mean=round(float(np.mean(v)), 4), mean_gap=round(float(np.mean(gaps)), 4),
                                    hits=int(sum(abs(x) < 1e-6 for x in gaps)), median_time=round(float(np.median(tt)), 4))
        wtl = {}
        for base in ("RQAOA1", "RQAOA1+descent"):
            for m in methods:
                if m in (base, "best_known", "SDP_lowerbound"):
                    continue
                a = np.array([cost[(i, base)] for i in ii]); b = np.array([cost[(i, m)] for i in ii])
                wtl[f"{base} vs {m}"] = [int((a < b - 1e-6).sum()), int((abs(a - b) <= 1e-6).sum()), int((a > b + 1e-6).sum())]
        gr["win_tie_loss"] = wtl
        res[f"n={g[0]},{g[1]}"] = gr
    out[s] = res
print(json.dumps(out, indent=1, ensure_ascii=False))
json.dump(out, open("results/analysis.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)

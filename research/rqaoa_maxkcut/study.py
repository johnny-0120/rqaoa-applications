"""Main study: qutrit RQAOA for 3-channel Wi-Fi planning + planted 3-colouring reproduction.

usage: python study.py <set>   set in {planted, wifi}
Rows are appended to results/<set>.csv (one row per instance x method). Seeds for the confirmation
sets were fixed before any result was seen; nothing is tuned on them.
"""
import csv, json, os, sys, time
import numpy as np
from rqk import maxkcut_problem, rqaoa, Sim, optimize_qaoa1
from baselines import (bravyi_graph, wifi_instance, cocolour, greedy, tabu, descent,
                       sdp_vectors, fj_round)

k = 3
SET = sys.argv[1]
OUT = f"results/{SET}.csv"
TRACE = f"results/{SET}_traces.jsonl"
FIELDS = ["set", "n", "param", "seed", "m", "total_w", "method", "cost", "time_s", "extra"]

if SET == "planted":
    CASES = [(n, d, s) for n in (30, 60) for d in (4, 6) for s in range(5)]
elif SET == "wifi":
    CASES = [(n, None, s) for n in (24, 48, 96) for s in range(100, 110)]
else:
    raise SystemExit("unknown set")

done = set()
if os.path.exists(OUT):
    for r in csv.DictReader(open(OUT, encoding="utf-8")):
        done.add((int(r["n"]), r["param"], int(r["seed"])))
f = open(OUT, "a", newline="", encoding="utf-8")
w = csv.DictWriter(f, FIELDS)
if not done:
    w.writeheader()


def row(n, param, seed, m, tw, method, cost, t, extra=""):
    w.writerow(dict(set=SET, n=n, param=param, seed=seed, m=m, total_w=round(tw, 6), method=method,
                    cost=round(float(cost), 9), time_s=round(t, 4), extra=extra))
    f.flush()


for n, d, seed in CASES:
    param = str(d) if d else "wifi"
    if (n, param, seed) in done:
        continue
    rng = np.random.default_rng(seed)
    if SET == "planted":
        E, W = bravyi_graph(n, d, rng)
    else:
        E, W, pts = wifi_instance(n, rng)
    m, tw = len(E), float(W.sum())
    pr = maxkcut_problem(n, E, W, k)

    # QAOA level 1 (expected co-channel weight of the optimised state)
    t = time.perf_counter()
    e1, g1, b1 = optimize_qaoa1(Sim(pr), np.pi)
    row(n, param, seed, m, tw, "QAOA1_expected", e1, time.perf_counter() - t, f"gamma={g1:.4f}")

    # RQAOA level 1, re-optimised every step, n_c = 6 brute force
    t = time.perf_counter()
    x, c, tr = rqaoa(pr, nc=6)
    trq = time.perf_counter() - t
    assert abs(cocolour(x, E, W) - c) < 1e-6
    row(n, param, seed, m, tw, "RQAOA1", c, trq)
    xd = descent(n, E, W, k, x, rng)
    row(n, param, seed, m, tw, "RQAOA1+descent", cocolour(xd, E, W), trq)
    with open(TRACE, "a", encoding="utf-8") as tf:
        tf.write(json.dumps(dict(n=n, param=param, seed=seed, x=[int(v) for v in x], trace=tr)) + "\n")

    # practical heuristics
    t = time.perf_counter(); xg = greedy(n, E, W, k); tg = time.perf_counter() - t
    row(n, param, seed, m, tw, "DSatur_greedy", cocolour(xg, E, W), tg)
    row(n, param, seed, m, tw, "DSatur+descent", cocolour(descent(n, E, W, k, xg, rng), E, W), tg)
    # LCCS-style: every AP switches to its least-interfered channel until stable, from all-on-one-channel
    t = time.perf_counter(); xl = descent(n, E, W, k, np.zeros(n, int), rng)
    row(n, param, seed, m, tw, "LCCS_best_response", cocolour(xl, E, W), time.perf_counter() - t)

    # SDP + Frieze-Jerrum rounding (100 samples)
    t = time.perf_counter()
    vec, lb = sdp_vectors(n, E, W, k)
    S = fj_round(vec, k, rng, 100)
    cs = np.array([cocolour(s, E, W) for s in S])
    tsdp = time.perf_counter() - t
    row(n, param, seed, m, tw, "SDP_lowerbound", lb, tsdp)
    row(n, param, seed, m, tw, "SDP_FJ_mean", cs.mean(), tsdp)
    row(n, param, seed, m, tw, "SDP_FJ_best100", cs.min(), tsdp)
    cd = [cocolour(descent(n, E, W, k, s, rng), E, W) for s in S]
    row(n, param, seed, m, tw, "SDP_FJ+descent_best100", min(cd), time.perf_counter() - t)

    # multi-start best-response and tabu with RQAOA's CPU time (3 seeds each, averaged later)
    for s in range(3):
        r2 = np.random.default_rng(1000 + s)
        t = time.perf_counter(); best = np.inf; starts = 0
        while time.perf_counter() - t < trq:
            best = min(best, cocolour(descent(n, E, W, k, r2.integers(k, size=n), r2), E, W)); starts += 1
        row(n, param, seed, m, tw, "multistart_descent_sameCPU", best, trq, f"rep={s};starts={starts}")
        for frac, name in ((1.0, "tabu_sameCPU"), (0.01, "tabu_1pctCPU")):
            xb, cb, it, tt = tabu(n, E, W, k, time_limit=trq * frac, rng=np.random.default_rng(2000 + s))
            row(n, param, seed, m, tw, name, cb, tt, f"rep={s};iters={it}")
    # best known: long tabu (10 s x 2) unless planted (optimum 0 known)
    if SET == "wifi":
        bk = min(tabu(n, E, W, k, time_limit=10, rng=np.random.default_rng(3000 + s))[1] for s in range(2))
        row(n, param, seed, m, tw, "best_known", min(bk, c), 20.0)
    else:
        row(n, param, seed, m, tw, "best_known", 0.0, 0.0, "planted optimum")
    print(SET, n, param, seed, "RQ", round(c, 3), "time", round(trq, 1), flush=True)
f.close()

import sys, time, json
import numpy as np
from rqk import maxkcut_problem, rqaoa, qaoa1_expectation
from baselines import bravyi_graph, wifi_instance, cocolour, greedy, tabu, sdp_vectors, fj_round, descent

k = 3
kind = sys.argv[1]
n = int(sys.argv[2])
seed = int(sys.argv[3])
rng = np.random.default_rng(seed)
if kind == "bravyi":
    d = int(sys.argv[4])
    E, W = bravyi_graph(n, d, rng)
else:
    E, W, _ = wifi_instance(n, rng)
print("n", n, "m", len(E), "W", W.sum().round(2))
pr = maxkcut_problem(n, E, W, k)
t = time.perf_counter()
e, g, b = qaoa1_expectation(pr)
print("QAOA1 expected co-colour", round(e, 3), "gamma", round(g, 3), time.perf_counter() - t)
t = time.perf_counter()
x, cost, tr = rqaoa(pr, nc=6)
trq = time.perf_counter() - t
print("RQAOA", cost, "time", round(trq, 2))
xg = greedy(n, E, W, k)
print("greedy", cocolour(xg, E, W), "greedy+descent", cocolour(descent(n, E, W, k, xg, rng), E, W))
t = time.perf_counter()
vec, lb = sdp_vectors(n, E, W, k)
S = fj_round(vec, k, rng, 100)
c = np.array([cocolour(s, E, W) for s in S])
print("SDP lb", round(lb, 3), "FJ mean", c.mean().round(3), "best", c.min(), "time", round(time.perf_counter() - t, 2))
cd = [cocolour(descent(n, E, W, k, s, rng), E, W) for s in S[:20]]
print("FJ+descent best of 20", min(cd))
for frac in (0.01, 0.1, 1.0):
    res = [tabu(n, E, W, k, time_limit=trq * frac, rng=np.random.default_rng(s))[1:3] for s in range(3)]
    print("tabu", frac, res)

"""Check the light-cone p=1 simulator against a dense qudit statevector, and elimination against brute force."""
import itertools, json, sys
import numpy as np
from rqk import Problem, Sim, fourier, brute

rng = np.random.default_rng(7)
k = 3
report = {}


def random_problem(n, p_edge=0.6):
    pr = Problem(k)
    pr.nodes = set(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < p_edge:
                pr.add(i, j, rng.normal(size=k))
    return pr


def dense_state(pr, gamma, beta):
    n = len(pr.nodes)
    X = np.array(list(itertools.product(range(k), repeat=n)))
    E = np.zeros(len(X))
    for (i, j), f in pr.F.items():
        E += f[(X[:, j] - X[:, i]) % k]
    psi = np.exp(-1j * gamma * E) / np.sqrt(k ** n)
    U = fourier(k)
    B = U @ np.diag(np.exp(1j * beta)) @ U.conj().T
    psi = psi.reshape([k] * n)
    for q in range(n):
        psi = np.moveaxis(np.tensordot(B, psi, axes=([1], [q])), 0, q)
    return psi.reshape(-1), X, E


maxerr_M = maxerr_E = 0.0
for trial in range(6):
    n = 6
    pr = random_problem(n)
    pr.const = 0.3
    g = rng.uniform(0, 2)
    beta = np.r_[0, rng.uniform(0, 2 * np.pi, 2)]
    psi, X, E = dense_state(pr, g, beta)
    prob = np.abs(psi) ** 2
    sim = Sim(pr)
    sig = sim.sigmas(g)
    e_sim = sim.energy_from_T(sim.energy_tensor(sig), beta)
    e_dense = float(prob @ E) + pr.const
    maxerr_E = max(maxerr_E, abs(e_sim - e_dense))
    M = sim.marginals(sig, beta)
    for pidx, (i, j) in enumerate(sim.pairs):
        for b in range(k):
            md = prob[(X[:, j] - X[:, i]) % k == b].sum()
            maxerr_M = max(maxerr_M, abs(md - M[pidx, b]))
report["max_err_energy"] = maxerr_E
report["max_err_pair_marginal"] = maxerr_M

# elimination: min over constrained space equals min of reduced problem
maxerr_elim = 0.0
for trial in range(20):
    n = 6
    pr = random_problem(n, 0.7)
    i, j = rng.choice(n, 2, replace=False)
    b = int(rng.integers(k))
    X = np.array(list(itertools.product(range(k), repeat=n)))
    E = np.array([pr.energy(x) for x in X])
    mask = (X[:, j] - X[:, i]) % k == b
    red = pr.copy()
    red.eliminate(int(i), int(j), b)
    e_red, sol = brute(red)
    maxerr_elim = max(maxerr_elim, abs(E[mask].min() - e_red))
    # also check that reconstruction reproduces energy
    x = np.zeros(n, int)
    for v, c in sol.items():
        x[v] = c
    x[i] = (x[j] - b) % k
    maxerr_elim = max(maxerr_elim, abs(pr.energy(x) - e_red))
report["max_err_elimination"] = maxerr_elim
print(json.dumps(report, indent=1))
ok = all(v < 1e-9 for v in report.values())
json.dump(dict(report, ok=ok), open("results/validation.json", "w"), indent=1)
sys.exit(0 if ok else 1)

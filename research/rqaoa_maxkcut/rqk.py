"""Level-1 RQAOA for MAX-k-CUT / weighted k-coloring (qudit version, Bravyi et al. 2022).

Cost (to MINIMIZE):  H(x) = const + sum_{i<j} F_ij[(x_j - x_i) mod k]
MAX-k-CUT conflict cost: F_ij = w_ij * [1, 0, ..., 0].

Exact p=1 classical simulation: two-qudit reduced density matrices only need the
neighbourhood of the pair (light cone of a single layer).
"""
import itertools
import numpy as np


class Problem:
    def __init__(self, k):
        self.k = k
        self.F = {}          # (i,j) with i<j -> np.array(k), function of (x_j - x_i)
        self.const = 0.0
        self.nodes = set()

    def copy(self):
        p = Problem(self.k)
        p.F = {e: v.copy() for e, v in self.F.items()}
        p.const = self.const
        p.nodes = set(self.nodes)
        return p

    def add(self, i, j, f):
        """add term f[(x_j - x_i) mod k]"""
        if i == j:
            self.const += f[0]
            return
        f = np.asarray(f, float)
        if i > j:
            i, j = j, i
            f = f[(-np.arange(self.k)) % self.k]
        if (i, j) in self.F:
            self.F[(i, j)] = self.F[(i, j)] + f
        else:
            self.F[(i, j)] = f.copy()

    def energy(self, x):
        e = self.const
        for (i, j), f in self.F.items():
            e += f[(x[j] - x[i]) % self.k]
        return e

    def neighbors(self):
        nb = {v: [] for v in self.nodes}
        for (i, j) in self.F:
            nb[i].append(j)
            nb[j].append(i)
        return nb

    def eliminate(self, i, j, b):
        """impose x_j = x_i + b (mod k), remove x_i (keep j)."""
        k = self.k
        newF = {}
        const = self.const
        for (a, c), f in self.F.items():
            if {a, c} == {i, j}:
                # f[(x_c - x_a)]; with x_j - x_i = b
                d = b if (a, c) == (i, j) else (-b) % k
                const += f[d]
            elif a == i or c == i:
                h = c if a == i else a
                # term as function of (x_h - x_i)
                g = f if a == i else f[(-np.arange(k)) % k]
                # x_i = x_j - b  ->  x_h - x_i = (x_h - x_j) + b
                g2 = g[(np.arange(k) + b) % k]     # function of (x_h - x_j)
                key = (j, h) if j < h else (h, j)
                gg = g2 if j < h else g2[(-np.arange(k)) % k]
                newF[key] = newF.get(key, 0) + gg
            else:
                newF[(a, c)] = newF.get((a, c), 0) + f
        self.F = {e: np.asarray(v, float) for e, v in newF.items()}
        self.const = const
        self.nodes.discard(i)


def maxkcut_problem(n, edges, weights, k):
    p = Problem(k)
    p.nodes = set(range(n))
    for (u, v), w in zip(edges, weights):
        f = np.zeros(k)
        f[0] = w
        p.add(u, v, f)
    return p


# ---------------------------------------------------------------- p=1 simulation
def fourier(k):
    om = np.exp(2j * np.pi / k)
    return np.array([[om ** (a * c) for a in range(k)] for c in range(k)]) / np.sqrt(k)  # col a = |phi_a>


class Sim:
    """Precomputes structure for a fixed Problem."""

    def __init__(self, prob):
        self.p = prob
        k = prob.k
        self.k = k
        self.pairs = list(prob.F.keys())
        self.Fp = np.array([prob.F[e] for e in self.pairs])  # (m,k)
        nb = prob.neighbors()
        termidx = {e: t for t, e in enumerate(self.pairs)}
        # incidences: for pair p=(i,j) and neighbour h (h != i,j): term index for (i,h),(j,h) or -1
        inc_p, inc_ti, inc_tj, inc_si, inc_sj = [], [], [], [], []
        for pidx, (i, j) in enumerate(self.pairs):
            hs = (set(nb[i]) | set(nb[j])) - {i, j}
            for h in hs:
                inc_p.append(pidx)
                for (u, lst_t, lst_s) in ((i, inc_ti, inc_si), (j, inc_tj, inc_sj)):
                    key = (u, h) if u < h else (h, u)
                    if key in termidx:
                        lst_t.append(termidx[key])
                        lst_s.append(1 if u < h else -1)   # +1: term is f(x_h - x_u)
                    else:
                        lst_t.append(-1)
                        lst_s.append(1)
        self.inc_p = np.array(inc_p, int)
        self.inc_ti = np.array(inc_ti, int)
        self.inc_tj = np.array(inc_tj, int)
        self.inc_si = np.array(inc_si, int)
        self.inc_sj = np.array(inc_sj, int)
        order = np.argsort(self.inc_p, kind="stable")
        for name in ("inc_p", "inc_ti", "inc_tj", "inc_si", "inc_sj"):
            setattr(self, name, getattr(self, name)[order])
        self.U = fourier(k)
        d = np.arange(k)
        # D[p, (xi,xj)] = F_p[(xj - xi)]
        self.diff = (d[None, :] - d[:, None]) % k            # [xi,xj] -> xj - xi
        self.D = self.Fp[:, self.diff.reshape(-1)]           # (m, k*k)

    def sigmas(self, gamma):
        """reduced 2-site density matrices (before mixer) for all pairs: (m, k,k, k,k) [xi,xj,xi',xj']"""
        k = self.k
        m = len(self.pairs)
        d = np.arange(k)
        # phase for term t as function of (xh - xu) with orientation sign
        ph = np.exp(-1j * gamma * self.Fp)                    # (m,k)
        sig = np.ones((m, k, k, k, k), complex)
        # pair's own term
        own = ph[:, self.diff]                                # (m, xi, xj)
        sig *= own[:, :, :, None, None] * np.conj(own)[:, None, None, :, :]
        if len(self.inc_p):
            # E_u[inc, xh, xu] = exp(-i g F[(xh - xu)]) or with reversed orientation
            def emat(t, s):
                diffhu = (d[:, None] - d[None, :]) % k        # [xh, xu] -> xh - xu
                idx = np.where(s[:, None, None] > 0, diffhu[None], (-diffhu[None]) % k)
                out = np.ones((len(t), k, k), complex)
                mask = t >= 0
                out[mask] = ph[t[mask]][np.arange(mask.sum())[:, None, None], idx[mask]]
                return out
            Ei = emat(self.inc_ti, self.inc_si)
            Ej = emat(self.inc_tj, self.inc_sj)
            # A[inc, xi, xj, xi', xj'] = 1/k sum_h Ei[h,xi] Ej[h,xj] conj(Ei[h,xi']) conj(Ej[h,xj'])
            A = np.einsum("nha,nhb,nhc,nhd->nabcd", Ei, Ej, np.conj(Ei), np.conj(Ej), optimize=True) / k
            starts = np.r_[0, np.flatnonzero(np.diff(self.inc_p)) + 1]
            prod = np.multiply.reduceat(A, starts, axis=0)
            sig[self.inc_p[starts]] *= prod
        return sig.reshape(m, k * k, k * k) / (k * k)

    def mixer2(self, beta):
        k = self.k
        B = self.U @ np.diag(np.exp(1j * np.asarray(beta))) @ self.U.conj().T
        return np.kron(B, B)

    def energy_tensor(self, sig):
        # T[r,s,t] = sum_p D[p,r] sig[p,s,t]
        return np.einsum("pr,pst->rst", self.D, sig, optimize=True)

    def energy_from_T(self, T, beta):
        BB = self.mixer2(beta)
        return float(np.real(np.einsum("rs,rst,rt->", BB, T, np.conj(BB), optimize=True))) + self.p.const

    def marginals(self, sig, beta):
        BB = self.mixer2(beta)
        rho = np.einsum("rs,pst,ut->pru", BB, sig, np.conj(BB), optimize=True)
        diagp = np.real(np.einsum("prr->pr", rho))            # (m, k*k)  prob of (xi,xj)
        k = self.k
        M = np.zeros((len(self.pairs), k))
        flat = self.diff.reshape(-1)
        for b in range(k):
            M[:, b] = diagp[:, flat == b].sum(1)
        return M


def optimize_beta(sim, T, rng, grid=10):
    k = sim.k
    best = (np.inf, None)
    g = np.linspace(0, 2 * np.pi, grid, endpoint=False)
    for bs in itertools.product(g, repeat=k - 1):
        beta = np.r_[0.0, bs]
        e = sim.energy_from_T(T, beta)
        if e < best[0]:
            best = (e, beta)
    from scipy.optimize import minimize
    r = minimize(lambda b: sim.energy_from_T(T, np.r_[0.0, b]), best[1][1:], method="Nelder-Mead",
                 options={"xatol": 1e-4, "fatol": 1e-9, "maxiter": 400})
    if r.fun < best[0]:
        best = (r.fun, np.r_[0.0, r.x])
    return best


def optimize_qaoa1(sim, gmax, ngrid=24, nrefine=12, rng=None):
    gs = np.linspace(gmax / ngrid, gmax, ngrid)
    res = []
    for g in gs:
        T = sim.energy_tensor(sim.sigmas(g))
        e, beta = optimize_beta(sim, T, rng, grid=8)
        res.append((e, g, beta))
    e0, g0, b0 = min(res, key=lambda t: t[0])
    step = gs[1] - gs[0]
    lo, hi = g0 - step, g0 + step
    from scipy.optimize import minimize_scalar
    cache = {}

    def f(g):
        T = sim.energy_tensor(sim.sigmas(g))
        e, beta = optimize_beta(sim, T, rng, grid=6)
        cache[g] = (e, beta)
        return e
    r = minimize_scalar(f, bounds=(max(lo, 1e-6), hi), method="bounded", options={"maxiter": nrefine, "xatol": 1e-4})
    if r.fun < e0:
        return r.fun, r.x, cache[r.x][1]
    return e0, g0, b0


def brute(prob):
    nodes = sorted(prob.nodes)
    k = prob.k
    best = (np.inf, None)
    idx = {v: t for t, v in enumerate(nodes)}
    pairs = list(prob.F.items())
    if not nodes:
        return prob.const, {}
    X = np.array(list(itertools.product(range(k), repeat=len(nodes))))
    E = np.full(len(X), prob.const)
    for (i, j), f in pairs:
        E += f[(X[:, idx[j]] - X[:, idx[i]]) % k]
    t = int(np.argmin(E))
    return float(E[t]), {v: int(X[t, idx[v]]) for v in nodes}


def rqaoa(prob, nc=6, gmax=np.pi, ngrid=24, nrefine=12, transfer=None, log=None):
    """returns assignment dict, final energy on original problem, trace"""
    orig = prob
    prob = prob.copy()
    k = prob.k
    constraints = []
    trace = []
    while len(prob.nodes) > nc:
        if not prob.F:
            break
        sim = Sim(prob)
        if transfer is not None:
            g, beta = transfer(prob)
            T = sim.energy_tensor(sim.sigmas(g))
            e = sim.energy_from_T(T, beta)
        else:
            e, g, beta = optimize_qaoa1(sim, gmax, ngrid, nrefine)
        M = sim.marginals(sim.sigmas(g), beta)
        pidx, b = np.unravel_index(np.argmax(M), M.shape)
        i, j = sim.pairs[pidx]
        # eliminate the endpoint with smaller degree? keep deterministic: eliminate i
        trace.append(dict(n=len(prob.nodes), m=len(sim.pairs), gamma=float(g), beta=[float(x) for x in beta],
                          E=float(e), pair=(int(i), int(j)), b=int(b), M=float(M[pidx, b])))
        constraints.append((i, j, int(b)))
        prob.eliminate(i, j, int(b))
        if log:
            log(trace[-1])
    # isolated nodes (no terms) get color 0
    _, sol = brute(prob) if prob.nodes else (prob.const, {})
    x = dict(sol)
    for i, j, b in reversed(constraints):
        x[i] = (x[j] - b) % k
    for v in orig.nodes:
        x.setdefault(v, 0)
    arr = np.array([x[v] for v in range(max(orig.nodes) + 1)])
    return arr, orig.energy(arr), trace


def qaoa1_expectation(prob, gmax=np.pi):
    sim = Sim(prob)
    e, g, beta = optimize_qaoa1(sim, gmax)
    return e, g, beta


def rms_weight(prob):
    if not prob.F:
        return 1.0
    v = np.array([f.max() - f.min() for f in prob.F.values()])
    return float(np.sqrt(np.mean(v ** 2))) or 1.0


def make_transfer(prob, gmax=np.pi):
    """optimise QAOA1 once on the full problem; later steps reuse beta and rescale gamma by RMS term range."""
    sim = Sim(prob)
    e, g, beta = optimize_qaoa1(sim, gmax)
    a = g * rms_weight(prob)

    def rule(p):
        return a / rms_weight(p), beta
    return rule, (e, g, beta)

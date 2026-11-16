"""Instance generators and classical baselines for weighted MAX-k-CUT (minimise co-colour weight)."""
import time
import numpy as np


# ------------------------------------------------------------------ instances
def bravyi_graph(n, d, rng, max_tries=1000, reject_triangles=False):
    """Planted 3-colourable graph, close to G[d,n] of Bravyi et al. 2022: union of three random
    (d/2)-regular bipartite graphs between planted colour classes, total degree d (d even);
    reject if disconnected. (The paper also rejects triangles; at these sizes that almost never
    succeeds, so it is optional here and off by default.)"""
    assert n % 3 == 0 and d % 2 == 0
    d = d // 2
    for _ in range(max_tries):
        perm = rng.permutation(n)
        V = [perm[:n // 3], perm[n // 3:2 * n // 3], perm[2 * n // 3:]]
        E = set()
        ok = True
        for r in range(3):
            for s in range(r + 1, 3):
                for _attempt in range(200):
                    deg = {int(w): 0 for w in V[s]}
                    sub = []
                    fail = False
                    for v in V[r]:
                        cand = [w for w in deg if deg[w] < d]
                        if len(cand) < d:
                            fail = True
                            break
                        ws = rng.choice(cand, d, replace=False)
                        for w in ws:
                            deg[int(w)] += 1
                            sub.append((int(min(v, w)), int(max(v, w))))
                    if not fail:
                        break
                else:
                    ok = False
                E.update(sub)
        if not ok:
            continue
        E = sorted(E)
        adj = [set() for _ in range(n)]
        for u, v in E:
            adj[u].add(v); adj[v].add(u)
        if reject_triangles and any(adj[u] & adj[v] for u, v in E):
            continue
        seen = {0}; st = [0]
        while st:
            u = st.pop()
            for w in adj[u]:
                if w not in seen:
                    seen.add(w); st.append(w)
        if len(seen) < n:
            continue
        return E, np.ones(len(E))
    raise RuntimeError("generation failed")


def wifi_instance(n_ap, rng, floors=None, area=(60.0, 40.0), floor_att=15.0,
                  ptx=20.0, pl_exp=3.5, wall_db_per_m=0.4, thresh_dbm=-82.0):
    """Indoor 2.4 GHz AP channel-assignment instance (channels 1/6/11 -> k=3).
    APs placed on a jittered grid over `floors`; interference weight = received power (mW, linear)
    of AP j at AP i above the CCA threshold, symmetrised (mean). Log-distance path loss, 40 dB @1 m,
    plus distance-proportional wall loss and per-floor attenuation. Weights normalised to mean 1."""
    if floors is None:
        floors = max(1, int(round(n_ap / 20)))
    per = int(np.ceil(n_ap / floors))
    pts = []
    for f in range(floors):
        nx = int(np.ceil(np.sqrt(per * area[0] / area[1])))
        ny = int(np.ceil(per / nx))
        xs = (np.arange(nx) + 0.5) * area[0] / nx
        ys = (np.arange(ny) + 0.5) * area[1] / ny
        grid = [(x, y) for x in xs for y in ys]
        sel = rng.choice(len(grid), per, replace=False)
        for s in sel:
            x, y = grid[s]
            pts.append((x + rng.normal(0, 0.18 * area[0] / nx), y + rng.normal(0, 0.18 * area[1] / ny), f))
    pts = np.array(pts[:n_ap])
    n = len(pts)
    E, W = [], []
    for i in range(n):
        for j in range(i + 1, n):
            dxy = np.hypot(pts[i, 0] - pts[j, 0], pts[i, 1] - pts[j, 1])
            df = abs(pts[i, 2] - pts[j, 2])
            dist = max(1.0, np.sqrt(dxy ** 2 + (3.5 * df) ** 2))
            pl = 40.0 + 10 * pl_exp * np.log10(dist) + wall_db_per_m * dxy + floor_att * df
            shadow = rng.normal(0, 4.0)
            rx = ptx - pl - shadow
            if rx > thresh_dbm:
                E.append((i, j)); W.append(10 ** (rx / 10))
    W = np.array(W)
    W = W / W.mean()
    return E, W, pts


# ------------------------------------------------------------------ evaluation
def cocolour(x, E, W):
    x = np.asarray(x)
    a = np.array([e[0] for e in E]); b = np.array([e[1] for e in E])
    return float(W[x[a] == x[b]].sum())


def adjacency(n, E, W):
    A = np.zeros((n, n))
    for (u, v), w in zip(E, W):
        A[u, v] += w; A[v, u] += w
    return A


# ------------------------------------------------------------------ greedy (DSatur-like, weighted)
def greedy(n, E, W, k, rng=None):
    A = adjacency(n, E, W)
    x = -np.ones(n, int)
    strength = A.sum(1)
    colored_w = np.zeros((n, k))     # weight to already-coloured neighbours of each colour
    for _ in range(n):
        unc = np.flatnonzero(x < 0)
        # saturation: number of distinct neighbour colours, tie-break by strength
        sat = (colored_w[unc] > 0).sum(1)
        v = unc[np.lexsort((-strength[unc], -sat))[0]]
        c = int(np.argmin(colored_w[v]))
        x[v] = c
        colored_w[:, c] += A[:, v]
    return x


# ------------------------------------------------------------------ local descent + tabu
def tabu(n, E, W, k, time_limit=None, max_iters=None, rng=None, x0=None, tenure_a=10, tenure_b=0.6):
    """Weighted Tabucol (Hertz & de Werra 1987; tenure rule of Galinier & Hao 1999). Moves: recolour
    a vertex currently in conflict. Aspiration on new best. Stops on time/iteration budget or 0 cost."""
    rng = rng or np.random.default_rng()
    A = adjacency(n, E, W)
    x = rng.integers(k, size=n) if x0 is None else np.array(x0)
    G = np.stack([A[:, x == c].sum(1) for c in range(k)], 1)   # G[v,c] weight of nbrs with colour c
    cost = 0.5 * G[np.arange(n), x].sum()
    best, bestx = cost, x.copy()
    tabu_until = np.zeros((n, k))
    it = 0
    t0 = time.perf_counter()
    while True:
        if best <= 1e-12:
            break
        if max_iters is not None and it >= max_iters:
            break
        if time_limit is not None and (it & 31) == 0 and time.perf_counter() - t0 > time_limit:
            break
        it += 1
        cur = G[np.arange(n), x]
        conf = np.flatnonzero(cur > 1e-12)
        if len(conf) == 0:
            break
        delta = G[conf] - cur[conf, None]
        delta[np.arange(len(conf)), x[conf]] = np.inf
        istabu = tabu_until[conf] > it
        aspir = cost + delta < best - 1e-12
        delta = np.where(istabu & ~aspir, np.inf, delta)
        mn = delta.min()
        if not np.isfinite(mn):
            continue
        cand = np.argwhere(delta <= mn + 1e-12)
        r, c = cand[rng.integers(len(cand))]
        v = conf[r]
        old = x[v]
        x[v] = c
        G[:, old] -= A[:, v]
        G[:, c] += A[:, v]
        cost += mn
        tabu_until[v, old] = it + tenure_a * rng.random() + tenure_b * len(conf)
        if cost < best - 1e-12:
            best, bestx = cost, x.copy()
    return bestx, cocolour(bestx, E, W), it, time.perf_counter() - t0


def descent(n, E, W, k, x, rng):
    """steepest 1-opt descent to a local minimum."""
    A = adjacency(n, E, W)
    x = np.array(x)
    G = np.stack([A[:, x == c].sum(1) for c in range(k)], 1)
    while True:
        cur = G[np.arange(n), x]
        delta = G - cur[:, None]
        v, c = np.unravel_index(np.argmin(delta), delta.shape)
        if delta[v, c] >= -1e-12:
            return x
        old = x[v]; x[v] = c
        G[:, old] -= A[:, v]; G[:, c] += A[:, v]


# ------------------------------------------------------------------ SDP (Frieze-Jerrum) rounding
def sdp_vectors(n, E, W, k):
    import cvxpy as cp
    X = cp.Variable((n, n), PSD=True)
    cons = [cp.diag(X) == 1]
    a = np.array([e[0] for e in E]); b = np.array([e[1] for e in E])
    cons.append(X[a, b] >= -1.0 / (k - 1))
    # minimise co-colour relaxation: maximise sum w (k-1)/k (1 - X_ij)
    obj = cp.Minimize(cp.sum(cp.multiply(W, X[a, b])))
    prob = cp.Problem(obj, cons)
    prob.solve(solver=cp.SCS, eps=1e-6, max_iters=20000)
    Xv = (X.value + X.value.T) / 2
    w, V = np.linalg.eigh(Xv)
    w = np.clip(w, 0, None)
    vecs = V * np.sqrt(w)
    lower = float(np.sum(W * (1 + (k - 1) * Xv[a, b]) / k))   # relaxation value of co-colour weight
    return vecs, lower


def fj_round(vecs, k, rng, samples=100):
    n, r = vecs.shape
    G = rng.normal(size=(samples, k, r))
    return np.argmax(np.einsum("nr,skr->snk", vecs, G), axis=2)   # (samples, n)

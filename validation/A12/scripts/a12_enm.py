import numpy as np

def pair_edges(X, cutoff, gamma=1.0, block=512):
    X = np.asarray(X, float)
    N = len(X)
    ii, jj = ([], [])
    for a in range(0, N, block):
        blk = X[a:a + block]
        d2 = ((blk[:, None, :] - X[None, :, :]) ** 2).sum(-1)
        r, c = np.nonzero(d2 <= cutoff * cutoff)
        r = r + a
        m = c > r
        ii.append(r[m])
        jj.append(c[m])
    i = np.concatenate(ii) if ii else np.zeros(0, int)
    j = np.concatenate(jj) if jj else np.zeros(0, int)
    return (i, j, np.full(len(i), float(gamma)))

def hessian(X, i, j, g):
    X = np.asarray(X, float)
    N = len(X)
    v = X[j] - X[i]
    d2 = (v * v).sum(1)
    blk = -(g / d2)[:, None, None] * v[:, :, None] * v[:, None, :]
    H = np.zeros(9 * N * N)
    a = np.arange(3)

    def flat(r, c):
        return ((3 * r[:, None, None] + a[None, :, None]) * (3 * N) + (3 * c[:, None, None] + a[None, None, :])).ravel()
    np.add.at(H, flat(i, j), blk.ravel())
    np.add.at(H, flat(j, i), blk.ravel())
    np.add.at(H, flat(i, i), (-blk).ravel())
    np.add.at(H, flat(j, j), (-blk).ravel())
    return H.reshape(3 * N, 3 * N)

def anm(X, cutoff=15.0, gamma=1.0):
    X = np.asarray(X, float)
    i, j, g = pair_edges(X, cutoff, gamma)
    return hessian(X, i, j, g)

def modes(H, n_zero=6):
    w, V = np.linalg.eigh(H)
    scale = w[-1]
    if not w[n_zero - 1] < 1e-08 * scale:
        raise ValueError('fewer than six zero modes: not a free (unrestrained) elastic network')
    if not w[n_zero] > 1e-08 * scale:
        raise ValueError('more than six zero modes: network not connected at this cutoff')
    return (w[n_zero:], V[:, n_zero:])

def covariance(w, V, n_modes=None):
    if n_modes is not None:
        w, V = (w[:n_modes], V[:, :n_modes])
    return V / w @ V.T

def msf(G):
    N = G.shape[0] // 3
    return np.einsum('iaib->i', G.reshape(N, 3, N, 3))

def prs_amplitude(G):
    N = G.shape[0] // 3
    G4 = G.reshape(N, 3, N, 3)
    return np.sqrt((G4 ** 2).sum(axis=(1, 3)) / 3.0)

def _masked_mean(ratio, src, axis_nodes):
    S = ratio.shape[0]
    mask = np.ones(ratio.shape, bool)
    mask[np.arange(S), src] = False
    return (ratio * mask).sum(0) / mask.sum(0)

def sensitivity(A, src_idx):
    src = np.asarray(src_idx, int)
    ratio = A[src, :] / np.diag(A)[src][:, None]
    return _masked_mean(ratio, src, None)

def effectiveness(A, tgt_idx):
    tgt = np.asarray(tgt_idx, int)
    ratio = (A[:, tgt] / np.diag(A)[:, None]).T
    return _masked_mean(ratio, tgt, None)

def path_score(A, hub_idx, tgt_idx):
    return np.log10(sensitivity(A, hub_idx)) + np.log10(effectiveness(A, tgt_idx))

def min_distance_to(X, idx_set):
    X = np.asarray(X, float)
    Y = X[np.asarray(idx_set, int)]
    d = np.sqrt(((X[:, None, :] - Y[None, :, :]) ** 2).sum(-1))
    d[np.asarray(idx_set, int), np.arange(len(idx_set))] = np.inf
    out = d.min(axis=1)
    out[np.asarray(idx_set, int)] = np.where(np.isfinite(out[np.asarray(idx_set, int)]), out[np.asarray(idx_set, int)], 0.0)
    return out

def condensation(P, S, cut_ps=10.0, cut_ss=10.0, g_ps=1.0, g_ss=1.0, ridge=1e-09):
    P = np.asarray(P, float)
    S = np.asarray(S, float)
    K = len(S)
    d2 = ((P[:, None, :] - S[None, :, :]) ** 2).sum(-1)
    ip, ks = np.nonzero(d2 <= cut_ps * cut_ps)
    if len(ip) == 0:
        return (np.zeros(0, int), np.zeros((0, 0)))
    M_idx = np.unique(ip)
    mi = np.searchsorted(M_idx, ip)
    nm = len(M_idx)
    v = S[ks] - P[ip]
    dd2 = (v * v).sum(1)
    blk = (g_ps / dd2)[:, None, None] * v[:, :, None] * v[:, None, :]
    D = np.zeros((3 * nm, 3 * nm))
    B = np.zeros((3 * nm, 3 * K))
    Hss = np.zeros((3 * K, 3 * K))
    for a in range(3):
        for b in range(3):
            np.add.at(D, (3 * mi + a, 3 * mi + b), blk[:, a, b])
            np.add.at(B, (3 * mi + a, 3 * ks + b), -blk[:, a, b])
            np.add.at(Hss, (3 * ks + a, 3 * ks + b), blk[:, a, b])
    ds = ((S[:, None, :] - S[None, :, :]) ** 2).sum(-1)
    k1, k2 = np.nonzero(np.triu(ds <= cut_ss * cut_ss, 1))
    if len(k1):
        vs = S[k2] - S[k1]
        bs = (g_ss / (vs * vs).sum(1))[:, None, None] * vs[:, :, None] * vs[:, None, :]
        for a in range(3):
            for b in range(3):
                np.add.at(Hss, (3 * k1 + a, 3 * k1 + b), bs[:, a, b])
                np.add.at(Hss, (3 * k2 + a, 3 * k2 + b), bs[:, a, b])
                np.add.at(Hss, (3 * k1 + a, 3 * k2 + b), -bs[:, a, b])
                np.add.at(Hss, (3 * k2 + a, 3 * k1 + b), -bs[:, a, b])
    Hss = Hss + ridge * np.eye(3 * K)
    dH = D - B @ np.linalg.solve(Hss, B.T)
    dH = 0.5 * (dH + dH.T)
    return (M_idx, dH)

def msf_after(G, M_idx, dH, nodes):
    nodes = np.asarray(nodes, int)
    rows = (3 * nodes[:, None] + np.arange(3)[None, :]).ravel()
    if len(M_idx) == 0:
        N = G.shape[0] // 3
        return np.einsum('iaib->i', G[np.ix_(rows, rows)].reshape(len(nodes), 3, len(nodes), 3))
    dof = (3 * np.asarray(M_idx)[:, None] + np.arange(3)[None, :]).ravel()
    GMM = G[np.ix_(dof, dof)]
    X = np.linalg.solve(np.eye(len(dof)) + GMM @ dH, G[np.ix_(dof, rows)])
    Gtt = G[np.ix_(rows, rows)] - G[np.ix_(rows, dof)] @ dH @ X
    return np.einsum('iaib->i', Gtt.reshape(len(nodes), 3, len(nodes), 3))

def relative_msf_change(G, M_idx, dH, nodes):
    nodes = np.asarray(nodes, int)
    m0 = msf(G)[nodes]
    m1 = msf_after(G, M_idx, dH, nodes)
    return float(np.mean(np.log(m1 / m0)))

def covariance_direct_augmented(P, S, cutoff_p=15.0, cut_ps=10.0, cut_ss=10.0, g=1.0):
    P = np.asarray(P, float)
    S = np.asarray(S, float)
    X = np.vstack([P, S])
    N = len(P)
    ei, ej, eg = ([], [], [])
    i, j, gg = pair_edges(P, cutoff_p, g)
    ei.append(i)
    ej.append(j)
    eg.append(gg)
    d2 = ((P[:, None, :] - S[None, :, :]) ** 2).sum(-1)
    a, b = np.nonzero(d2 <= cut_ps ** 2)
    ei.append(a)
    ej.append(b + N)
    eg.append(np.full(len(a), g))
    ds = ((S[:, None, :] - S[None, :, :]) ** 2).sum(-1)
    a, b = np.nonzero(np.triu(ds <= cut_ss ** 2, 1))
    ei.append(a + N)
    ej.append(b + N)
    eg.append(np.full(len(a), g))
    H = hessian(X, np.concatenate(ei), np.concatenate(ej), np.concatenate(eg))
    w, V = modes(H)
    G = covariance(w, V)
    return G[:3 * N, :3 * N]

def lipid_anchors(ca_xyz_opm, cov_rows, slab=15.2, rsa_min=0.05, neigh=12.0):
    res = np.array(sorted(ca_xyz_opm))
    X = np.array([ca_xyz_opm[int(r)] for r in res])
    keep = []
    out = []
    for k, r in enumerate(res):
        c = cov_rows.get(int(r))
        if c is None or c['pool'] != 'TM' or (not c['helix']) or (c['rsa'] < rsa_min) or (abs(X[k, 2]) > slab):
            continue
        d = np.linalg.norm(X - X[k], axis=1)
        nb = X[(d < neigh) & (d > 0)]
        if len(nb) < 3:
            continue
        v = X[k, :2] - nb[:, :2].mean(0)
        n = np.linalg.norm(v)
        if n < 1e-06:
            continue
        keep.append(X[k])
        out.append(np.array([v[0] / n, v[1] / n, 0.0]))
    return (np.array(keep), np.array(out))

def place_templates(templates, anchors, outward, heavy_xyz, rng, n_mol=None, contact=(3.0, 4.2), min_inter=3.0, start=14.0, step=0.5, max_tries=200):
    from scipy.spatial import cKDTree
    tree = cKDTree(heavy_xyz)
    from a12_lib import random_rotation
    n_mol = len(templates) if n_mol is None else n_mol
    placed, mol = ([], [])
    tries = 0
    while len(placed) < n_mol and tries < max_tries * max(n_mol, 1):
        tries += 1
        tpl = templates[int(rng.integers(len(templates)))]
        a = int(rng.integers(len(anchors)))
        Rm = random_rotation(rng)
        local = tpl @ Rm.T
        L = start
        ok = False
        while L > 0:
            pos = local + anchors[a] + outward[a] * L
            dmin = tree.query(pos)[0].min()
            if dmin < contact[1]:
                ok = dmin >= contact[0]
                break
            L -= step
        if not ok:
            continue
        if placed and any((np.linalg.norm(pos[:, None, :] - q[None, :, :], axis=2).min() < min_inter for q in placed)):
            continue
        placed.append(pos)
        mol.append(np.full(len(pos), len(placed) - 1))
    if len(placed) < n_mol:
        raise RuntimeError('could not place all molecules')
    return (np.vstack(placed), np.concatenate(mol))

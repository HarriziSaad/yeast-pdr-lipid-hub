import numpy as np
from scipy.linalg import cho_factor, solve_triangular, LinAlgError
from scipy.spatial import cKDTree
from a12_lib import random_rotation

def make_tree(heavy_xyz):
    return cKDTree(heavy_xyz)

def place_templates_fast(templates, anchors, outward, tree, rng, n_mol=None, contact=(3.0, 4.2), min_inter=3.0, start=14.0, step=0.5, max_tries=200):
    n_mol = len(templates) if n_mol is None else n_mol
    Ls = []
    L = start
    while L > 0:
        Ls.append(L)
        L -= step
    nL = len(Ls)
    placed = []
    n_placed = 0
    q_all = None
    tries = 0
    while n_placed < n_mol and tries < max_tries * max(n_mol, 1):
        tries += 1
        tpl = templates[int(rng.integers(len(templates)))]
        a = int(rng.integers(len(anchors)))
        Rm = random_rotation(rng)
        local = tpl @ Rm.T
        base = local + anchors[a]
        ou = outward[a]
        i = 0
        ok = False
        while i < nL:
            pos = base + ou * Ls[i]
            dmin = tree.query(pos)[0].min()
            if dmin < contact[1]:
                ok = dmin >= contact[0]
                break
            i += 1 + max(int((dmin - contact[1] - 1e-06) / step), 0)
        if not ok:
            continue
        if q_all is not None and np.linalg.norm(pos[:, None, :] - q_all[None, :, :], axis=2).min() < min_inter:
            continue
        placed.append(pos)
        n_placed += 1
        q_all = pos if q_all is None else np.vstack([q_all, pos])
    if n_placed < n_mol:
        raise RuntimeError('could not place all molecules')
    return (np.vstack(placed), np.concatenate([np.full(len(p), j) for j, p in enumerate(placed)]))

def _diag_sums(idx, blk, n):
    out = np.empty((n, 3, 3))
    for a in range(3):
        for b in range(3):
            out[:, a, b] = np.bincount(idx, weights=blk[:, a, b], minlength=n)
    return out

def sterol_block(S, cut_ss=10.0, g_ss=1.0):
    S = np.asarray(S, float)
    K = len(S)
    ds = ((S[:, None, :] - S[None, :, :]) ** 2).sum(-1)
    k1, k2 = np.nonzero(np.triu(ds <= cut_ss * cut_ss, 1))
    H4 = np.zeros((K, 3, K, 3))
    if len(k1):
        vs = S[k2] - S[k1]
        bs = (g_ss / (vs * vs).sum(1))[:, None, None] * vs[:, :, None] * vs[:, None, :]
        H4[k1, :, k2, :] = -bs
        H4[k2, :, k1, :] = -bs
        diag = _diag_sums(k1, bs, K) + _diag_sums(k2, bs, K)
        ar = np.arange(K)
        H4[ar, :, ar, :] += diag
    return H4.reshape(3 * K, 3 * K)

def candidate_rows(P, S, cut_max, tree_S=None):
    if tree_S is None:
        tree_S = cKDTree(S)
    dist = tree_S.query(P)[0]
    return np.nonzero(dist <= cut_max + 1e-06)[0]

def condensation_fast(P, S, cut_ps, Hss_ss, cand=None, g_ps=1.0, ridge=1e-09, stats=None):
    P = np.asarray(P, float)
    S = np.asarray(S, float)
    K = len(S)
    Pc = P if cand is None else P[cand]
    d2 = ((Pc[:, None, :] - S[None, :, :]) ** 2).sum(-1)
    ip, ks = np.nonzero(d2 <= cut_ps * cut_ps)
    if len(ip) == 0:
        return (np.zeros(0, int), np.zeros((0, 0)))
    rows = np.unique(ip)
    M_idx = rows if cand is None else cand[rows]
    nm = len(rows)
    mi = np.searchsorted(rows, ip)
    v = S[ks] - Pc[ip]
    dd2 = (v * v).sum(1)
    blk = (g_ps / dd2)[:, None, None] * v[:, :, None] * v[:, None, :]
    B4 = np.zeros((nm, 3, K, 3))
    B4[mi, :, ks, :] = -blk
    B = B4.reshape(3 * nm, 3 * K)
    D4 = np.zeros((nm, 3, nm, 3))
    arm = np.arange(nm)
    D4[arm, :, arm, :] = _diag_sums(mi, blk, nm)
    D = D4.reshape(3 * nm, 3 * nm)
    Hss = Hss_ss.copy()
    Hv = Hss.reshape(K, 3, K, 3)
    ark = np.arange(K)
    Hv[ark, :, ark, :] += _diag_sums(ks, blk, K)
    Hss[np.diag_indices(3 * K)] += ridge
    try:
        c = cho_factor(Hss, lower=True, overwrite_a=True, check_finite=False)[0]
        Y = solve_triangular(c, B.T, lower=True, check_finite=False)
        dH = D - Y.T @ Y
    except LinAlgError:
        if stats is not None:
            stats['chol_fallback'] = stats.get('chol_fallback', 0) + 1
        Hss = Hss_ss.copy()
        Hv = Hss.reshape(K, 3, K, 3)
        Hv[ark, :, ark, :] += _diag_sums(ks, blk, K)
        Hss[np.diag_indices(3 * K)] += ridge
        dH = D - B @ np.linalg.solve(Hss, B.T)
    dH = 0.5 * (dH + dH.T)
    return (M_idx, dH)

def msf_prime(G, M_idx, dH, nodes):
    nodes = np.asarray(nodes, int)
    rows = (3 * nodes[:, None] + np.arange(3)[None, :]).ravel()
    if len(M_idx) == 0:
        return np.einsum('iaib->i', G[np.ix_(rows, rows)].reshape(len(nodes), 3, len(nodes), 3))
    dof = (3 * np.asarray(M_idx)[:, None] + np.arange(3)[None, :]).ravel()
    GMM = G[np.ix_(dof, dof)]
    X = np.linalg.solve(np.eye(len(dof)) + GMM @ dH, G[np.ix_(dof, rows)])
    Gtt = G[np.ix_(rows, rows)] - G[np.ix_(rows, dof)] @ dH @ X
    return np.einsum('iaib->i', Gtt.reshape(len(nodes), 3, len(nodes), 3))

def msf_diag(G):
    N = G.shape[0] // 3
    return np.einsum('iaib->i', np.asarray(G).reshape(N, 3, N, 3))

def lnratio_nodes(G, msf0_nodes, S, P, cuts, nodes_union, Hss_ss=None, cut_ss=10.0, g_ss=1.0, tree_S=None, cand=None, stats=None):
    S = np.asarray(S, float)
    if Hss_ss is None:
        Hss_ss = sterol_block(S, cut_ss, g_ss)
    if cand is None:
        cand = candidate_rows(P, S, max(cuts), tree_S)
    out = np.empty((len(cuts), len(nodes_union)))
    nm = np.empty(len(cuts), int)
    for j, cps in enumerate(cuts):
        M, dH = condensation_fast(P, S, cps, Hss_ss, cand, stats=stats)
        out[j] = np.log(msf_prime(G, M, dH, nodes_union) / msf0_nodes)
        nm[j] = len(M)
    return (out, nm)
_W = {}

def init_worker(G_path, P, nodes_union, m0u, cuts, cut_ss, anchors, outward, heavy_xyz):
    _W['G'] = np.load(G_path, mmap_mode='r') if isinstance(G_path, str) else G_path
    _W.update(P=P, nodes_union=nodes_union, m0u=m0u, cuts=list(cuts), cut_ss=cut_ss, anchors=anchors, outward=outward, tree=make_tree(heavy_xyz))

def task_place(args):
    from a12_lib import make_rng
    entry, seed_offset, templates, nd, n_mol = args
    rng = make_rng(seed_offset)
    out = [place_templates_fast(templates, _W['anchors'], _W['outward'], _W['tree'], rng, n_mol=n_mol)[0] for _ in range(nd)]
    return (entry, out)

def task_eval(args):
    entry, start, Slist = args
    G = np.asarray(_W['G'])
    stats = {}
    L_out = np.empty((len(Slist), len(_W['cuts']), len(_W['nodes_union'])))
    N_out = np.empty((len(Slist), len(_W['cuts'])), int)
    for i, S in enumerate(Slist):
        L_out[i], N_out[i] = lnratio_nodes(G, _W['m0u'], S, _W['P'], _W['cuts'], _W['nodes_union'], cut_ss=_W['cut_ss'], stats=stats)
    return (entry, start, L_out, N_out, stats)

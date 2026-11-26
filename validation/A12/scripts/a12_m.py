import numpy as np
from a12_lib import SEED, bootstrap_delta, candidate_sets, draw_sets, effect_summary, make_rng, match_report, perm_p

def _hit_vector(cov_res_sorted, hits):
    return np.array([1.0 if r in hits else 0.0 for r in cov_res_sorted])

def _index(cov):
    res = np.array(sorted(cov))
    pos = {int(r): i for i, r in enumerate(res)}
    return (res, pos)

def m1a(region, hits, cov, n_draws=10000, seed_offset=0, tested=None, n_boot=2000, tiers=None):
    reg = [r for r in sorted(set(region)) if r in cov and (tested is None or r in tested)]
    uni = [r for r in cov if tested is None or r in tested]
    if len(reg) < 10:
        nan = float('nan')
        return dict(obs=nan, null_mean=nan, null_sd=nan, delta=nan, z=nan, mde80=nan, n_eff=len(reg), p=nan, hit_fraction=nan, null_fraction=nan, fold=nan, match=dict(frac_fallback=nan, weakly_matched=False, n=len(reg)), ci95_delta_fraction=(nan, nan), null_dist=np.zeros(0), balance={}, not_testable=True)
    res, pos = _index(cov)
    hv = _hit_vector(res, hits)
    cs = candidate_sets(reg, uni, cov, tiers=tiers)
    rng = make_rng(seed_offset)
    draws, forced = draw_sets({t: cs[t]['cands'] for t in reg}, reg, n_draws, rng)
    idx = np.vectorize(pos.get, otypes=[int])(draws)
    ctrl = hv[idx]
    obs_vec = np.array([1.0 if r in hits else 0.0 for r in reg])
    obs = float(obs_vec.sum())
    null = ctrl.sum(axis=1)
    out = effect_summary(obs, null)
    out.update(n_eff=len(reg), p=perm_p(obs, null, 'upper'), hit_fraction=obs / max(len(reg), 1), null_fraction=float(null.mean()) / max(len(reg), 1), fold=obs / null.mean() if null.mean() > 0 else float('inf'), match=match_report(cs, forced))
    out['ci95_delta_fraction'] = bootstrap_delta(obs_vec, ctrl, n_boot) if len(reg) else (float('nan'),) * 2
    out['null_dist'] = null
    out['balance'] = balance_report(reg, draws, cov) if len(reg) else {}
    return out

def m1b(region, cavity, hits, cov, n_perm=10000, seed_offset=0, n_boot=2000, tiers=None):
    R0, C0 = (set(region), set(cavity))
    ov = R0 & C0
    R = [r for r in sorted(R0 - ov) if r in cov]
    C = [r for r in sorted(C0 - ov) if r in cov]
    excl = R0 | C0
    uni = [r for r in cov if r not in excl]
    hv = {r: 1.0 if r in hits else 0.0 for r in cov}
    csR = candidate_sets(R, uni, cov, tiers=tiers)
    csC = candidate_sets(C, uni, cov, tiers=tiers)
    bR = np.array([np.mean([hv[c] for c in csR[r]['cands']]) for r in R])
    bC = np.array([np.mean([hv[c] for c in csC[r]['cands']]) for r in C])
    yR = np.array([hv[r] for r in R]) - bR
    yC = np.array([hv[r] for r in C]) - bC
    D = float(yR.mean() - yC.mean())
    y = np.concatenate([yR, yC])
    lab = np.array([1] * len(R) + [0] * len(C))
    keys = [(cov[r]['pool'], cov[r]['helix'], cov[r]['exposed']) for r in R + C]
    groups = {}
    for i, k in enumerate(keys):
        groups.setdefault(k, []).append(i)
    groups = [np.array(g) for g in groups.values()]
    rng = make_rng(seed_offset)
    null = np.empty(n_perm)
    for b in range(n_perm):
        lp = lab.copy()
        for g in groups:
            lp[g] = lab[g][rng.permutation(len(g))]
        null[b] = y[lp == 1].mean() - y[lp == 0].mean()
    sd = float(np.std(null, ddof=1))
    rb = np.random.default_rng(SEED + 1)
    boots = []
    for _ in range(n_boot):
        boots.append(yR[rb.integers(0, len(R), len(R))].mean() - yC[rb.integers(0, len(C), len(C))].mean())
    ci = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))) if n_boot else (float('nan'),) * 2
    return dict(n_R=len(R), n_C=len(C), n_overlap_removed=len(ov), hit_frac_R=float(np.mean([hv[r] for r in R])), hit_frac_C=float(np.mean([hv[r] for r in C])), expected_frac_R=float(bR.mean()), expected_frac_C=float(bC.mean()), excess_R=float(yR.mean()), excess_C=float(yC.mean()), delta_RC=D, null_sd=sd, z=D / sd if sd > 0 else float('nan'), p=perm_p(D, null, 'upper'), ci95_delta=ci, mde80=2.486 * sd, n_strata=len(groups), matchR=match_report(csR), matchC=match_report(csC), null_dist=null)

def label_permutation_check(region, cavity, hits, cov, n_perm=10000, seed_offset=0):
    R0, C0 = (set(region), set(cavity))
    ov = R0 & C0
    R = [r for r in sorted(R0 - ov) if r in cov]
    C = [r for r in sorted(C0 - ov) if r in cov]
    allr = R + C
    rsa = np.array([cov[r]['rsa'] for r in allr])
    cuts = np.quantile(rsa, [1 / 3, 2 / 3])
    key = [(cov[r]['pool'], cov[r]['exposed'], cov[r]['helix'], int(np.searchsorted(cuts, cov[r]['rsa']))) for r in allr]
    lab = np.array([1] * len(R) + [0] * len(C))
    h = np.array([1.0 if r in hits else 0.0 for r in allr])
    groups = {}
    for i, k in enumerate(key):
        groups.setdefault(k, []).append(i)
    obs = h[lab == 1].mean() - h[lab == 0].mean()
    rng = make_rng(seed_offset + 11)
    null = np.empty(n_perm)
    for b in range(n_perm):
        lp = lab.copy()
        for g in groups.values():
            lp[g] = lab[g][rng.permutation(len(g))]
        null[b] = h[lp == 1].mean() - h[lp == 0].mean()
    return dict(obs=float(obs), p=perm_p(obs, null, 'upper'), n_strata=len(groups))

def _design(cov, residues):
    pools = ('TM', 'TMDx', 'NBD', 'ECD', 'NTE')
    rows = []
    for r in residues:
        c = cov[r]
        zz = c['z'] / 20.0
        rows.append([1.0] + [1.0 if c['pool'] == p else 0.0 for p in pools[1:]] + [float(c['helix']), float(c['exposed']), c['rsa'], c['rsa'] ** 2, zz if c['pool'] in ('TM', 'TMDx') else 0.0, zz * zz if c['pool'] in ('TM', 'TMDx') else 0.0])
    return np.array(rows)

def fit_logistic(Xd, y, ridge=0.001, n_iter=60):
    beta = np.zeros(Xd.shape[1])
    for _ in range(n_iter):
        eta = np.clip(Xd @ beta, -30, 30)
        p = 1.0 / (1.0 + np.exp(-eta))
        W = np.maximum(p * (1 - p), 1e-09)
        H = Xd.T @ (Xd * W[:, None]) + ridge * np.eye(Xd.shape[1])
        g = Xd.T @ (y - p) - ridge * beta
        step = np.linalg.solve(H, g)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-08:
            break
    return beta

def poisson_binomial_upper(p, k):
    pmf = np.zeros(len(p) + 1)
    pmf[0] = 1.0
    for pi in p:
        pmf[1:] = pmf[1:] * (1 - pi) + pmf[:-1] * pi
        pmf[0] *= 1 - pi
    return float(pmf[int(k):].sum())

def m1a_logistic(region, hits, cov, tested=None):
    reg = [r for r in sorted(set(region)) if r in cov and (tested is None or r in tested)]
    uni = [r for r in sorted(cov) if r not in set(reg) and (tested is None or r in tested)]
    Xu = _design(cov, uni)
    yu = np.array([1.0 if r in hits else 0.0 for r in uni])
    beta = fit_logistic(Xu, yu)
    pr = 1.0 / (1.0 + np.exp(-np.clip(_design(cov, reg) @ beta, -30, 30)))
    obs = sum((1 for r in reg if r in hits))
    return dict(n_eff=len(reg), obs=obs, expected=float(pr.sum()), p=poisson_binomial_upper(pr, obs))

def balance_report(targets, draws, cov):
    tg = list(targets)
    allr = list(cov)
    sd_rsa = float(np.std([cov[r]['rsa'] for r in allr]))
    sd_z = float(np.std([cov[r]['z'] for r in allr]))
    rsa_t = np.mean([cov[r]['rsa'] for r in tg])
    z_t = np.mean([cov[r]['z'] for r in tg])
    rsa_c = float(np.mean(np.vectorize(lambda r: cov[r]['rsa'])(draws)))
    z_c = float(np.mean(np.vectorize(lambda r: cov[r]['z'])(draws)))
    out = dict(smd_rsa=(rsa_t - rsa_c) / sd_rsa, smd_depth=(z_t - z_c) / sd_z)
    out['imbalanced'] = bool(abs(out['smd_rsa']) > 0.2 or abs(out['smd_depth']) > 0.2)
    return out

def pair_counts(res_matrix, dist, pos, r_cut, min_sep=5):
    D, n = res_matrix.shape
    idx = np.vectorize(pos.get, otypes=[int])(res_matrix)
    iu, ju = np.triu_indices(n, 1)
    a, b = (idx[:, iu], idx[:, ju])
    close = dist[a, b] <= r_cut
    sep = np.abs(res_matrix[:, iu] - res_matrix[:, ju]) >= min_sep
    return (close & sep).sum(axis=1)

def mean_nn(res_matrix, dist, pos, min_sep=5):
    D, n = res_matrix.shape
    idx = np.vectorize(pos.get, otypes=[int])(res_matrix)
    out = np.empty(D)
    for d in range(D):
        sub = dist[np.ix_(idx[d], idx[d])].copy()
        rs = res_matrix[d]
        bad = np.abs(rs[:, None] - rs[None, :]) < min_sep
        sub[bad] = np.inf
        out[d] = sub.min(axis=1).mean()
    return out

def m2(hit_residues, ca_xyz, cov, n_draws=10000, seed_offset=0, radii=(10.0,), primary_r=10.0, min_sep=5):
    A = [r for r in sorted(set(hit_residues)) if r in cov and r in ca_xyz]
    if len(A) < 3:
        nan = dict(obs=float('nan'), null_mean=float('nan'), null_sd=float('nan'), delta=float('nan'), z=float('nan'), mde80=float('nan'), p=float('nan'))
        return dict(n_eff=len(A), match=dict(frac_fallback=float('nan'), weakly_matched=False), radii={r: dict(nan) for r in sorted(set(radii) | {primary_r})}, mean_nn=dict(nan), primary=dict(nan), null_dist_primary=np.zeros(0))
    uni = [r for r in cov if r in ca_xyz]
    allr = np.array(sorted((r for r in cov if r in ca_xyz)))
    pos = {int(r): i for i, r in enumerate(allr)}
    X = np.array([ca_xyz[int(r)] for r in allr], float)
    dist = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    cs = candidate_sets(A, uni, cov)
    draws, forced = draw_sets({t: cs[t]['cands'] for t in A}, A, n_draws, make_rng(seed_offset))
    obs_mat = np.array([A])
    out = dict(n_eff=len(A), match=match_report(cs, forced), radii={})
    for r_cut in sorted(set(radii) | {primary_r}):
        obs = float(pair_counts(obs_mat, dist, pos, r_cut, min_sep)[0])
        null = pair_counts(draws, dist, pos, r_cut, min_sep)
        s = effect_summary(obs, null)
        s['p'] = perm_p(obs, null, 'upper')
        out['radii'][r_cut] = s
    o = float(mean_nn(obs_mat, dist, pos, min_sep)[0])
    nn = mean_nn(draws, dist, pos, min_sep)
    s = effect_summary(o, nn)
    s['p'] = perm_p(o, nn, 'lower')
    out['mean_nn'] = s
    out['primary'] = out['radii'][primary_r]
    out['null_dist_primary'] = pair_counts(draws, dist, pos, primary_r, min_sep)
    return out

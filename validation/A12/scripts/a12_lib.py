import csv
import hashlib
import math
import os
import numpy as np
SEED = 20261118
RSA_EXPOSED = 0.05
TIERS_BASE = ((0.1, 3.0), (0.15, 4.5), (0.2, 6.0))
MIN_CAND = 10
WEAK_FALLBACK_FRACTION = 0.25
TM_POOLS = ('TM', 'TMDx')

def sha256_file(path, block=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        while True:
            b = fh.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def read_covariates(path):
    cov = {}
    with open(path, encoding='utf-8', newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            res = int(r['residue'])
            rsa = float(r['RSA'])
            cov[res] = dict(rsa=rsa, z=float(r['z_CA_OPM_A']), helix=int(r['in_HELX_P']), pool=r['pool_label'], domain=r['domain'], exposed=int(rsa >= RSA_EXPOSED))
    return cov

def read_set_table(path):
    out = {}
    with open(path, encoding='utf-8', newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            out.setdefault(r['set'], set()).add(int(r['residue']))
    return out

class CovArrays:

    def __init__(self, cov, residues):
        self.res = np.array(sorted(residues), dtype=int)
        self.rsa = np.array([cov[r]['rsa'] for r in self.res])
        self.z = np.array([cov[r]['z'] for r in self.res])
        self.helix = np.array([cov[r]['helix'] for r in self.res])
        self.exposed = np.array([cov[r]['exposed'] for r in self.res])
        self.pool = np.array([cov[r]['pool'] for r in self.res])
        self.domain = np.array([cov[r]['domain'] for r in self.res])

def candidate_sets(targets, universe, cov, extra=(), min_cand=MIN_CAND, include_self=True, tiers=None):
    tiers = TIERS_BASE if tiers is None else tiers
    uni = [u for u in universe if u in cov]
    ca = CovArrays(cov, uni)
    ex = []
    for vals, cal in extra:
        ex.append((vals, np.array([vals.get(int(r), np.nan) for r in ca.res]), cal))
    out = {}
    for t in targets:
        ct = cov[t]
        base = (ca.pool == ct['pool']) & (ca.helix == ct['helix'])
        if not include_self:
            base &= ca.res != t
        strict = base & (ca.exposed == ct['exposed'])
        if ct['pool'] in TM_POOLS:
            geo = lambda tier: np.abs(ca.z - ct['z']) <= tiers[tier][1]
        else:
            geo = lambda tier: ca.domain == ct['domain']
        chosen, tier_used = (None, 3)
        for tier in range(3):
            m = strict & geo(tier) & (np.abs(ca.rsa - ct['rsa']) <= tiers[tier][0])
            for vals, arr, cal in ex:
                m &= np.abs(arr - vals.get(int(t), np.nan)) <= cal[tier]
            if m.sum() >= min_cand:
                chosen, tier_used = (ca.res[m], tier)
                break
        if chosen is None:
            idx = np.where(base)[0]
            d2 = ((ca.rsa[idx] - ct['rsa']) / tiers[0][0]) ** 2
            if ct['pool'] in TM_POOLS:
                d2 = d2 + ((ca.z[idx] - ct['z']) / tiers[0][1]) ** 2
            else:
                d2 = d2 + (ca.domain[idx] != ct['domain']) * 100.0
            for vals, arr, cal in ex:
                d2 = d2 + ((arr[idx] - vals.get(int(t), np.nan)) / cal[0]) ** 2
            d2 = np.where(np.isnan(d2), np.inf, d2)
            order = np.argsort(d2, kind='stable')[:min_cand]
            chosen = ca.res[idx[order]]
        out[t] = dict(cands=np.asarray(chosen, dtype=int), tier=tier_used, n=len(chosen), fallback=tier_used == 3)
    return out

def draw_sets(cands_by_target, order, n_draws, rng):
    n = len(order)
    cl = [np.asarray(cands_by_target[t], dtype=int) for t in order]
    out = np.zeros((n_draws, n), dtype=int)
    forced = 0
    for d in range(n_draws):
        used = set()
        for i in rng.permutation(n):
            c = cl[i]
            free = [x for x in c.tolist() if x not in used]
            if not free:
                free = c.tolist()
                forced += 1
            pick = free[int(rng.integers(len(free)))]
            used.add(pick)
            out[d, i] = pick
    return (out, forced)

def match_report(cs, n_forced=0):
    n = len(cs)
    nfb = sum((1 for v in cs.values() if v['fallback']))
    sizes = [v['n'] for v in cs.values()]
    tiers = [v['tier'] for v in cs.values()]
    return dict(n=n, n_fallback=nfb, frac_fallback=nfb / n if n else float('nan'), weakly_matched=bool(n and nfb / n > WEAK_FALLBACK_FRACTION), cand_min=int(min(sizes)) if sizes else 0, cand_median=float(np.median(sizes)) if sizes else 0.0, tiers=[tiers.count(k) for k in range(4)], forced_reuse=int(n_forced))

def perm_p(obs, null, tail='upper'):
    null = np.asarray(null, dtype=float)
    B = len(null)
    if tail == 'upper':
        b = int(np.sum(null >= obs - 1e-12))
    elif tail == 'lower':
        b = int(np.sum(null <= obs + 1e-12))
    else:
        raise ValueError(tail)
    return (1.0 + b) / (1.0 + B)

def effect_summary(obs, null):
    null = np.asarray(null, dtype=float)
    sd = float(np.std(null, ddof=1)) if len(null) > 1 else float('nan')
    delta = float(obs - np.mean(null))
    return dict(obs=float(obs), null_mean=float(np.mean(null)), null_sd=sd, delta=delta, z=delta / sd if sd and sd > 0 else float('nan'), mde80=2.486 * sd)

def bootstrap_delta(class_vals, ctrl_vals, n_boot=2000, seed=SEED + 1):
    class_vals = np.asarray(class_vals, dtype=float)
    ctrl_vals = np.asarray(ctrl_vals, dtype=float)
    rng = np.random.default_rng(seed)
    n = len(class_vals)
    colmean = ctrl_vals.mean(axis=0)
    res = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        res[b] = class_vals[idx].mean() - colmean[idx].mean()
    return (float(np.percentile(res, 2.5)), float(np.percentile(res, 97.5)))

def bh_qvalues(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.minimum(q, 1.0)
    out = np.empty(n)
    out[order] = q
    return out

def make_rng(seed_offset=0):
    return np.random.default_rng(SEED + seed_offset)

def kabsch(P, Q):
    P = np.asarray(P, float)
    Q = np.asarray(Q, float)
    pc, qc = (P.mean(0), Q.mean(0))
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = qc - R @ pc
    rmsd = float(np.sqrt(np.mean(np.sum((P @ R.T + t - Q) ** 2, axis=1))))
    return (R, t, rmsd)

def random_rotation(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)], [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)], [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d]])

def pairwise_dist(X):
    X = np.asarray(X, float)
    d = X[:, None, :] - X[None, :, :]
    return np.sqrt((d * d).sum(-1))

def guard_amendment_appended(amendments_path, token='## A12:'):
    if not os.path.exists(amendments_path):
        raise SystemExit('AMENDMENTS.md not found: refusing to run on real data')
    with open(amendments_path, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith(token):
                return True
    raise SystemExit('A12 is not appended to AMENDMENTS.md: refusing to run on real data (pre-registration rule)')

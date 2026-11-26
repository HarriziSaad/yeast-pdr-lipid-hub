import csv
import os
import sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.path.insert(0, 'C:\\Users\\MSI\\Documents\\Advisor\\project\\phase5_validation\\A12\\scripts')
import collections
import time
import numpy as np
from scipy.spatial import cKDTree
import a12_enm as E
import a12_lib as L
import a12_run_b as RB
import a12_b2_fast as F
ENTRIES = ('9S2H', '9S4S', '9S4T', '9S4U', '9S2C')
NSUB = 300
B_REG = 1000

def main(root, b2dir):
    out = []

    def log(m=''):
        print(m, flush=True)
        out.append(m)
    t0 = time.time()
    cov, sets, structs, hits = RB.load_b_inputs(root)
    qc, ster = RB.build_sterols(root, None)
    s = structs['7P04']
    P, res = (s['X'], s['res'])
    gpath = 'C:/Users/MSI/AppData/Local/Temp/claude/6ff21cdddf57ef83/b10be529-0824-468f-9fda-f7f9216985a8/scratchpad/G_real_timing.npy'
    if os.path.exists(gpath):
        G = np.load(gpath)
        log(f'G loaded from {gpath}')
    else:
        G = E.covariance(*E.modes(E.anm(P, 15.0)))
    pos = {r: i for i, r in enumerate(res)}
    gate = np.array([pos[r] for r in sorted(sets['GATE_ALL']) if r in pos])
    m0 = F.msf_diag(G)[gate]
    anchors, outward = E.lipid_anchors(s['ca'], cov)
    treeP = cKDTree(P)
    treeH = F.make_tree(s['heavy'])
    obs, n_obs = ({}, {})
    for r in csv.DictReader(open(os.path.join(b2dir, 'a12_B2_results.tsv'), encoding='utf-8'), delimiter='\t'):
        if r['variant'] in ('PRIMARY', 'outward-facing entry'):
            obs[r['structure']] = float(r['obs'])
            n_obs[r['structure']] = int(r['note'].split()[2].rstrip(';'))
    nul = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(b2dir, 'a12_B2_null_draws_nd2000_clean.tsv'), encoding='utf-8'), delimiter='\t'):
        if r['variant'] in ('PRIMARY', 'outward-facing entry') and int(r['draw']) < B_REG:
            nul[r['structure']].append((int(r['n_contact_nodes']), float(r['mean_ln_msf_ratio'])))
    log("B2 contact-extent-matched POST HOC sensitivity (Step 6c); statistic = mean ln(MSF'/MSF) over GATE-ALL, sterol-protein cutoff 10 A, lower tail")
    log(f'(i) feasibility: null draws (first {B_REG}) with n10 within +/-10 % of the observed n10')
    res_ii, res_iii = ({}, {})
    for e in ENTRIES:
        n = np.array([a for a, b in nul[e]], float)
        sv = np.array([b for a, b in nul[e]], float)
        lo, hi = (0.9 * n_obs[e], 1.1 * n_obs[e])
        k_match = int(np.sum((n >= lo) & (n <= hi)))
        log(f'  {e}: observed n10 {n_obs[e]}; null n10 mean {n.mean():.1f}, min {int(n.min())}, max {int(n.max())}; draws within +/-10 % of observed: {k_match}; extrapolation factor n10_obs / max null n10 = {n_obs[e] / n.max():.2f}')
        x, y = (np.log(n), np.log(-sv))
        b, a = np.polyfit(x, y, 1)
        r = y - (a + b * x)
        ro = np.log(-obs[e]) - (a + b * np.log(n_obs[e]))
        p_log = (1 + np.sum(r >= ro)) / (1 + len(r))
        pred = -np.exp(a + b * np.log(n_obs[e]))
        z_log = ro / r.std(ddof=1)
        bl, al = np.polyfit(n, sv, 1)
        rl = sv - (al + bl * n)
        rol = obs[e] - (al + bl * n_obs[e])
        p_lin = (1 + np.sum(rl <= rol)) / (1 + len(rl))
        pred_l = al + bl * n_obs[e]
        res_ii[e] = dict(obs=obs[e], pred_log=pred, p_log=p_log, z_log=z_log, slope=b, pred_lin=pred_l, p_lin=p_lin, z_lin=rol / rl.std(ddof=1), r2_log=1 - r.var() / y.var())
    log('(ii) covariate-adjusted null (OLS of ln(-s) on ln(n10) over the registered null draws; EXTRAPOLATED to the observed n10)')
    log('entry | observed s | fitted null mean at observed n10 (log-log) | z of residual | p (log-log) | fitted null mean (linear) | p (linear) | slope of ln(-s) on ln(n10) | R2 of null fit')
    for e in ENTRIES:
        d = res_ii[e]
        log(f"  {e} | {d['obs']:.5f} | {d['pred_log']:.5f} | {d['z_log']:.1f} | {d['p_log']:.6f} | {d['pred_lin']:.5f} | {d['p_lin']:.6f} | {d['slope']:.2f} | {d['r2_log']:.2f}")
    log('(iii) interpolation: random sub-arrangements of the observed molecules inside the null n10 range, compared with null draws of the same n10 (+/-10 %)')
    log('entry | subsets drawn | kept (n10 inside null range, >= 20 matched null draws) | median molecules per kept subset | median p_sub | fraction p_sub < 0.1 | min p_sub | max p_sub')
    for e in ENTRIES:
        xyz, mol = (ster[e]['xyz'], ster[e]['mol'])
        mols = np.unique(mol)
        n = np.array([a for a, b in nul[e]], float)
        sv = np.array([b for a, b in nul[e]], float)
        rng = np.random.default_rng(L.SEED + sum(map(ord, e)) % 1000 + 9000)
        kept = []
        for _ in range(NSUB):
            k = int(rng.integers(2, len(mols)))
            sel = rng.choice(mols, size=k, replace=False)
            S = xyz[np.isin(mol, sel)]
            o, nm = F.lnratio_nodes(G, m0, S, P, [10.0], gate)
            s_sub = float(np.mean(o[0]))
            n10 = int(nm[0])
            if not n.min() <= n10 <= n.max():
                continue
            mt = (n >= 0.9 * n10) & (n <= 1.1 * n10)
            if mt.sum() < 20:
                continue
            kept.append((k, n10, s_sub, int(mt.sum()), (1 + np.sum(sv[mt] <= s_sub)) / (1 + mt.sum())))
        ps = np.array([q[4] for q in kept]) if kept else np.array([np.nan])
        res_iii[e] = dict(kept=len(kept), med=float(np.median(ps)), frac=float(np.mean(ps < 0.1)) if kept else float('nan'))
        log(f"  {e} | {NSUB} | {len(kept)} | {(np.median([q[0] for q in kept]) if kept else float('nan'))} | {np.median(ps):.4f} | {(np.mean(ps < 0.1) if kept else float('nan')):.3f} | {np.nanmin(ps):.4f} | {np.nanmax(ps):.4f}")
        with open(os.path.join(b2dir, f'B2_contact_matched_subsets_{e}.tsv'), 'w', encoding='utf-8', newline='') as fh:
            fh.write('n_molecules\tn10\ts_sub\tn_matched_null\tp_sub\n')
            for q in kept:
                fh.write(f'{q[0]}\t{q[1]}\t{q[2]!r}\t{q[3]}\t{q[4]:.6f}\n')
    log('Conclusion by the declared rule (survives for an entry if (ii) log-log p < 0.1 and, where >= 20 kept subsets, median p_sub < 0.1; overall if in >= 3 of the 4 primary entries)')
    surv = {}
    for e in ENTRIES:
        ok2 = res_ii[e]['p_log'] < 0.1
        ok3 = res_iii[e]['med'] < 0.1 if res_iii[e]['kept'] >= 20 else True
        surv[e] = bool(ok2 and ok3)
        log(f"  {e}: (ii) p_log {res_ii[e]['p_log']:.6f} -> {ok2}; (iii) kept {res_iii[e]['kept']}, median p_sub {res_iii[e]['med']:.4f} -> {(ok3 if res_iii[e]['kept'] >= 20 else 'n/a (too few subsets)')}; survives: {surv[e]}")
    prim = ('9S2H', '9S4S', '9S4T', '9S4U')
    log(f"  survives overall (>= 3 of 4 primary entries): {sum((surv[e] for e in prim)) >= 3} ({sum((surv[e] for e in prim))} of 4); 9S2C: {surv['9S2C']}")
    log('Descriptive (no test): per-molecule contact nodes at 10 A, observed molecules vs 300 single random placements (registered placement mechanism, stream seed 777), and observed sterol atoms closer than 3.0 A to a protein heavy atom')
    for e in ENTRIES:
        xyz, mol = (ster[e]['xyz'], ster[e]['mol'])
        cnt = lambda S: len(set(np.concatenate(treeP.query_ball_point(S, 10.0))))
        per = [cnt(xyz[mol == k]) for k in np.unique(mol)]
        tpl = [xyz[mol == k] - xyz[mol == k].mean(0) for k in np.unique(mol)]
        rng = L.make_rng(777)
        single = [cnt(F.place_templates_fast(tpl, anchors, outward, treeH, rng, n_mol=1)[0]) for _ in range(300)]
        dmin = treeH.query(xyz)[0]
        log(f'  {e}: observed per-molecule median {np.median(per):.0f} (range {min(per)}-{max(per)}); random single molecule median {np.median(single):.0f} (90th pct {np.percentile(single, 90):.0f}, max {max(single)}); observed atoms with nearest protein heavy atom < 3.0 A: {int(np.sum(dmin < 3.0))} of {len(xyz)} (min distance {dmin.min():.2f} A)')
    log(f'elapsed {time.time() - t0:.0f} s')
    with open(os.path.join(b2dir, 'B2_contact_matched_posthoc.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(out) + '\n')
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

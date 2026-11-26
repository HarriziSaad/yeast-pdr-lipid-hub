import numpy as np
import a12_enm as E
from a12_lib import bootstrap_delta, candidate_sets, draw_sets, effect_summary, make_rng, match_report, perm_p
CAL_DIST = (0.1, 0.15, 0.25)
CAL_SELF = (0.15, 0.25, 0.4)

def b1_features(residues, X, hub, targets, cutoff=15.0, gamma=1.0, n_modes=None, G=None):
    res = list(residues)
    pos = {r: i for i, r in enumerate(res)}
    if G is None:
        H = E.anm(X, cutoff, gamma)
        w, V = E.modes(H)
        G = E.covariance(w, V, n_modes)
    A = E.prs_amplitude(G)
    h = [pos[r] for r in hub if r in pos]
    t = [pos[r] for r in targets if r in pos]
    sH = E.sensitivity(A, h)
    eT = E.effectiveness(A, t)
    dH = E.min_distance_to(X, h)
    dT = E.min_distance_to(X, t)
    diag = np.diag(A)
    return dict(res=res, A=A, G=G, sH=dict(zip(res, sH)), eT=dict(zip(res, eT)), c=dict(zip(res, np.log10(sH) + np.log10(eT))), logdH=dict(zip(res, np.log10(np.maximum(dH, 0.001)))), logdT=dict(zip(res, np.log10(np.maximum(dT, 0.001)))), logself=dict(zip(res, np.log10(diag))))

def _extras(feat, with_distance=True):
    ex = [(feat['logself'], CAL_SELF)]
    if with_distance:
        ex = [(feat['logdH'], CAL_DIST), (feat['logdT'], CAL_DIST)] + ex
    return ex

def b1_as_test(as_residues, hub, targets, feat, cov, n_draws=10000, seed_offset=0, with_distance=True, n_boot=2000):
    excl = set(hub) | set(targets)
    A_ = [r for r in sorted(set(as_residues)) if r in cov and r in feat['c'] and (r not in excl)]
    uni = [r for r in cov if r in feat['c'] and r not in excl]
    cs = candidate_sets(A_, uni, cov, extra=_extras(feat, with_distance))
    draws, forced = draw_sets({t: cs[t]['cands'] for t in A_}, A_, n_draws, make_rng(seed_offset))
    cvec = np.vectorize(feat['c'].get)(draws)
    obs_vals = np.array([feat['c'][r] for r in A_])
    obs = float(obs_vals.mean())
    null = cvec.mean(axis=1)
    out = effect_summary(obs, null)
    out.update(n_eff=len(A_), p=perm_p(obs, null, 'upper'), match=match_report(cs, forced), ci95_delta=bootstrap_delta(obs_vals, cvec, n_boot) if len(A_) else (float('nan'),) * 2)
    out['null_dist'] = null
    return out

def b1_hub_test(hub, target, feat_fn_values, cov, X_by_res, residues_all, n_draws=10000, seed_offset=0, n_boot=2000):
    excl = set(hub) | set(target)
    H_ = [r for r in sorted(set(hub)) if r in cov and r in feat_fn_values['eT']]
    uni = [r for r in cov if r in feat_fn_values['eT'] and r not in excl]
    extra = [(feat_fn_values['logdT'], CAL_DIST), (feat_fn_values['logself'], CAL_SELF)]
    cs = candidate_sets(H_, uni, cov, extra=extra)
    draws, forced = draw_sets({t: cs[t]['cands'] for t in H_}, H_, n_draws, make_rng(seed_offset))
    vals = np.vectorize(lambda r: np.log10(feat_fn_values['eT'][r]))(draws)
    obs_vals = np.array([np.log10(feat_fn_values['eT'][r]) for r in H_])
    obs = float(obs_vals.mean())
    null = vals.mean(axis=1)
    out = effect_summary(obs, null)
    out.update(n_eff=len(H_), p=perm_p(obs, null, 'upper'), match=match_report(cs, forced), ci95_delta=bootstrap_delta(obs_vals, vals, n_boot) if len(H_) else (float('nan'),) * 2)
    out['null_dist'] = null
    return out

def target_features(residues, X, target, G):
    res = list(residues)
    pos = {r: i for i, r in enumerate(res)}
    A = E.prs_amplitude(G)
    t = [pos[r] for r in target if r in pos]
    eT = E.effectiveness(A, t)
    dT = E.min_distance_to(X, t)
    return dict(eT=dict(zip(res, eT)), logdT=dict(zip(res, np.log10(np.maximum(dT, 0.001)))), logself=dict(zip(res, np.log10(np.diag(A)))))

def b2_test(P, G, gate_nodes, observed_sterol_xyz, templates, anchors, outward, heavy_xyz, n_draws=2000, seed_offset=0, cut_ps=10.0, cut_ss=10.0, g_ps=1.0, g_ss=1.0, n_mol=None):
    M, dH = E.condensation(P, observed_sterol_xyz, cut_ps, cut_ss, g_ps, g_ss)
    obs = E.relative_msf_change(G, M, dH, gate_nodes)
    rng = make_rng(seed_offset)
    n_mol = len(templates) if n_mol is None else n_mol
    null = np.empty(n_draws)
    for d in range(n_draws):
        S, _ = E.place_templates(templates, anchors, outward, heavy_xyz, rng, n_mol=n_mol)
        Md, dHd = E.condensation(P, S, cut_ps, cut_ss, g_ps, g_ss)
        null[d] = E.relative_msf_change(G, Md, dHd, gate_nodes)
    out = effect_summary(obs, null)
    out.update(p=perm_p(obs, null, 'lower'), n_contact_nodes=int(len(M)))
    out['null_dist'] = null
    return out

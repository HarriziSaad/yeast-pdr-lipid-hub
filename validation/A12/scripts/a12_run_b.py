import csv
import os
import sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import a12_b as B
import a12_enm as E
import a12_lib as L
import a12_select as S
FIELDS = ['test', 'structure', 'variant', 'n_eff', 'obs', 'null_mean', 'delta', 'z', 'p', 'q', 'ci_lo', 'ci_hi', 'mde80', 'frac_fallback', 'note']
INWARD_ENTRIES = ('9S2H', '9S4S', '9S4T', '9S4U')
OUTWARD_ENTRIES = ('9S2C',)

def _row(**k):
    r = {f: '' for f in FIELDS}
    r.update(k)
    return r

def _ci(o):
    ci = o.get('ci95_delta', (float('nan'), float('nan')))
    return (ci[0], ci[1])

def run_b1(structures, cov, sets, hits_by_name, nd=10000, nb=2000, cutoffs=(15.0,), n_modes_list=(None,), log=print):
    H = sorted(sets['R_HUB'])
    gate = sorted(sets['GATE_CORE'])
    nbs = sorted(sets['NBS'])
    T_all = sorted(set(gate) | set(nbs))
    out = []
    primary_tag = next(iter(structures))
    for tag, st in structures.items():
        res, X = (st['res'], st['X'])
        for cutoff in cutoffs:
            Hm = E.anm(X, cutoff)
            w, V = E.modes(Hm)
            for nm in n_modes_list:
                G = E.covariance(w, V, nm)
                is_primary = tag == primary_tag and cutoff == 15.0 and (nm is None)
                var = 'PRIMARY' if is_primary else f"cutoff {cutoff:g}, modes {('all' if nm is None else nm)}"
                feat = B.b1_features(res, X, [r for r in H if r in set(res)], [r for r in T_all if r in set(res)], cutoff, G=G)
                hits_main = hits_by_name['U1 pooled']
                o = B.b1_as_test(hits_main, H, T_all, feat, cov, n_draws=nd, n_boot=nb)
                ci = _ci(o)
                out.append(_row(test='B1-AS', structure=tag, variant=var, n_eff=o['n_eff'], obs=o['obs'], null_mean=o['null_mean'], delta=o['delta'], z=o['z'], p=o['p'], ci_lo=ci[0], ci_hi=ci[1], mde80=o['mde80'], frac_fallback=o['match']['frac_fallback'], note='statistic = mean log10(s_H * e_T) over the altered-specificity residues; T = GATE-core U NBS-MOTIF'))
                for tname, tset in (('GATE-core', gate), ('NBS-MOTIF', nbs)):
                    tf = B.target_features(res, X, [r for r in tset if r in set(res)], G)
                    oh = B.b1_hub_test(H, tset, tf, cov, None, res, n_draws=nd, n_boot=nb)
                    ci = _ci(oh)
                    out.append(_row(test=f'B1-HUB-{tname}', structure=tag, variant=var, n_eff=oh['n_eff'], obs=oh['obs'], null_mean=oh['null_mean'], delta=oh['delta'], z=oh['z'], p=oh['p'], ci_lo=ci[0], ci_hi=ci[1], mde80=oh['mde80'], frac_fallback=oh['match']['frac_fallback'], note='statistic = mean log10 e_target over the hub residues (effector direction)'))
                if is_primary:
                    for hn in ('U1 Pdr5-own', 'U1 Cdr1-mapped', 'U2 pooled'):
                        o2 = B.b1_as_test(hits_by_name[hn], H, T_all, feat, cov, n_draws=nd, n_boot=nb)
                        out.append(_row(test='B1-AS', structure=tag, variant=f'hits {hn}', n_eff=o2['n_eff'], obs=o2['obs'], null_mean=o2['null_mean'], delta=o2['delta'], z=o2['z'], p=o2['p'], mde80=o2['mde80'], frac_fallback=o2['match']['frac_fallback']))
                    o3 = B.b1_as_test(hits_main, H, T_all, feat, cov, n_draws=nd, n_boot=nb, with_distance=False)
                    out.append(_row(test='B1-AS', structure=tag, variant='no distance calipers (shows the proximity effect)', n_eff=o3['n_eff'], obs=o3['obs'], null_mean=o3['null_mean'], delta=o3['delta'], z=o3['z'], p=o3['p'], frac_fallback=o3['match']['frac_fallback']))
                    for tname, tset in (('GATE-core only', gate), ('NBS-MOTIF only', nbs)):
                        f2 = B.b1_features(res, X, [r for r in H if r in set(res)], [r for r in tset if r in set(res)], cutoff, G=G)
                        o4 = B.b1_as_test(hits_main, H, tset, f2, cov, n_draws=nd, n_boot=nb)
                        out.append(_row(test='B1-AS', structure=tag, variant=f'target {tname}', n_eff=o4['n_eff'], obs=o4['obs'], null_mean=o4['null_mean'], delta=o4['delta'], z=o4['z'], p=o4['p'], frac_fallback=o4['match']['frac_fallback']))
                    for part, key in (('effector only (log10 e_T)', 'eT'), ('sensor only (log10 s_H)', 'sH')):
                        f3 = dict(feat)
                        f3['c'] = {r: np.log10(feat[key][r]) for r in feat['c']}
                        o5 = B.b1_as_test(hits_main, H, T_all, f3, cov, n_draws=nd, n_boot=nb)
                        out.append(_row(test='B1-AS', structure=tag, variant=part, n_eff=o5['n_eff'], obs=o5['obs'], null_mean=o5['null_mean'], delta=o5['delta'], z=o5['z'], p=o5['p'], frac_fallback=o5['match']['frac_fallback']))
    prim = [r for r in out if r['variant'] == 'PRIMARY' and r['structure'] == primary_tag]
    q = L.bh_qvalues([r['p'] for r in prim])
    for r, qq in zip(prim, q):
        r['q'] = qq
    return out

def run_b2(P, res, cov, ca_opm, sets, heavy_xyz, sterol_sets, nd=1000, cut_ps=10.0, cut_ss=10.0, cutoff=15.0, log=print):
    H = E.anm(P, cutoff)
    G = E.covariance(*E.modes(H))
    pos = {r: i for i, r in enumerate(res)}
    gate_all = [pos[r] for r in sorted(sets['GATE_ALL']) if r in pos]
    gate_core = [pos[r] for r in sorted(sets['GATE_CORE']) if r in pos]
    anchors, outward = E.lipid_anchors(ca_opm, cov)
    out = []
    for entry, d in sterol_sets.items():
        xyz, mol = (d['xyz'], d['mol'])
        templates = [xyz[mol == k] - xyz[mol == k].mean(0) for k in np.unique(mol)]
        for variant, nodes, cps in (('PRIMARY' if entry in INWARD_ENTRIES else 'outward-facing entry', gate_all, cut_ps), ('gate-core only', gate_core, cut_ps), ('sterol-protein cutoff 8', gate_all, 8.0), ('sterol-protein cutoff 12', gate_all, 12.0)):
            o = B.b2_test(P, G, nodes, xyz, templates, anchors, outward, heavy_xyz, n_draws=nd, seed_offset=sum(map(ord, entry)) % 1000, cut_ps=cps, cut_ss=cut_ss)
            out.append(_row(test='B2', structure=entry, variant=variant, n_eff=len(templates), obs=o['obs'], null_mean=o['null_mean'], delta=o['delta'], z=o['z'], p=o['p'], mde80=o['mde80'], note=f"contact nodes {o['n_contact_nodes']}; statistic = mean ln(MSF'/MSF) over gate nodes (lower tail)"))
    return out

def b2_consistency(rows):
    prim = [r for r in rows if r['test'] == 'B2' and r['variant'] == 'PRIMARY']
    k = sum((1 for r in prim if float(r['p']) < 0.1))
    neg = all((float(r['delta']) < 0 for r in prim))
    return dict(n_entries=len(prim), n_p_lt_0_1=k, all_delta_negative=neg, label='consistent' if len(prim) >= 4 and k >= 3 and neg else 'not consistent')

def write_rows(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow(r)

def selftest():
    import test_a12_toy as T
    nres = 120
    rng = np.random.default_rng(31)
    X = T.synth_protein(nres, 32, 30.0)
    res = list(range(1, nres + 1))
    cov = {r: dict(rsa=float(np.clip(rng.beta(1.2, 3), 0, 1)), z=float(X[r - 1, 2]), helix=int(rng.random() < 0.7), pool='TM' if abs(X[r - 1, 2]) < 14 else 'TMDx', domain='TMD1', exposed=0) for r in res}
    for r in cov:
        cov[r]['exposed'] = int(cov[r]['rsa'] >= 0.05)
    sets = dict(R_HUB=set(range(5, 25)), GATE_CORE={30, 31, 32, 33}, GATE_ALL=set(range(30, 38)), NBS=set(range(60, 75)))
    hits = {'U1 pooled': set(rng.choice(range(40, 120), 12, replace=False).tolist())}
    hits['U1 Pdr5-own'] = set(list(hits['U1 pooled'])[:6])
    hits['U1 Cdr1-mapped'] = set(list(hits['U1 pooled'])[6:])
    hits['U2 pooled'] = hits['U1 pooled']
    out = run_b1({'7P04': dict(res=res, X=X)}, cov, sets, hits, nd=200, nb=20, cutoffs=(12.0,), n_modes_list=(None,))
    print(f'selftest B1: {len(out)} rows')
    heavy = np.vstack([X + rng.normal(size=X.shape) * 0.8 for _ in range(3)])
    ca = {r: tuple(X[r - 1]) for r in res}
    tpl = rng.normal(size=(8, 3)) * 1.4
    tpl -= tpl.mean(0)
    mols = [tpl @ L.random_rotation(rng).T + X[40] + np.array([6.0, 0, 0]), tpl @ L.random_rotation(rng).T + X[70] + np.array([0, 6.0, 0])]
    sterol = {'9S4S': dict(xyz=np.vstack(mols), mol=np.repeat(np.arange(2), 8))}
    import a12_enm
    o2 = run_b2(X, res, cov, ca, sets, heavy, sterol, nd=20, cut_ps=8.0, cut_ss=8.0, cutoff=12.0)
    print(f'selftest B2: {len(o2)} rows', b2_consistency(o2))
    return (out, o2)

def P(root, *parts):
    return os.path.join(root, *parts)

def build_sterols(root, write_dir=None, entries=INWARD_ENTRIES + OUTWARD_ENTRIES):
    import a12_inputs as I
    ca7, _, _ = I.read_chain(P(root, 'phase0_verify', 'residues', 'downloads', '7P04.cif'), 'A')
    caO, _, _ = I.read_chain(P(root, 'preregistration', 'downloads', 'opm_7p04', 'opm_7p04.pdb'), 'A')
    Ro, to, rms_o = I.to_frame_from_ca(ca7, caO)
    qc, sets = ({}, {})
    for e in entries:
        cif = P(root, 'phase3_structure', 'downloads', 'rcsb_cdr1', f'{e}.cif')
        ca_m, _, _ = I.read_chain(cif, 'A')
        mols = I.read_ligand_molecules(cif, 'ERG')
        trans, pairs = ({}, {})
        for dom in ('D_TMD1', 'D_TMD2'):
            Rm, tm, pr, meta = I.read_usalign_json(P(root, 'phase3_structure', 'outputs', 'alignments', 'pairs_domainwise', f'{e}__{dom}.json'))
            direction, rm = I.check_transform(Rm, tm, ca_m, ca7, pr, max_rmsd=2.0)
            trans[dom] = I.model_to_ref(Rm, tm, direction)
            pairs[dom] = pr
            qc.setdefault(e, {})[dom] = dict(direction=direction, rmsd=round(rm, 3), stored_rmsd=meta['RMSD'])
        xyz, mol, dom = I.sterol_nodes_in_ref_frame(mols, ca_m, trans, pairs, to_opm=(Ro, to))
        qc[e].update(n_mol=len(mols), n_atoms=len(xyz), domains=dom, opm_fit_rmsd=round(rms_o, 4))
        sets[e] = dict(xyz=xyz, mol=mol)
        if write_dir:
            os.makedirs(write_dir, exist_ok=True)
            with open(os.path.join(write_dir, f'a12_sterol_nodes_{e}.tsv'), 'w', encoding='utf-8') as fh:
                fh.write('mol' + chr(9) + 'x' + chr(9) + 'y' + chr(9) + 'z' + chr(10))
                for m, p in zip(mol, xyz):
                    fh.write(f'{m}' + chr(9) + f'{p[0]:.3f}' + chr(9) + f'{p[1]:.3f}' + chr(9) + f'{p[2]:.3f}' + chr(10))
    return (qc, sets)

def load_b_inputs(root):
    import a12_inputs as I
    import a12_run_m as RM
    cov, sets, rows, caO = RM.load_real(root)
    structs = {}
    for tag, path in (('7P04', P(root, 'preregistration', 'downloads', 'opm_7p04', 'opm_7p04.pdb')), ('7P03', P(root, 'phase0_verify', 'residues', 'downloads', '7P03.cif')), ('7P06', P(root, 'phase0_verify', 'residues', 'downloads', '7P06.cif'))):
        ca, _, heavy = I.read_chain(path, 'A')
        res = sorted(ca)
        structs[tag] = dict(res=res, X=np.array([ca[r] for r in res]), heavy=heavy, ca=ca)
    hits = {'U1 pooled': None}
    U1, U2 = (S.hit_sets(rows, 'U1'), S.hit_sets(rows, 'U2'))
    hits = {'U1 pooled': U1['pooled'], 'U1 Pdr5-own': U1['pdr5'], 'U1 Cdr1-mapped': U1['cdr1'], 'U2 pooled': U2['pooled']}
    return (cov, sets, structs, hits)
if __name__ == '__main__':
    if '--selftest' in sys.argv:
        selftest()
    elif '--check-sterols' in sys.argv:
        root = sys.argv[sys.argv.index('--check-sterols') + 1]
        qc, sets = build_sterols(root)
        for e, v in qc.items():
            print(e, v)
    elif '--build-sterols' in sys.argv:
        k = sys.argv.index('--build-sterols')
        root, outdir = (sys.argv[k + 1], sys.argv[k + 2])
        L.guard_amendment_appended(P(root, 'preregistration', 'AMENDMENTS.md'))
        qc, _ = build_sterols(root, outdir)
        print(qc)
    elif '--real-b1' in sys.argv:
        k = sys.argv.index('--real-b1')
        root, outdir = (sys.argv[k + 1], sys.argv[k + 2])
        L.guard_amendment_appended(P(root, 'preregistration', 'AMENDMENTS.md'))
        cov, sets, structs, hits = load_b_inputs(root)
        order = {t: structs[t] for t in ('7P04', '7P03', '7P06')}
        rows = run_b1(order, cov, sets, hits, cutoffs=(15.0, 12.0, 18.0), n_modes_list=(None, 50))
        write_rows(rows, P(outdir, 'a12_B1_results.tsv'))
        print('B1 written', len(rows))
    elif '--real-b2' in sys.argv:
        k = sys.argv.index('--real-b2')
        root, outdir = (sys.argv[k + 1], sys.argv[k + 2])
        L.guard_amendment_appended(P(root, 'preregistration', 'AMENDMENTS.md'))
        cov, sets, structs, hits = load_b_inputs(root)
        qc, ster = build_sterols(root, outdir)
        s = structs['7P04']
        rows = run_b2(s['X'], s['res'], cov, s['ca'], sets, s['heavy'], ster, nd=1000)
        write_rows(rows, P(outdir, 'a12_B2_results.tsv'))
        print('B2 written', len(rows), b2_consistency(rows))
    else:
        print(__doc__)

import argparse
import json
import os
import platform
import sys
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a11_lib as lib
REF_ID_DEFAULT = 'Saccharomyces_cerevisiae|g002312.m1'
WINDOW_HALF = 3
RUNS_REFERENCE_ALIGNER = [('PRIMARY', 'primary', 'primary', 'H6_PRIMARY'), ('S_unweighted', 'unweighted', 'primary', 'H6_PRIMARY'), ('S_nopenalty', 'nopenalty', 'primary', 'H6_PRIMARY'), ('S_window', 'window', 'primary', 'H6_PRIMARY'), ('S_empbg', 'empbg', 'primary', 'H6_PRIMARY'), ('C_comp', 'primary', 'comp', 'H6_PRIMARY'), ('C_exposed', 'primary', 'exposed', 'H6_PRIMARY'), ('T_all', 'primary', 'primary', 'H6_ALL'), ('T_indep', 'primary', 'primary', 'H6_INDEP'), ('T_novel', 'primary', 'primary', 'H6_NOVEL'), ('T_outside_descriptive', 'primary', 'primary', 'H6_OUTSIDE_HUB')]

def residue_scores(seqs, ref_seq, scorable_cols, gap, weights_primary, mode):
    rc = lib.reference_columns(ref_seq)
    cols = {r: c for r, c in rc.items() if c in scorable_cols}
    n = len(seqs)
    if mode == 'unweighted':
        w, bg, pen = ([1.0 / n] * n, lib.BG_BLOSUM62, True)
    elif mode == 'nopenalty':
        w, bg, pen = (weights_primary, lib.BG_BLOSUM62, False)
    elif mode == 'empbg':
        w = weights_primary
        bg, pen = (lib.empirical_background(seqs, w, sorted(scorable_cols)), True)
    else:
        w, bg, pen = (weights_primary, lib.BG_BLOSUM62, True)
    col_scores = lib.js_conservation_columns(seqs, sorted(set(cols.values())), w, bg, pen)
    sc = {r: col_scores[c] for r, c in cols.items()}
    if mode == 'window':
        sc = lib.window_scores(sc, WINDOW_HALF, 0.5)
    return sc

def analyse_alignment(ids, seqs, ref_id, cov, truth_sets, exclude, runs, n_draws=lib.N_DRAWS, seed=lib.SEED_BASE, n_boot=lib.N_BOOT, label=''):
    pos = ids.index(ref_id)
    ref_seq = seqs[pos]
    bad = lib.check_reference(ref_seq, {r: v['aa'] for r, v in cov.items()})
    if bad:
        raise RuntimeError(f'reference row does not match the structure covariates ({len(bad)}): {bad[:5]}')
    gcol = lib.gap_fractions(seqs)
    scorable = {c for c, g in enumerate(gcol) if g <= lib.GAP_NA_THRESHOLD}
    w = lib.henikoff_weights(seqs, sorted(scorable))
    rc = lib.reference_columns(ref_seq)
    gap_res = {r: gcol[rc[r]] for r in cov if r in rc}
    universe = sorted((r for r in cov if r not in exclude))
    cache = {}
    results = {}
    for name, smode, cmode, tset in runs:
        if smode not in cache:
            cache[smode] = residue_scores(seqs, ref_seq, scorable, gcol, w, smode)
        sc = cache[smode]
        res = lib.permutation_test(sorted(truth_sets[tset]), universe, cov, sc, gap=gap_res, n_draws=n_draws, seed=seed, n_boot=n_boot, comp_match=cmode == 'comp', exposed_only=cmode == 'exposed')
        res.update({'run': name, 'score_mode': smode, 'control_mode': cmode, 'truth_set': tset, 'alignment': label, 'n_rows': len(seqs)})
        results[name] = res
    info = {'n_rows': len(seqs), 'n_columns': len(gcol), 'n_scorable_columns': len(scorable), 'ref_residues_scorable': sum((1 for r in rc if rc[r] in scorable))}
    return (results, info)

def conf_variant(ids, seqs, ref_id, cov, truth, frozen, conf_res, n_draws=lib.N_DRAWS, n_boot=lib.N_BOOT, label=''):
    conf = {r for r in conf_res if r in cov}
    t2 = dict(truth)
    t2['H6_UNRELEASED'] = conf
    base_excl = set(frozen) | set(truth['H6_ALL']) | conf
    r1, info = analyse_alignment(ids, seqs, ref_id, cov, t2, base_excl, [('H6_UNRELEASED', 'primary', 'primary', 'H6_UNRELEASED')], n_draws, lib.SEED_BASE, n_boot, label)
    r2, _ = analyse_alignment(ids, seqs, ref_id, cov, t2, base_excl, [('H6_exclUNRELEASED', 'primary', 'primary', 'H6_PRIMARY')], n_draws, lib.SEED_BASE, n_boot, label)
    r1.update(r2)
    return (r1, info)

def write_results(out_dir, all_results, fname='A11_H6_results.tsv'):
    cols = ['alignment', 'run', 'truth_set', 'score_mode', 'control_mode', 'n_rows', 'n_requested', 'n_eff', 'status', 'observed_mean', 'control_mean', 'control_sd', 'delta', 'z', 'ci95_lo', 'ci95_hi', 'b', 'p', 'n_fallback_residues', 'weakly_matched', 'mde_80pct_power_alpha05_onesided', 'dropped_NA']
    with open(os.path.join(out_dir, fname), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\t'.join(cols) + '\n')
        for lab, rs in all_results.items():
            for name, r in rs.items():
                row = dict(r)
                ci = row.get('ci95_delta', ('', ''))
                row['ci95_lo'], row['ci95_hi'] = ci
                row['dropped_NA'] = ','.join((str(x) for x in row.get('dropped_NA', [])))
                fh.write('\t'.join((str(row.get(c, '')) for c in cols)) + '\n')

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', required=True, help='project directory')
    ap.add_argument('--aln-mafft', required=True, help='full untrimmed MAFFT alignment (Windows copy)')
    ap.add_argument('--aln-muscle', required=True, help='full untrimmed MUSCLE super5 alignment (Windows copy)')
    ap.add_argument('--expect-sha-mafft', default='77e7a5fddd2365af9854bbaaa25f21ced1a906be97eff3871bd82c4e60d63d25', help='sha256 recorded in mafft.aln_qc.json (PHASE1B)')
    ap.add_argument('--expect-sha-muscle', default='841f438c9f754c9c1d87bc203e50c9d5b6162cef1b0fe4678db1a295ee142ab3', help='sha256 recorded in muscle.aln_qc.json (PHASE1B)')
    ap.add_argument('--truth-sets', required=True, help='a11_truth_sets.tsv from a11_build_truth_sets.py')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--ref-id', default=REF_ID_DEFAULT)
    ap.add_argument('--n-draws', type=int, default=lib.N_DRAWS)
    ap.add_argument('--check-only', action='store_true', help='integrity checks only; no score is computed')
    ap.add_argument('--conf-set', default=None, help="optional residue list (one number per line)")
    a = ap.parse_args()
    R = a.root
    p_man = os.path.join(R, 'phase1_orthology', 'phase1b', 'outputs', 'set', 'analysis_set_manifest.tsv')
    p_sub = os.path.join(R, 'phase1_orthology', 'phase1c', 'outputs', 'subset_manifest.tsv')
    p_cov = os.path.join(R, 'preregistration', 'outputs', 'structure_covariates_7P04.tsv')
    p_frz = os.path.join(R, 'preregistration', 'site_classes_FROZEN.tsv')
    expect = {'analysis_set_manifest.tsv': '735b5231a4b9cdde5372dd16f679d96f1664c07269897e8c6e7d166f920c4f85', 'subset_manifest.tsv': 'eaf72750d4ae21052471bedb65b72b7fa4c3c6093e7b8ac1e81794e2fbcce245', 'site_classes_FROZEN.tsv': 'e409356b6e46c1b0b10fdaa79a1f1b1c835cfbe038b939ee03ceeb8d853fb0a3'}
    for p in (p_man, p_sub, p_frz):
        got = lib.sha256_file(p)
        if got != expect[os.path.basename(p)]:
            sys.exit(f'HASH MISMATCH {p}: {got}')
    os.makedirs(a.out_dir, exist_ok=True)
    cov = lib.load_covariates(p_cov)
    frozen = {int(r['residue']) for r in lib.read_tsv(p_frz)}
    truth = {}
    for r in lib.read_tsv(a.truth_sets):
        truth.setdefault(r['set'], set()).add(int(r['residue']))
    exclude = frozen | truth['H6_ALL']
    man_all = [r['seq_id'] for r in lib.read_tsv(p_man) if r['role'] == 'INGROUP']
    man_sub = [r['seq_id'] for r in lib.read_tsv(p_sub) if r['role'] == 'INGROUP']
    if len(man_all) != 1253 or len(man_sub) != 530:
        sys.exit(f'unexpected row counts: ingroup {len(man_all)} (expect 1253), subset ingroup {len(man_sub)} (expect 530)')
    manifest = {'python': sys.version, 'platform': platform.platform(), 'start_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'seed': lib.SEED_BASE, 'n_draws': a.n_draws, 'inputs': {}}
    for tag, path, pref in (('mafft', a.aln_mafft, a.expect_sha_mafft), ('muscle', a.aln_muscle, a.expect_sha_muscle)):
        h = lib.sha256_file(path)
        manifest['inputs'][tag] = {'path': path, 'sha256': h}
        if pref and (not h.startswith(pref)):
            sys.exit(f'alignment {tag}: sha256 {h} does not start with the recorded prefix {pref}')
    for p in (p_man, p_sub, p_cov, p_frz, a.truth_sets):
        manifest['inputs'][os.path.basename(p)] = {'path': p, 'sha256': lib.sha256_file(p)}
    if a.conf_set:
        import a11_h7_overlap as h7mod
        h7mod.guard(a.out_dir, False)
        with open(a.conf_set, 'r', encoding='utf-8') as fh:
            B = [int(x.split()[0]) for x in fh if x.strip() and (not x.startswith('#'))]
        conf_all = {}
        for tag, path in (('mafft', a.aln_mafft), ('muscle', a.aln_muscle)):
            ids, seqs = lib.read_fasta(path)
            sel_ids, sel_seqs = lib.select_rows(ids, seqs, man_all)
            conf_all[f'{tag}_ingroup'], _ = conf_variant(sel_ids, sel_seqs, a.ref_id, cov, truth, frozen, B, n_draws=a.n_draws, label=f'{tag}_ingroup')
        write_results(a.out_dir, conf_all, 'A11_H6_UNRELEASED_results_CONF.tsv')
        print('conditional variant written')
        return
    all_results = {}
    for tag, path in (('mafft', a.aln_mafft), ('muscle', a.aln_muscle)):
        ids, seqs = lib.read_fasta(path)
        for rowset, want in (('ingroup', man_all), ('sacch', man_sub)):
            sel_ids, sel_seqs = lib.select_rows(ids, seqs, want)
            lab = f'{tag}_{rowset}'
            if a.check_only:
                bad = lib.check_reference(sel_seqs[sel_ids.index(a.ref_id)], {r: v['aa'] for r, v in cov.items()})
                print(f'{lab}: {len(sel_ids)} rows, {len(sel_seqs[0])} columns, reference mismatches: {len(bad)}')
                if bad:
                    sys.exit('reference mismatch: ' + '; '.join(bad[:5]))
                continue
            runs = RUNS_REFERENCE_ALIGNER if rowset == 'ingroup' else [RUNS_REFERENCE_ALIGNER[0]]
            res, info = analyse_alignment(sel_ids, sel_seqs, a.ref_id, cov, truth, exclude, runs, n_draws=a.n_draws, label=lab)
            all_results[lab] = res
            manifest.setdefault('alignment_info', {})[lab] = info
            print(lab, info, {k: (v.get('p'), v.get('delta')) for k, v in res.items()})
    if a.check_only:
        return
    write_results(a.out_dir, all_results)
    M, U = (all_results['mafft_ingroup'], all_results['muscle_ingroup'])
    runs = {'R0': M['PRIMARY'], 'R1': U['PRIMARY'], 'R2': all_results['mafft_sacch']['PRIMARY'], 'R3': all_results['muscle_sacch']['PRIMARY']}
    for k in ('S_unweighted', 'S_nopenalty', 'S_window', 'S_empbg', 'C_comp', 'C_exposed', 'T_indep'):
        runs[k] = M[k]
    lab = lib.robust_label_h6(runs)
    with open(os.path.join(a.out_dir, 'A11_H6_label.json'), 'w', encoding='utf-8') as fh:
        json.dump(lab, fh, indent=2)
    manifest['end_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    with open(os.path.join(a.out_dir, 'A11_H6_manifest.json'), 'w', encoding='utf-8') as fh:
        json.dump(manifest, fh, indent=2)
    print('H6 label:', lab['label'], 'failed checks:', lab['failed'])
if __name__ == '__main__':
    main()

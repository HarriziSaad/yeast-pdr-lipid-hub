import csv
import os
import sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import a12_lib as L
import a12_m as M
import a12_select as S
NARROW = ((0.05, 2.0), (0.075, 3.0), (0.1, 4.0))
FIELDS = ['test', 'variant', 'region', 'hits', 'n_hit_residues', 'n_eff', 'obs', 'null_mean', 'delta', 'fold', 'z', 'p', 'q', 'ci_lo', 'ci_hi', 'mde80', 'frac_fallback', 'weakly_matched', 'smd_rsa', 'smd_depth', 'balance_flag', 'note']

def _row(test, variant, region, hits_name, n_hits, res, note=''):
    r = {k: '' for k in FIELDS}
    r.update(test=test, variant=variant, region=region, hits=hits_name, n_hit_residues=n_hits, note=note)
    r.update(res)
    return r

def _m1a_row(variant, region_name, region, hits_name, hits, cov, nd, nb, tiers=None, tested=None, seed=0):
    o = M.m1a(region, hits, cov, n_draws=nd, n_boot=nb, tiers=tiers, tested=tested, seed_offset=seed)
    ci = o['ci95_delta_fraction']
    return _row('M1a', variant, region_name, hits_name, len(hits), dict(n_eff=o['n_eff'], obs=o['obs'], null_mean=o['null_mean'], delta=o['delta'], fold=o['fold'], z=o['z'], p=o['p'], ci_lo=ci[0] * o['n_eff'], ci_hi=ci[1] * o['n_eff'], mde80=o['mde80'], frac_fallback=o['match']['frac_fallback'], weakly_matched=int(o['match']['weakly_matched']), smd_rsa=o['balance'].get('smd_rsa', ''), smd_depth=o['balance'].get('smd_depth', ''), balance_flag=int(o['balance'].get('imbalanced', False))))

def _m1b_row(variant, region_name, region, cav, hits_name, hits, cov, nperm, nb, tiers=None, seed=0):
    o = M.m1b(region, cav, hits, cov, n_perm=nperm, n_boot=nb, tiers=tiers, seed_offset=seed)
    return _row('M1b', variant, region_name, hits_name, len(hits), dict(n_eff=f"{o['n_R']}/{o['n_C']}", obs=o['excess_R'], null_mean=o['excess_C'], delta=o['delta_RC'], z=o['z'], p=o['p'], ci_lo=o['ci95_delta'][0], ci_hi=o['ci95_delta'][1], mde80=o['mde80'], frac_fallback=max(o['matchR']['frac_fallback'], o['matchC']['frac_fallback']), weakly_matched=int(o['matchR']['weakly_matched'] or o['matchC']['weakly_matched'])), note='obs=excess of region, null_mean=excess of cavity (adjusted hit fractions)')

def _m2_row(variant, hits_name, hits, ca_xyz, cov, nd, radius=10.0, seed=0, note=''):
    o = M.m2(hits, ca_xyz, cov, n_draws=nd, radii=(8.0, 10.0, 12.0, 15.0), primary_r=radius, seed_offset=seed)
    pr = o['radii'][radius]
    return _row('M2', variant, '-', hits_name, len(hits), dict(n_eff=o['n_eff'], obs=pr['obs'], null_mean=pr['null_mean'], delta=pr['delta'], z=pr['z'], p=pr['p'], mde80=pr['mde80'], frac_fallback=o['match']['frac_fallback'], weakly_matched=int(o['match']['weakly_matched'])), note=note + ' r8 p=%.4f r12 p=%.4f r15 p=%.4f meanNN p=%.4f' % (o['radii'][8.0]['p'], o['radii'][12.0]['p'], o['radii'][15.0]['p'], o['mean_nn']['p']))

def run_m(cov, sets, rows, ca_xyz, n_draws=10000, n_boot=2000, n_perm=10000, log=print):
    out = []
    U1, U2, UA = (S.hit_sets(rows, t) for t in ('U1', 'U2', 'UANY'))
    TAS = set(S.residue_table(rows, 'TAS'))
    TESTED = set(S.residue_table(rows, 'TESTED')) | TAS
    R, CAV = (sets['R_HUB'], sets['CAV'])
    out.append(_m1a_row('PRIMARY', 'R_HUB', R, 'U1 pooled', U1['pooled'], cov, n_draws, n_boot))
    out.append(_m1b_row('PRIMARY', 'R_HUB', R, CAV, 'U1 pooled', U1['pooled'], cov, n_perm, n_boot))
    out.append(_m2_row('PRIMARY', 'U1 pooled', U1['pooled'], ca_xyz, cov, n_draws))
    for name in ('R_HUBGATE', 'R_LHCONF', 'R_H6'):
        out.append(_m1a_row('region', name, sets[name], 'U1 pooled', U1['pooled'], cov, n_draws, n_boot))
    for hn, hs in (('U1 Pdr5-own', U1['pdr5']), ('U1 Cdr1-mapped', U1['cdr1']), ('U2 pooled', U2['pooled']), ('UANY pooled', UA['pooled'])):
        out.append(_m1a_row('hits', 'R_HUB', R, hn, hs, cov, n_draws, n_boot))
    for p, hs in U1['loo'].items():
        out.append(_m1a_row('leave-one-paper-out', 'R_HUB', R, f'U1 minus {p}', hs, cov, n_draws, n_boot))
    out.append(_m1a_row('narrow calipers', 'R_HUB', R, 'U1 pooled', U1['pooled'], cov, n_draws, n_boot, tiers=NARROW))
    lg = M.m1a_logistic(R, U1['pooled'], cov)
    out.append(_row('M1a', 'logistic-adjusted (exact Poisson-binomial)', 'R_HUB', 'U1 pooled', len(U1['pooled']), dict(n_eff=lg['n_eff'], obs=lg['obs'], null_mean=lg['expected'], delta=lg['obs'] - lg['expected'], p=lg['p'])))
    out.append(_m1a_row('targeted mutants, tested-residue universe (BIASED design, sensitivity)', 'R_HUB', R, 'TAS (targeted AS)', TAS, cov, n_draws, n_boot, tested=TESTED))
    for name in ('R_HUBGATE', 'R_LHCONF'):
        out.append(_m1b_row('region', name, sets[name], CAV, 'U1 pooled', U1['pooled'], cov, n_perm, n_boot))
    for hn, hs in (('U1 Pdr5-own', U1['pdr5']), ('U1 Cdr1-mapped', U1['cdr1']), ('U2 pooled', U2['pooled'])):
        out.append(_m1b_row('hits', 'R_HUB', R, CAV, hn, hs, cov, n_perm, n_boot))
    for p, hs in U1['loo'].items():
        out.append(_m1b_row('leave-one-paper-out', 'R_HUB', R, CAV, f'U1 minus {p}', hs, cov, n_perm, n_boot))
    out.append(_m1b_row('narrow calipers', 'R_HUB', R, CAV, 'U1 pooled', U1['pooled'], cov, n_perm, n_boot, tiers=NARROW))
    for hn, hs in (('U1 Pdr5-own', U1['pdr5']), ('U1 Cdr1-mapped', U1['cdr1']), ('U2 pooled', U2['pooled'])):
        out.append(_m2_row('hits', hn, hs, ca_xyz, cov, n_draws))
    out.append(_m2_row('hits outside R_HUB', 'U1 pooled minus R_HUB', U1['pooled'] - set(R), ca_xyz, cov, n_draws))
    for p, hs in U1['loo'].items():
        out.append(_m2_row('leave-one-paper-out', f'U1 minus {p}', hs, ca_xyz, cov, n_draws))
    prim = [r for r in out if r['variant'] == 'PRIMARY']
    fm1 = [r for r in prim if r['test'] in ('M1a', 'M1b')]
    q = L.bh_qvalues([r['p'] for r in fm1])
    for r, qq in zip(fm1, q):
        r['q'] = qq
    for r in prim:
        if r['test'] == 'M2':
            r['q'] = r['p']
    return out

def _f(x):
    return float(x) if x not in ('', None) else float('nan')

def consistency_labels(rows):
    labels = {}
    for test in ('M1a', 'M1b'):
        prim = next((r for r in rows if r['test'] == test and r['variant'] == 'PRIMARY'))
        fails = []
        if not _f(prim['q']) < 0.1:
            labels[test] = 'no evidence (q >= 0.1 at the primary)'
            continue
        sub = [r for r in rows if r['test'] == test and r['variant'] != 'PRIMARY']
        for r in sub:
            if r['variant'] in ('hits', 'leave-one-paper-out', 'region', 'narrow calipers') and (not _f(r['delta']) > 0) and (not r['hits'].startswith('U2')) and (not r['hits'].startswith('UANY')):
                fails.append(f"sign: {r['variant']} / {r['hits']} / {r['region']}")
        if test == 'M1a':
            lg = next((r for r in sub if r['variant'].startswith('logistic')))
            if not _f(lg['p']) < 0.1:
                fails.append('logistic-adjusted p >= 0.1')
            if prim['balance_flag'] == 1:
                fails.append('caliper imbalance flagged')
        labels[test] = 'consistent' if not fails else 'not consistent: ' + '; '.join(fails)
    m2 = next((r for r in rows if r['test'] == 'M2' and r['variant'] == 'PRIMARY'))
    if not _f(m2['q']) < 0.1:
        labels['M2'] = 'no evidence (q >= 0.1 at the primary)'
    else:
        fails = []
        for r in rows:
            if r['test'] == 'M2' and r['variant'] == 'hits' and (r['hits'] in ('U1 Pdr5-own', 'U1 Cdr1-mapped')) and (not _f(r['delta']) > 0):
                fails.append(f"sign: {r['hits']}")
        labels['M2'] = 'consistent' if not fails else 'not consistent: ' + '; '.join(fails)
    return labels

def write_results(rows, labels, outdir):
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'a12_M_results.tsv'), 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(outdir, 'a12_M_labels.txt'), 'w', encoding='utf-8') as fh:
        for k, v in labels.items():
            fh.write(f'{k}\t{v}\n')

def selftest():
    import test_a12_toy as T
    cov = T.synth_cov(700, seed=3)
    ca = T.synth_coords(cov, seed=4)
    tm = [r for r in cov if cov[r]['pool'] == 'TM']
    sets = dict(R_HUB=set(tm[10:46]), R_HUBGATE=set(tm[10:56]), R_LHCONF=set(tm[5:60]), R_H6=set(tm[10:37]), CAV=set(tm[60:78]) | {tm[20]})
    rows = []
    rng = np.random.default_rng(5)
    for k, (paper, prot) in enumerate([('tanabe2019', 'Pdr5'), ('tanabe2019', 'Cdr1'), ('niimi2022', 'Cdr1'), ('egner1998', 'Pdr5'), ('niimi2012', 'Cdr1')]):
        for r in rng.choice(list(cov), 14, replace=False):
            rows.append(dict(paper_key=paper, protein=prot, residue_native=str(int(r)), pdr5_equivalent_residue=str(int(r)), evidence_use='single', phenotype_category='altered specificity', expr_loc_defect_flag='N', evidence_access='full text', other_flags='', residue_verified_vs_UniProt='Y (A)', phenotype_summary='Spontaneous isolate'))
    out = run_m(cov, sets, rows, ca, n_draws=200, n_boot=30, n_perm=200)
    labels = consistency_labels(out)
    print(f'selftest: {len(out)} result rows; labels: {labels}')
    return (out, labels)
EXPECTED_HASHES = {'preregistration/outputs/structure_covariates_7P04.tsv': '5a32e3fd32c93aee4bd0093da69caea657a3fa2f0aa2345c4ede52d0c9644eba', 'preregistration/site_classes_FROZEN.tsv': 'e409356b6e46c1b0b10fdaa79a1f1b1c835cfbe038b939ee03ceeb8d853fb0a3', 'reframing/outputs/a11_truth_sets.tsv': 'a2e1781c76792e5aa1454533334f46dc5288a8ec14d7052ae6a57080a77eae1c', 'phase5_validation/mutants/mutant_table.tsv': '02666b17e66bb3b9d1c733ce1388060d2c8ac16e8abd30d1ae4650932d1c942f', 'preregistration/downloads/opm_7p04/opm_7p04.pdb': '6080fbc9f328cd75652b4a31a29064ef8718d8e654fe2807a182f40454084f52', 'phase5_validation/A12/outputs/a12_sets.tsv': '6ec1056716615cc959959fb6c1bd381df774302c89c1dc64069c528d2ad84e1a'}

def load_real(root):
    import a12_sets as ST
    for rel, h in EXPECTED_HASHES.items():
        got = L.sha256_file(os.path.join(root, *rel.split('/')))
        if got != h:
            raise SystemExit(f'hash mismatch for {rel}: {got}')
    cov = L.read_covariates(os.path.join(root, 'preregistration', 'outputs', 'structure_covariates_7P04.tsv'))
    sets = ST.build_sets(ST.read_classes(os.path.join(root, 'preregistration', 'site_classes_FROZEN.tsv')), ST.read_truth(os.path.join(root, 'reframing', 'outputs', 'a11_truth_sets.tsv')))
    stored = L.read_set_table(os.path.join(root, 'phase5_validation', 'A12', 'outputs', 'a12_sets.tsv'))
    if any((stored[k] != sets[k] for k in sets)):
        raise SystemExit('region sets differ from the stored a12_sets.tsv')
    rows = S.read_table(os.path.join(root, 'phase5_validation', 'mutants', 'mutant_table.tsv'))
    import a12_inputs as I
    ca, _, _ = I.read_chain(os.path.join(root, 'preregistration', 'downloads', 'opm_7p04', 'opm_7p04.pdb'), 'A')
    return (cov, sets, rows, ca)
if __name__ == '__main__':
    if '--selftest' in sys.argv:
        selftest()
    elif '--check-only' in sys.argv:
        root = sys.argv[sys.argv.index('--check-only') + 1]
        cov, sets, rows, ca = load_real(root)
        print(f'inputs verified: {len(cov)} covariate rows, {len(rows)} mutant rows, {len(ca)} CA atoms, sets {sorted(((k, len(v)) for k, v in sets.items()))}')
    elif '--real' in sys.argv:
        k = sys.argv.index('--real')
        root, outdir = (sys.argv[k + 1], sys.argv[k + 2])
        L.guard_amendment_appended(os.path.join(root, 'preregistration', 'AMENDMENTS.md'))
        cov, sets, rows, ca = load_real(root)
        res = run_m(cov, sets, rows, ca)
        labels = consistency_labels(res)
        write_results(res, labels, outdir)
        print('M block written to', outdir, labels)
    else:
        print(__doc__)

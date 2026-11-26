import argparse
import collections
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a11_lib as lib
HUB_HELICES = ('TH10', 'TH11a', 'TH11b', 'TH11c', 'TH12')
NAMED4 = {1329: 'L', 1478: 'Y', 1365: 'M', 1383: 'W'}

def harris_helix(res: int, ranges):
    for name, a, b in ranges:
        if a <= res <= b:
            return name
    return 'none'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', required=True, help='project directory (contains preregistration, phase3_structure, ...)')
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    R = a.root
    p_frozen = os.path.join(R, 'preregistration', 'site_classes_FROZEN.tsv')
    p_proxy = os.path.join(R, 'phase3_structure', 'outputs', 'c4b_proxy_LH-erg-CgCdr1proxy.tsv')
    p_harris = os.path.join(R, 'phase0_verify', 'residues', 'scripts', 'harris_th_helices.tsv')
    p_cov = os.path.join(R, 'preregistration', 'outputs', 'structure_covariates_7P04.tsv')
    p_pool = os.path.join(R, 'preregistration', 'outputs', 'g3_calibration_pool_RSA005.tsv')
    os.makedirs(a.out_dir, exist_ok=True)
    cov = lib.load_covariates(p_cov)
    harris = [(r['name'], int(r['start']), int(r['end'])) for r in lib.read_tsv(p_harris)]
    frozen = collections.defaultdict(set)
    frozen_aa = {}
    for r in lib.read_tsv(p_frozen):
        frozen[int(r['residue'])].add(r['class'])
        frozen_aa[int(r['residue'])] = r['aa']
    proxy = {}
    for r in lib.read_tsv(p_proxy):
        proxy[int(r['P33302_residue'])] = (r['P33302_aa'], r['Harris_helix'])
    pool = {int(r['residue']) for r in lib.read_tsv(p_pool)}
    problems = []
    for res, (aa, hx) in proxy.items():
        if cov[res]['aa'] != aa:
            problems.append(f"proxy {res}: covariates aa {cov[res]['aa']} != proxy aa {aa}")
        if harris_helix(res, harris) != hx:
            problems.append(f'proxy {res}: Harris table gives {harris_helix(res, harris)}, proxy file {hx}')
    for res, aa in NAMED4.items():
        if cov[res]['aa'] != aa or frozen_aa.get(res) != aa:
            problems.append(f"named residue {res}{aa}: covariates {cov[res]['aa']}, frozen {frozen_aa.get(res)}")
        if not {'LH-ergC', 'LH-ergH'} & frozen[res]:
            problems.append(f'named residue {res}: not in LH-ergC/LH-ergH of the frozen table')
        if harris_helix(res, harris) not in HUB_HELICES:
            problems.append(f'named residue {res}: helix {harris_helix(res, harris)} is not a hub helix')
    if problems:
        print('CONSISTENCY PROBLEMS:\n  ' + '\n  '.join(problems))
        sys.exit(2)
    named = set(NAMED4)
    proxy_all = set(proxy)
    proxy_hub = {r for r, (_, hx) in proxy.items() if hx in HUB_HELICES}
    sets = {'H6_PRIMARY': named | proxy_hub, 'H6_ALL': named | proxy_all, 'H6_INDEP': proxy_hub - named, 'H6_NOVEL': (named | proxy_hub) - set(frozen), 'H6_OUTSIDE_HUB': proxy_all - proxy_hub - named}
    excl = set(frozen) | sets['H6_ALL']
    universe = sorted((r for r in cov if r not in excl))
    rows = []
    for name, s in sets.items():
        for r in sorted(s):
            rows.append({'set': name, 'residue': r, 'aa': cov[r]['aa'], 'harris_helix': harris_helix(r, harris), 'in_named4': int(r in named), 'in_cgcdr1_proxy': int(r in proxy_all), 'frozen_classes': ';'.join(sorted(frozen[r])) if r in frozen else '-', 'pool_label': cov[r]['pool'], 'rsa': f"{cov[r]['rsa']:.3f}", 'z_OPM': f"{cov[r]['z']:.2f}", 'in_helx_p': cov[r]['helix'], 'in_g3_pool': int(r in pool)})
    cols = list(rows[0])
    with open(os.path.join(a.out_dir, 'a11_truth_sets.tsv'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join((str(r[c]) for c in cols)) + '\n')
    L = []
    L.append('A11 truth sets (structure-only; no alignment, conservation, population or UNRELEASED input was read)')
    L.append(f'inputs: frozen table {lib.sha256_file(p_frozen)[:16]}..., proxy file {lib.sha256_file(p_proxy)}, harris table {lib.sha256_file(p_harris)[:16]}..., covariates {lib.sha256_file(p_cov)[:16]}..., G3 pool {lib.sha256_file(p_pool)[:16]}...')
    L.append('')
    for name, s in sets.items():
        pl = collections.Counter((cov[r]['pool'] for r in s))
        hx = collections.Counter((harris_helix(r, harris) for r in s))
        shared = sorted((r for r in s if r in frozen))
        L.append(f'{name}: n = {len(s)}; in frozen table: {len(shared)}; in G3 pool: {sum((1 for r in s if r in pool))}; pool labels {dict(pl)}; Harris helices {dict(sorted(hx.items()))}')
        L.append(f"""   residues: {' '.join((f"{cov[r]['aa']}{r}" for r in sorted(s)))}""")
        L.append(f"""   shared with frozen classes: {' '.join((f"{cov[r]['aa']}{r}[{','.join(sorted(frozen[r]))}]" for r in shared))}""")
    L.append('')
    L.append(f"control universe (all 7P04 modelled residues minus the 273 frozen-table residues minus H6_ALL): {len(universe)} residues; by pool {dict(collections.Counter((cov[r]['pool'] for r in universe)))}")
    cands, level = lib.candidate_sets(sorted(sets['H6_PRIMARY']), universe, cov, gap=None)
    lv = collections.Counter(level.values())
    L.append(f'structure-only matching feasibility for H6_PRIMARY (gap caliper ignored): levels {dict(lv)}; candidates per residue min {min((len(v) for v in cands.values()))}, median {sorted((len(v) for v in cands.values()))[len(cands) // 2]}, max {max((len(v) for v in cands.values()))}')
    L.append('residues needing fallback (structure-only): ' + (' '.join((f"{cov[r]['aa']}{r}" for r, v in level.items() if v.startswith('fallback'))) or 'none'))
    pa = len(proxy_all & pool)
    L.append(f'set sizes for H7 (public side only): proxy n = {len(proxy_all)}, proxy within the G3 pool (N = {len(pool)}) a_pool = {pa}')
    with open(os.path.join(a.out_dir, 'a11_truth_sets_summary.txt'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(L) + '\n')
    print('\n'.join(L))
if __name__ == '__main__':
    main()

import sys, os, json, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import p3lib as L
info = json.load(open(os.path.join(L.OUT, 'structure_info.json')))
cov = {int(r['residue']): r['domain'] for r in L.load_tsv(os.path.join(L.IN, 'structure_covariates_7P04.tsv'))}
grp = {'NBD1': 'NBD1', 'NBD2': 'NBD2', 'TMD1': 'TMD1', 'ECD1': 'TMD1', 'TMD2': 'TMD2', 'ECD2': 'TMD2'}
MODELS = {'Pdr5': 'P33302', 'Snq2': 'P32568', 'Pdr10': 'P51533', 'Pdr11': 'P40550', 'Pdr12': 'Q02785', 'Pdr15': 'Q04182', 'Pdr18': 'P53756', 'Aus1': 'Q08409', 'YOL075C': 'Q08234', 'Adp1': 'P25371', 'CaCdr1_AF': 'Q5ANA3'}
meta = {m['uniprotAccession']: m for m in json.load(open(os.path.join(L.IN, 'alphafold_meta.json')))}
rows, pae_rows = ([], [])
for name, acc in MODELS.items():
    st, _ = L.read_clean(os.path.join(L.OUT, 'pdb', f'{name}.pdb'))
    rl = [r for r in L.residues(st) if r['ca'] is not None]
    pl = np.array([r['b'] for r in rl])
    nums = [r['num'] for r in rl]
    a4 = json.load(open(os.path.join(L.OUT, 'align', f'{name}__to__7P04.json')))
    pair = {p[0]: p[1] for p in a4['pairs']}
    aligned = [n for n in nums if n in pair]
    core = np.array([pl[nums.index(n)] for n in aligned])
    rows.append({'model': name, 'UniProt': acc, 'evidence': 'C-AlphaFold', 'AFDB_version': meta[acc]['latestVersion'], 'model_date': meta[acc]['modelCreatedDate'][:10], 'n_res': len(rl), 'AFDB_global_pLDDT': meta[acc]['globalMetricValue'], 'mean_pLDDT_file': pl.mean(), 'frac_gt90': (pl > 90).mean(), 'frac_70_90': ((pl > 70) & (pl <= 90)).mean(), 'frac_50_70': ((pl > 50) & (pl <= 70)).mean(), 'frac_le50': (pl <= 50).mean(), 'n_aligned_to_7P04': len(aligned), 'mean_pLDDT_aligned_core': core.mean(), 'frac_core_ge80': (core >= 80).mean(), 'TM_7P04_normRef': a4['TM_norm_ref'], 'RMSD_7P04_A': a4['RMSD'], 'first_aligned_7P04_res': min(pair.values()), 'last_aligned_7P04_res': max(pair.values())})
    pj = json.load(open(os.path.join(L.IN, 'alphafold', f'{acc}_pae.json')))
    if isinstance(pj, list):
        pj = pj[0]
    pae = np.array(pj['predicted_aligned_error'], dtype=float)
    assert pae.shape[0] == len(rl), (name, pae.shape, len(rl))
    lab = {}
    for n in aligned:
        g = grp.get(cov.get(pair[n]))
        if g:
            lab[n] = g
    srt = sorted(lab)
    for i, n in enumerate(nums):
        if n in lab:
            continue
        prev = max([m for m in srt if m < n], default=None)
        nxt = min([m for m in srt if m > n], default=None)
        if prev is not None and nxt is not None and (lab[prev] == lab[nxt]) and (n - prev <= 15) and (nxt - n <= 15):
            lab[n] = lab[prev]
    idx = {g: np.array([i for i, n in enumerate(nums) if lab.get(n) == g], dtype=int) for g in ('NBD1', 'NBD2', 'TMD1', 'TMD2')}
    rec = {'model': name}
    for a in idx:
        sub = pae[np.ix_(idx[a], idx[a])]
        rec[f'within_{a}'] = sub.mean() if len(idx[a]) else np.nan
        rec[f'n_{a}'] = len(idx[a])
    for a, b in (('NBD1', 'NBD2'), ('NBD1', 'TMD1'), ('NBD1', 'TMD2'), ('NBD2', 'TMD1'), ('NBD2', 'TMD2'), ('TMD1', 'TMD2')):
        if len(idx[a]) and len(idx[b]):
            rec[f'{a}-{b}'] = 0.5 * (pae[np.ix_(idx[a], idx[b])].mean() + pae[np.ix_(idx[b], idx[a])].mean())
        else:
            rec[f'{a}-{b}'] = np.nan
    pae_rows.append(rec)
    print(name, 'meanpLDDT %.2f (API %s) core %.2f' % (pl.mean(), meta[acc]['globalMetricValue'], core.mean()), {k: round(v, 1) for k, v in rec.items() if '-' in k}, flush=True)
with open(os.path.join(L.OUT, 'model_qc.tsv'), 'w') as f:
    cols = list(rows[0])
    f.write('\t'.join(cols) + '\n')
    for r in rows:
        f.write('\t'.join(('%.3f' % r[c] if isinstance(r[c], (float, np.floating)) else str(r[c]) for c in cols)) + '\n')
with open(os.path.join(L.OUT, 'pae_domains.tsv'), 'w') as f:
    cols = list(pae_rows[0])
    f.write('\t'.join(cols) + '\n')
    for r in pae_rows:
        f.write('\t'.join(('%.2f' % r[c] if isinstance(r[c], (float, np.floating)) else str(r[c]) for c in cols)) + '\n')

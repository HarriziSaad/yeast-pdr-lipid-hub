import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import p3lib as L
info = json.load(open(os.path.join(L.OUT, 'structure_info.json')))
cov = {int(r['residue']): r for r in L.load_tsv(os.path.join(L.IN, 'structure_covariates_7P04.tsv'))}

def coords(name):
    st, _ = L.read_clean(os.path.join(L.OUT, 'pdb', f'{name}.pdb'))
    return {r['num']: r['ca'] for r in L.residues(st) if r['ca'] is not None}
sets = {'NBD1': [], 'NBD2': [], 'cyt1': [], 'cyt2': [], 'ext1': [], 'ext2': []}
for n, r in cov.items():
    d, z = (r['domain'], float(r['z_CA_OPM_A']))
    if d == 'NBD1':
        sets['NBD1'].append(n)
    elif d == 'NBD2':
        sets['NBD2'].append(n)
    if r['TM_helix'] == '1':
        tmd = '1' if d in ('TMD1', 'ECD1') else '2' if d in ('TMD2', 'ECD2') else None
        if tmd:
            if -15.2 <= z <= -5.0:
                sets['cyt' + tmd].append(n)
            elif 5.0 <= z <= 15.2:
                sets['ext' + tmd].append(n)
print({k: len(v) for k, v in sets.items()})
REFS = ['7P03', '7P04', '7P05', '7P06']
refxyz = {r: coords(r) for r in REFS}

def centroid(xyz, nums):
    pts = [xyz[n] for n in nums if n in xyz]
    return (np.mean(pts, axis=0), len(pts))

def descriptors(xyz, allowed):
    out = {}
    c = {}
    cnt = {}
    for k, v in sets.items():
        c[k], cnt[k] = centroid(xyz, [n for n in v if n in allowed])
    out['d_NBD'] = float(np.linalg.norm(c['NBD1'] - c['NBD2']))
    out['g_cyt'] = float(np.linalg.norm(c['cyt1'] - c['cyt2']))
    out['g_ext'] = float(np.linalg.norm(c['ext1'] - c['ext2']))
    gd = [np.linalg.norm(xyz[a] - xyz[b]) for a, b in ((683, 1373), (682, 1373), (684, 1373)) if a in allowed and b in allowed and (a in xyz) and (b in xyz)]
    out['d_gate'] = float(np.mean(gd)) if len(gd) == 3 else float('nan')
    out['n_set'] = cnt
    return out

def normalise(dv, dref):
    s = {}
    for k in ('d_NBD', 'g_cyt', 'd_gate'):
        fin = np.mean([dref[r][k] for r in ('7P03', '7P04', '7P05')])
        fout = dref['7P06'][k]
        s[k] = float((dv[k] - fin) / (fout - fin)) if abs(fout - fin) > 1e-06 else float('nan')
    return s

def label(s):
    if s['g_cyt'] < 0.33:
        lab = 'inward-facing'
    elif s['g_cyt'] > 0.67 and s['d_gate'] > 0.5:
        lab = 'outward-facing'
    elif s['g_cyt'] > 0.67 and (not s['d_gate'] > 0.5):
        lab = 'occluded' if s['d_gate'] == s['d_gate'] else 'closed-cytoplasm (gate undefined)'
    else:
        lab = 'intermediate'
    nb = 'NBDs closed' if s['d_NBD'] > 0.67 else 'NBDs open' if s['d_NBD'] < 0.33 else 'NBD intermediate'
    return (lab, nb)
common = set.intersection(*[set(refxyz[r]) for r in REFS])
cal = {r: descriptors(refxyz[r], common) for r in REFS}
print('Reference descriptors (A), common set of 4 refs:')
for r in REFS:
    print(' ', r, {k: round(v, 2) for k, v in cal[r].items() if k != 'n_set'})
print('  set sizes', cal['7P04']['n_set'])
STRUCTS = [k for k in info if k not in REFS]
rows = []
for name in STRUCTS:
    xyz = coords(name)
    tm = {}
    labs = {}
    for ref in ('7P04', '7P03', '7P06'):
        fp = os.path.join(L.OUT, 'align', f'{name}__to__{ref}.json')
        if not os.path.exists(fp):
            continue
        a = json.load(open(fp))
        tm[ref] = a['TM_norm_ref']
        m2n = {p[1]: p[0] for p in a['pairs']}
        sx = {rn: xyz[mn] for rn, mn in m2n.items() if mn in xyz}
        allowed = set(sx) & common
        dv = descriptors(sx, allowed)
        dr = {r: descriptors(refxyz[r], allowed) for r in REFS}
        s = normalise(dv, dr)
        lb = label(s)
        if min(dv['n_set'].values()) < 10:
            lb = ('NA (fewer than 10 aligned residues in a descriptor set; half-size or truncated model)', 'NA')
        labs[ref] = (lb, s, dv, dr['7P04']['n_set'])
    if not tm:
        continue
    best = max(tm, key=tm.get)
    (lab, nb), s, dv, nset = labs[best]
    consist = {ref: labs[ref][0][0] for ref in labs}
    rows.append({'name': name, 'evidence': info[name]['evidence'], 'TM_7P04': tm.get('7P04'), 'TM_7P03': tm.get('7P03'), 'TM_7P06': tm.get('7P06'), 'best_ref': best, 'd_NBD': dv['d_NBD'], 'g_cyt': dv['g_cyt'], 'g_ext': dv['g_ext'], 'd_gate': dv['d_gate'], 's_NBD': s['d_NBD'], 's_cyt': s['g_cyt'], 's_gate': s['d_gate'], 'state': lab, 'NBD': nb, 'label_by_ref': consist, 'n_residues_in_sets': nset})
    print(name, 'TM(7P04/7P03/7P06) =', {k: round(v, 3) for k, v in tm.items()}, 'best=', best, 's=', {k: round(v, 2) for k, v in s.items()}, '->', lab, '|', nb, '| by-ref:', consist, flush=True)
json.dump({'calibration': {r: {k: v for k, v in cal[r].items()} for r in REFS}, 'rows': rows}, open(os.path.join(L.OUT, 'state_classification.json'), 'w'), indent=1, default=str)
with open(os.path.join(L.OUT, 'state_classification.tsv'), 'w') as f:
    cols = ['name', 'evidence', 'TM_7P03', 'TM_7P04', 'TM_7P06', 'best_ref', 'd_NBD', 'g_cyt', 'g_ext', 'd_gate', 's_NBD', 's_cyt', 's_gate', 'state', 'NBD', 'label_by_ref']
    f.write('\t'.join(cols) + '\n')
    for r in rows:
        f.write('\t'.join(('%.3f' % r[c] if isinstance(r[c], float) else str(r[c]) for c in cols)) + '\n')

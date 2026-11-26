import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, gemmi
import p3lib as L
import pocklib as P
X_A, Y_FRAC = (5.0, 0.6)
TM_Z = 20.2
OUT = L.OUT
os.makedirs(OUT + '/pockets', exist_ok=True)
P.write_local_launcher(OUT + '/prank_local.sh')
st = gemmi.read_structure(L.IN + '/pdr5_copy/7P05.cif')
st.setup_entities()
st.remove_hydrogens()
st.remove_alternative_conformations()
rhq = [r for ch in st[0] for r in ch if r.name == 'RHQ'][0]
r6g = np.array([a.pos.tolist() for a in rhq])
geo = json.load(open(OUT + '/r6g_geometry_7P05.json'))
contacts = sorted((int(k) for k in geo['contacts']))
assert len(contacts) == 15
ca = {n: P.ca_by_num(f'{OUT}/pdb/{n}.pdb') for n in ('7P03', '7P04', '7P05')}
opm = P.opm_ca()
common = sorted(set(ca['7P04']) & set(opm))
R_to_opm, t_to_opm = P.kabsch([ca['7P04'][n] for n in common], [opm[n] for n in common])
rms = np.sqrt(np.mean(np.sum((np.array([ca['7P04'][n] for n in common]) @ R_to_opm.T + t_to_opm - np.array([opm[n] for n in common])) ** 2, axis=1)))
print('7P04 -> OPM frame RMS (A):', round(float(rms), 4), 'n', len(common))

def frame_transforms(name):
    if name == '7P04':
        R1, t1 = (np.eye(3), np.zeros(3))
    else:
        cm = sorted(set(ca[name]) & set(ca['7P04']))
        R1, t1 = P.kabsch([ca[name][n] for n in cm], [ca['7P04'][n] for n in cm])
    return (R1, t1)

def zopm(name, xyz):
    R1, t1 = frame_transforms(name)
    return float((R_to_opm @ (R1 @ np.asarray(xyz) + t1) + t_to_opm)[2])

def r6g_in(name):
    if name == '7P05':
        return r6g
    cm = sorted(set(ca['7P05']) & set(ca[name]))
    R, t = P.kabsch([ca['7P05'][n] for n in cm], [ca[name][n] for n in cm])
    return r6g @ R.T + t
METHODS = [('M1_P2Rank_default', None), ('M2_P2Rank_alphafold', 'alphafold'), ('M3_fpocket_default', 'fpocket')]
rows, allp = ([], {})
for name in ('7P05', '7P04', '7P03'):
    pdb = f'{OUT}/pdb/{name}.pdb'
    lig = r6g_in(name)
    cen = lig.mean(0)
    for mname, cfg in METHODS:
        if cfg == 'fpocket':
            od = P.run_fpocket(pdb, f'{OUT}/pockets/fpocket_{name}', log=f'{OUT}/pockets/fpocket_{name}.log')
            pk = P.parse_fpocket(od)
        else:
            csvp = P.run_p2rank(pdb, f'{OUT}/pockets/{mname}_{name}', config=cfg, log=f'{OUT}/pockets/{mname}_{name}.log')
            pk = P.parse_p2rank(csvp)
        allp[name, mname] = pk
        for p in pk[:60]:
            dcc = float(np.linalg.norm(p['center'] - cen))
            dca = float(np.min(np.linalg.norm(lig - p['center'], axis=1)))
            lin = set(p['residues'])
            cov = len(lin & set(contacts))
            z = zopm(name, p['center'])
            rec = {'structure': name, 'method': mname, 'rank': p['rank'], 'score': p['score'], 'n_lining': len(lin), 'DCC_A': dcc, 'DCA_A': dca, 'z_OPM_A': z, 'in_TM_zone': abs(z) <= TM_Z, 'R6G_contacts_in_pocket': cov, 'frac': cov / 15.0, 'precision': cov / max(len(lin), 1), 'L': dcc <= X_A, 'C': cov / 15.0 >= Y_FRAC, 'volume_A3': p['extra'].get('Volume', float('nan'))}
            rec['recovered'] = rec['L'] and rec['C']
            rows.append(rec)
        print(name, mname, 'n_pockets', len(pk), flush=True)
cols = ['structure', 'method', 'rank', 'score', 'n_lining', 'DCC_A', 'DCA_A', 'z_OPM_A', 'in_TM_zone', 'R6G_contacts_in_pocket', 'frac', 'precision', 'L', 'C', 'recovered', 'volume_A3']
with open(OUT + '/benchmark_pockets_top15.tsv', 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows:
        f.write('\t'.join(('%.3f' % r[c] if isinstance(r[c], float) else str(r[c]) for c in cols)) + '\n')
summ = []
for name in ('7P05', '7P04', '7P03'):
    for mname, _ in METHODS:
        rr = [r for r in rows if r['structure'] == name and r['method'] == mname]
        rank1 = rr[0]
        tm = [r for r in rr if r['in_TM_zone']]
        gate = tm[0] if tm else None
        rec_ranks = [r['rank'] for r in rr if r['recovered']]
        best_frac = max(rr, key=lambda r: r['frac'])
        summ.append({'structure': name, 'method': mname, 'rank1_overall_recovered': rank1['recovered'], 'rank1_DCC': round(rank1['DCC_A'], 1), 'rank1_frac': round(rank1['frac'], 2), 'rank1_in_TM': rank1['in_TM_zone'], 'TMfiltered_rank1_rank': gate['rank'] if gate else None, 'TMfiltered_rank1_DCC': round(gate['DCC_A'], 1) if gate else None, 'TMfiltered_rank1_DCA': round(gate['DCA_A'], 1) if gate else None, 'TMfiltered_rank1_frac': round(gate['frac'], 2) if gate else None, 'TMfiltered_rank1_n_lining': gate['n_lining'] if gate else None, 'TMfiltered_rank1_precision': round(gate['precision'], 2) if gate else None, 'GATE_pass': bool(gate and gate['recovered']), 'first_recovering_rank_in_top15': min(rec_ranks) if rec_ranks else None, 'best_frac_pocket_rank': best_frac['rank'], 'best_frac': round(best_frac['frac'], 2), 'best_frac_DCC': round(best_frac['DCC_A'], 1), 'sens_50pct_gate': bool(gate and gate['L'] and (gate['frac'] >= 0.5)), 'sens_75pct_gate': bool(gate and gate['L'] and (gate['frac'] >= 0.75))})
json.dump(summ, open(OUT + '/benchmark_summary.json', 'w'), indent=1)
for s in summ:
    print(s)

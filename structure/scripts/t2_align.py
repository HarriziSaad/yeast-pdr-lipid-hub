import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import p3lib as L
AF = {'Pdr5': 'P33302', 'Snq2': 'P32568', 'Pdr10': 'P51533', 'Pdr11': 'P40550', 'Pdr12': 'Q02785', 'Pdr15': 'Q04182', 'Pdr18': 'P53756', 'Aus1': 'Q08409', 'YOL075C': 'Q08234', 'Adp1': 'P25371', 'CaCdr1_AF': 'Q5ANA3'}
CA_ENTRIES = ['9IUK', '9IUL', '9IUM', '22FR', '9VT1', '9VTM', '9VTO', '9VSU', '9VSV', '9VSW', '9VSZ', '9VT2', '9VTN']
CG_ENTRIES = ['9S2C', '9S2H', '9S4S', '9S4T', '9S4U']
P5 = ['7P03', '7P04', '7P05', '7P06']
pdbdir = os.path.join(L.OUT, 'pdb')
aldir = os.path.join(L.OUT, 'align')
os.makedirs(pdbdir, exist_ok=True)
os.makedirs(aldir, exist_ok=True)
srcs = {}
for k, acc in AF.items():
    srcs[k] = (os.path.join(L.IN, 'alphafold', f'{acc}_model.cif'), 'C-AlphaFold')
for p in P5:
    srcs[p] = (os.path.join(L.IN, 'pdr5_copy', f'{p}.cif'), 'E-cryoEM')
for p in CA_ENTRIES + CG_ENTRIES:
    srcs[p] = (os.path.join(L.IN, 'rcsb_cdr1', f'{p}.cif'), 'E-cryoEM')
info = {}
for k, (path, ev) in srcs.items():
    st, orig = L.read_clean(path)
    rl = L.residues(st)
    icodes = [r['num'] for r in rl if r['icode']]
    noca = sum((1 for r in rl if r['ca'] is None))
    L.write_pdb(st, os.path.join(pdbdir, f'{k}.pdb'))
    nums = [r['num'] for r in rl if r['ca'] is not None]
    info[k] = {'evidence': ev, 'orig_chain': orig, 'n_res': len(rl), 'n_noCA': noca, 'n_icode': len(icodes), 'first': nums[0], 'last': nums[-1], 'seq': ''.join((r['aa'] for r in rl if r['ca'] is not None)), 'nums': nums}
json.dump(info, open(os.path.join(L.OUT, 'structure_info.json'), 'w'))
print({k: (v['n_res'], v['n_noCA'], v['n_icode']) for k, v in info.items()})

def align(mob, ref):
    out = os.path.join(aldir, f'{mob}__to__{ref}.json')
    if os.path.exists(out):
        return json.load(open(out))
    mfile = os.path.join(aldir, f'{mob}__to__{ref}.matrix')
    t0 = time.time()
    r = L.usalign(os.path.join(pdbdir, f'{mob}.pdb'), os.path.join(pdbdir, f'{ref}.pdb'), matrix_file=mfile)
    mob_ca = [(x['num'], x['ca']) for x in L.residues(L.read_clean(srcs[mob][0])[0]) if x['ca'] is not None]
    ref_ca = [(x['num'], x['ca']) for x in L.residues(L.read_clean(srcs[ref][0])[0]) if x['ca'] is not None]
    pairs = []
    for i, j in r['pairs']:
        d = float(np.linalg.norm(L.apply_rt(r['R'], r['t'], mob_ca[i][1]) - ref_ca[j][1]))
        pairs.append([mob_ca[i][0], ref_ca[j][0], round(d, 3)])
    d = np.array([p[2] for p in pairs])
    rec = {'mobile': mob, 'ref': ref, 'TM_norm_mobile': r.get('TM_norm_S1'), 'TM_norm_ref': r.get('TM_norm_S2'), 'RMSD': r['RMSD'], 'Lali': r['Lali'], 'len_mobile': r['len_S1'], 'len_ref': r['len_S2'], 'R': r['R'].tolist(), 't': r['t'].tolist(), 'n_pairs_lt3': int((d <= 3.0).sum()), 'n_pairs_lt5': int((d < 5.0).sum()), 'pairs': pairs, 'usalign_seconds': round(time.time() - t0, 1)}
    json.dump(rec, open(out, 'w'))
    return rec
rows = []
jobs = []
mobiles = list(AF) + CA_ENTRIES + CG_ENTRIES
for m in mobiles:
    for ref in ('7P04', '7P03', '7P06'):
        jobs.append((m, ref))
for m in ('7P03', '7P05', '7P06'):
    jobs.append((m, '7P04'))
jobs += [('7P03', '7P06'), ('7P04', '7P06'), ('7P05', '7P06'), ('7P04', '7P03'), ('7P05', '7P03')]
if len(sys.argv) > 1 and sys.argv[1] == 'reverse':
    jobs = jobs[::-1]
for m, ref in jobs:
    rec = align(m, ref)
    rows.append((m, ref, rec['TM_norm_ref'], rec['TM_norm_mobile'], rec['RMSD'], rec['Lali'], rec['len_mobile'], rec['len_ref'], rec['n_pairs_lt3'], rec['n_pairs_lt5']))
    print(m, '->', ref, 'TM(ref)=%.3f TM(mob)=%.3f RMSD=%.2f Lali=%d  t=%ss' % (rec['TM_norm_ref'], rec['TM_norm_mobile'], rec['RMSD'], rec['Lali'], rec['usalign_seconds']), flush=True)
with open(os.path.join(L.OUT, 'usalign_summary%s.tsv' % ('_rev' if len(sys.argv) > 1 else '')), 'w') as f:
    f.write('mobile\tref\tTM_norm_ref\tTM_norm_mobile\tRMSD_A\tLali\tlen_mobile\tlen_ref\tpairs_d_le3\tpairs_d_lt5\n')
    for r in rows:
        f.write('\t'.join((str(x) for x in r)) + '\n')

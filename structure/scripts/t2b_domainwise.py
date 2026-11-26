import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, gemmi
import p3lib as L
DOM = {'D_NBD1': (128, 490), 'D_TMD1': (491, 795), 'D_NBD2': (796, 1166), 'D_TMD2': (1167, 1511)}
MODELS = ['Pdr5', 'Snq2', 'Pdr10', 'Pdr11', 'Pdr12', 'Pdr15', 'Pdr18', 'Aus1', 'YOL075C', 'Adp1', 'CaCdr1_AF', '22FR', '9VT1', '9VTO', '9IUK', '9S2C', '9S2H', '9S4S', '9S4T', '9S4U']
odir = os.path.join(L.OUT, 'align_dom')
os.makedirs(odir, exist_ok=True)
os.makedirs(os.path.join(L.OUT, 'pdb_dom'), exist_ok=True)

def slice_pdb(src_pdb, nums_keep, dst):
    st, _ = L.read_clean(src_pdb)
    keep = set(nums_keep)
    new = gemmi.Structure()
    m = gemmi.Model('1')
    c = gemmi.Chain('A')
    for r in st[0]['A']:
        if r.seqid.num in keep:
            c.add_residue(r)
    m.add_chain(c)
    new.add_model(m)
    new.setup_entities()
    new.write_pdb(dst)
    return len(c)
ref_pdb = os.path.join(L.OUT, 'pdb', '7P04.pdb')
ref_st, _ = L.read_clean(ref_pdb)
ref_ca = {r['num']: r['ca'] for r in L.residues(ref_st) if r['ca'] is not None}
summary = []
for mname in MODELS:
    fp = os.path.join(L.OUT, 'align', f'{mname}__to__7P04.json')
    if not os.path.exists(fp):
        print('skip (no global alignment yet)', mname)
        continue
    g = json.load(open(fp))
    st, _ = L.read_clean(os.path.join(L.OUT, 'pdb', f'{mname}.pdb'))
    mob_ca = {r['num']: r['ca'] for r in L.residues(st) if r['ca'] is not None}
    mob_nums = sorted(mob_ca)
    for dname, (a, b) in DOM.items():
        out = os.path.join(odir, f'{mname}__{dname}.json')
        if os.path.exists(out):
            continue
        pr = [(p[0], p[1]) for p in g['pairs'] if a <= p[1] <= b]
        if len(pr) < 30:
            json.dump({'model': mname, 'domain': dname, 'pairs': [], 'note': 'fewer than 30 globally paired residues in the domain'}, open(out, 'w'))
            continue
        offs = np.array([m - r for m, r in pr])
        med = np.median(offs)
        pr = [(m, r) for m, r in pr if abs(m - r - med) <= 60]
        lo, hi = (min((m for m, r in pr)), max((m for m, r in pr)))
        mslice = [n for n in mob_nums if lo <= n <= hi]
        rslice = sorted((n for n in ref_ca if a <= n <= b))
        mp = os.path.join(L.OUT, 'pdb_dom', f'{mname}__{dname}.pdb')
        rp = os.path.join(L.OUT, 'pdb_dom', f'7P04__{dname}.pdb')
        slice_pdb(os.path.join(L.OUT, 'pdb', f'{mname}.pdb'), mslice, mp)
        if not os.path.exists(rp):
            slice_pdb(ref_pdb, rslice, rp)
        mat = os.path.join(odir, f'{mname}__{dname}.matrix')
        r = L.usalign(mp, rp, matrix_file=mat)
        pairs = []
        for i, j in r['pairs']:
            mnum, rnum = (mslice[i], rslice[j])
            d = float(np.linalg.norm(L.apply_rt(r['R'], r['t'], mob_ca[mnum]) - ref_ca[rnum]))
            pairs.append([mnum, rnum, round(d, 3)])
        json.dump({'model': mname, 'domain': dname, 'TM_norm_ref': r.get('TM_norm_S2'), 'RMSD': r['RMSD'], 'Lali': r['Lali'], 'len_model_slice': len(mslice), 'len_ref_slice': len(rslice), 'R': r['R'].tolist(), 't': r['t'].tolist(), 'pairs': pairs}, open(out, 'w'))
        summary.append((mname, dname, len(mslice), len(rslice), r['Lali'], r['RMSD'], r.get('TM_norm_S2')))
        print(mname, dname, 'slice', len(mslice), '/', len(rslice), 'Lali', r['Lali'], 'RMSD', r['RMSD'], 'TM(ref)', r.get('TM_norm_S2'), flush=True)
with open(os.path.join(L.OUT, 'domainwise_summary.tsv'), 'w') as f:
    f.write('model\tdomain\tlen_model_slice\tlen_ref_slice\tLali\tRMSD_A\tTM_norm_ref\n')
    for s in summary:
        f.write('\t'.join((str(x) for x in s)) + '\n')

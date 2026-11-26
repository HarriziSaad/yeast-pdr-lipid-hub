import os, re, subprocess, math, json
import numpy as np
import gemmi
W = '/home/pdr/work/phase3'
IN = os.path.join(W, 'in')
OUT = os.path.join(W, 'out')
THREE2ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

def read_clean(path):
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    model = st[0]
    best, bestn = (None, -1)
    for ch in model:
        n = sum((1 for r in ch if r.entity_type == gemmi.EntityType.Polymer))
        if n > bestn:
            best, bestn = (ch, n)
    new = gemmi.Structure()
    new.cell = st.cell
    new.spacegroup_hm = st.spacegroup_hm
    nm = gemmi.Model('1')
    nc = gemmi.Chain('A')
    for r in best:
        if r.entity_type != gemmi.EntityType.Polymer:
            continue
        if r.name not in THREE2ONE:
            pass
        nc.add_residue(r)
    nm.add_chain(nc)
    new.add_model(nm)
    new.setup_entities()
    return (new, best.name)

def residues(st):
    out = []
    for r in st[0]['A']:
        ca = r.find_atom('CA', '*')
        out.append({'num': r.seqid.num, 'icode': r.seqid.icode.strip(), 'name': r.name, 'aa': THREE2ONE.get(r.name, 'X'), 'ca': np.array(ca.pos.tolist()) if ca else None, 'b': ca.b_iso if ca else float('nan')})
    return out

def write_pdb(st, path):
    st.write_pdb(path)

def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, check=True, **kw)

def usalign(mobile_pdb, ref_pdb, matrix_file=None, extra=()):
    cmd = ['USalign', mobile_pdb, ref_pdb, '-mol', 'prot', '-outfmt', '0', *extra]
    if matrix_file:
        cmd += ['-m', matrix_file]
    p = run(cmd)
    txt = p.stdout
    res = {'raw': txt}
    m = re.search('Aligned length=\\s*(\\d+), RMSD=\\s*([\\d.]+), Seq_ID=n_identical/n_aligned=\\s*([\\d.]+)', txt)
    res['Lali'], res['RMSD'], res['seqid_aligned'] = (int(m.group(1)), float(m.group(2)), float(m.group(3)))
    tms = re.findall('TM-score=\\s*([\\d.]+) \\(normalized by length of Structure_(\\d)', txt)
    for v, s in tms:
        res[f'TM_norm_S{s}'] = float(v)
    ml = re.findall('Length of Structure_(\\d): (\\d+) residues', txt)
    for s, n in ml:
        res[f'len_S{s}'] = int(n)
    lines = txt.splitlines()
    idx = [i for i, l in enumerate(lines) if l.startswith('(":" denotes')]
    a1 = lines[idx[0] + 1]
    mid = lines[idx[0] + 2]
    a2 = lines[idx[0] + 3]
    res['aln1'], res['aln_mid'], res['aln2'] = (a1, mid, a2)
    pairs, i1, i2 = ([], 0, 0)
    for c1, c2 in zip(a1, a2):
        if c1 != '-' and c2 != '-':
            pairs.append((i1, i2))
        if c1 != '-':
            i1 += 1
        if c2 != '-':
            i2 += 1
    res['pairs'] = pairs
    res['n1'], res['n2'] = (i1, i2)
    if matrix_file:
        rows = []
        for l in open(matrix_file):
            f = l.split()
            if len(f) == 5 and f[0] in ('0', '1', '2'):
                rows.append([float(x) for x in f[1:]])
        M = np.array(rows)
        res['t'] = M[:, 0]
        res['R'] = M[:, 1:4]
    return res

def apply_rt(R, t, xyz):
    return (R @ np.asarray(xyz).T).T + t

def ca_list(st):
    return [r for r in residues(st) if r['ca'] is not None]

def load_tsv(path):
    import csv
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))

import json
import numpy as np
from a12_lib import kabsch
MAXASA = {'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0, 'GLU': 223.0, 'GLN': 225.0, 'GLY': 104.0, 'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0, 'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0, 'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0}

def read_chain(path, chain='A'):
    import gemmi
    st = gemmi.read_structure(path)
    ca, names, heavy = ({}, {}, [])
    for res in st[0][chain]:
        a = res.find_atom('CA', '*')
        if a is not None:
            ca[res.seqid.num] = (a.pos.x, a.pos.y, a.pos.z)
            names[res.seqid.num] = res.name
        for at in res:
            if at.element.name != 'H':
                heavy.append((at.pos.x, at.pos.y, at.pos.z))
    return (ca, names, np.array(heavy))

def read_ligand_molecules(path, resname='ERG'):
    import gemmi
    st = gemmi.read_structure(path)
    out = []
    for ch in st[0]:
        for res in ch:
            if res.name == resname:
                xs = [(a.pos.x, a.pos.y, a.pos.z) for a in res if a.element.name != 'H']
                nm = [a.name for a in res if a.element.name != 'H']
                out.append(dict(chain=ch.name, seqid=res.seqid.num, atoms=nm, xyz=np.array(xs)))
    return out

def read_usalign_json(path):
    d = json.load(open(path, encoding='utf-8'))
    return (np.array(d['R'], float), np.array(d['t'], float), d['pairs'], d)

def apply_rt(R, t, xyz):
    return np.asarray(xyz, float) @ np.asarray(R, float).T + np.asarray(t, float)

def check_transform(R, t, model_ca, ref_ca, pairs, max_rmsd=None):
    pr = [(m, r) for m, r, _ in pairs if r in ref_ca and m in model_ca]
    M = np.array([model_ca[m] for m, r in pr])
    Rf = np.array([ref_ca[r] for m, r in pr])
    fwd = float(np.sqrt(np.mean(np.sum((apply_rt(R, t, M) - Rf) ** 2, axis=1))))
    inv = float(np.sqrt(np.mean(np.sum((apply_rt(R, t, Rf) - M) ** 2, axis=1))))
    direction = 'forward' if fwd <= inv else 'reverse'
    rm = min(fwd, inv)
    if max_rmsd is not None and rm > max_rmsd:
        raise ValueError(f'transform does not reproduce the stored pairing (rmsd {rm:.2f} A)')
    return (direction, rm)

def model_to_ref(R, t, direction):
    if direction == 'forward':
        return (np.asarray(R, float), np.asarray(t, float))
    Ri = np.asarray(R, float).T
    return (Ri, -Ri @ np.asarray(t, float))

def to_frame_from_ca(ca_src, ca_dst):
    common = sorted(set(ca_src) & set(ca_dst))
    P = np.array([ca_src[r] for r in common])
    Q = np.array([ca_dst[r] for r in common])
    return kabsch(P, Q)

def assign_molecule_domain(centroid, model_ca, pairs_by_domain):
    best, bd = (None, 1000000000.0)
    for name, pairs in pairs_by_domain.items():
        xs = np.array([model_ca[m] for m, _, _ in pairs if m in model_ca])
        if len(xs) == 0:
            continue
        d = np.linalg.norm(xs - np.asarray(centroid), axis=1).min()
        if d < bd:
            best, bd = (name, d)
    return (best, bd)

def sterol_nodes_in_ref_frame(mols, model_ca, domain_transforms, pairs_by_domain, to_opm=None):
    out, mid, dom = ([], [], [])
    for k, m in enumerate(mols):
        name, _ = assign_molecule_domain(m['xyz'].mean(0), model_ca, pairs_by_domain)
        R, t = domain_transforms[name]
        x = apply_rt(R, t, m['xyz'])
        if to_opm is not None:
            x = apply_rt(to_opm[0], to_opm[1], x)
        out.append(x)
        mid.append(np.full(len(x), k))
        dom.append(name)
    return (np.vstack(out), np.concatenate(mid), dom)

def plddt_trim(path, chain='A', min_plddt=50.0):
    import gemmi
    st = gemmi.read_structure(path)
    keep = []
    for res in st[0][chain]:
        a = res.find_atom('CA', '*')
        if a is not None and a.b_iso >= min_plddt:
            keep.append(res.seqid.num)
    return keep

def write_pdb_in_frame(src_path, dst_path, R, t, keep=None, chain='A'):
    import gemmi
    st = gemmi.read_structure(src_path)
    st.setup_entities()
    st.remove_ligands_and_waters()
    st.remove_hydrogens()
    for ch in list(st[0]):
        if ch.name != chain:
            st[0].remove_chain(ch.name)
    ch = st[0][chain]
    if keep is not None:
        ks = set(keep)
        for i in reversed(range(len(ch))):
            if ch[i].seqid.num not in ks:
                del ch[i]
    R = np.asarray(R, float)
    t = np.asarray(t, float)
    for res in ch:
        for a in res:
            p = R @ np.array([a.pos.x, a.pos.y, a.pos.z]) + t
            a.pos = gemmi.Position(*p)
    st.write_pdb(dst_path)

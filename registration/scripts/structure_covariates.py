import sys, math
import numpy as np
import gemmi
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
cif_path, opm_path, out_path = sys.argv[1:4]
ann_path = sys.argv[4] if len(sys.argv) > 4 else None
HALF_THICKNESS = 15.2
MARGIN = 5.0
MAXASA = {'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0, 'GLU': 223.0, 'GLN': 225.0, 'GLY': 104.0, 'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0, 'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0, 'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0}
AA3TO1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
st = gemmi.read_structure(cif_path)
helix_res = set()
for h in st.helices:
    if h.start.chain_name != 'A':
        continue
    helix_res.update(range(h.start.res_id.seqid.num, h.end.res_id.seqid.num + 1))
opm = gemmi.read_structure(opm_path)
opm_z, opm_ca = ({}, {})
for r in opm[0]['A']:
    a = r.find_atom('CA', '*')
    if a:
        opm_z[r.seqid.num] = a.pos.z
        opm_ca[r.seqid.num] = (a.pos.x, a.pos.y, a.pos.z)
orig_ca = {}
for r in st[0]['A']:
    a = r.find_atom('CA', '*')
    if a:
        orig_ca[r.seqid.num] = (a.pos.x, a.pos.y, a.pos.z)
common = sorted(set(opm_ca) & set(orig_ca))
rng = np.random.default_rng(1)
sel = rng.choice(common, size=400, replace=False)
maxdiff = max((abs(math.dist(opm_ca[sel[i]], opm_ca[sel[i + 1]]) - math.dist(orig_ca[sel[i]], orig_ca[sel[i + 1]])) for i in range(0, 400, 2)))
assert maxdiff < 0.01, maxdiff
s = MMCIFParser(QUIET=True).get_structure('7P04', cif_path)
chain = s[0]['A']
for res in list(chain):
    if res.id[0] != ' ':
        chain.detach_child(res.id)
for res in chain:
    for atom in list(res):
        if atom.element == 'H' or atom.get_id().startswith('H'):
            res.detach_child(atom.get_id())
ShrakeRupley(probe_radius=1.4, n_points=960).compute(chain, level='R')
DOMAINS = [('NTE', [(1, 127)]), ('NBD1', [(128, 140), (169, 490)]), ('TMD1', [(491, 716), (774, 795)]), ('ECD1', [(717, 773)]), ('sensor', [(141, 168), (796, 861)]), ('NBD2', [(862, 1166)]), ('TMD2', [(1167, 1404), (1471, 1511)]), ('ECD2', [(1405, 1470)])]

def domain_of(n):
    for name, rngs in DOMAINS:
        if any((a <= n <= b for a, b in rngs)):
            return name
    return 'NA'

def pool_label(n, zone):
    d = domain_of(n)
    if zone:
        return 'TM'
    if d in ('NBD1', 'NBD2', 'sensor'):
        return 'NBD'
    if d in ('ECD1', 'ECD2'):
        return 'ECD'
    if d in ('TMD1', 'TMD2'):
        return 'TMDx'
    return 'NTE'
n_rows = 0
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('residue\taa\tSASA_A2\tRSA\tz_CA_OPM_A\tin_hydrophobic_slab\tin_matching_zone\tin_HELX_P\tTM_helix\tdomain\tpool_label\n')
    for res in chain:
        nm = res.get_resname()
        n = res.id[1]
        if nm not in MAXASA or n not in opm_z:
            continue
        z = opm_z[n]
        slab = abs(z) <= HALF_THICKNESS
        zone = abs(z) <= HALF_THICKNESS + MARGIN
        hx = n in helix_res
        f.write(f'{n}\t{AA3TO1[nm]}\t{res.sasa:.2f}\t{res.sasa / MAXASA[nm]:.4f}\t{z:.2f}\t{int(slab)}\t{int(zone)}\t{int(hx)}\t{int(slab and hx)}\t{domain_of(n)}\t{pool_label(n, zone)}\n')
        n_rows += 1
print(f'wrote {n_rows} residues to {out_path}; OPM frame vs 7P04 max pairwise-distance difference {maxdiff:.4f} A')
if ann_path:
    import csv
    xs, ys = ([], [])
    for r in csv.DictReader(open(ann_path, encoding='utf-8'), delimiter='\t'):
        try:
            zr = float(r['z_rel_A_7P05_CA'])
        except (ValueError, KeyError):
            continue
        n = int(r['pos'])
        if n in opm_z:
            xs.append(zr)
            ys.append(opm_z[n])
    print(f'sign check: correlation of OPM z with Phase-0 z_rel (+ = extracellular) over {len(xs)} residues = {np.corrcoef(xs, ys)[0, 1]:.4f}')
print(f"software: gemmi {gemmi.__version__}, numpy {np.__version__}, biopython {__import__('Bio').__version__}")

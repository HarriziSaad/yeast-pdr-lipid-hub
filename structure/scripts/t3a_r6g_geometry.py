import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, gemmi
import p3lib as L
st = gemmi.read_structure(L.IN + '/pdr5_copy/7P05.cif')
st.setup_entities()
st.remove_hydrogens()
st.remove_alternative_conformations()
rhq = [r for ch in st[0] for r in ch if r.name == 'RHQ']
assert len(rhq) == 1
lig = np.array([a.pos.tolist() for a in rhq[0]])
print('RHQ atoms:', len(lig), 'centroid:', lig.mean(0).round(2), 'max extent (A):', round(float(np.max(np.linalg.norm(lig[:, None] - lig[None], axis=2))), 1))
ns = gemmi.NeighborSearch(st[0], st.cell, 6).populate()
contacts = {}
for a in rhq[0]:
    for m in ns.find_atoms(a.pos, '\x00', radius=4.0):
        cra = m.to_cra(st[0])
        if cra.residue.entity_type == gemmi.EntityType.Polymer:
            contacts.setdefault(cra.residue.seqid.num, cra.residue.name)
print('R6G 4.0 A contact residues:', len(contacts))
print(sorted(contacts.items()))
cav = {int(r['residue']) for r in L.load_tsv(L.IN + '/site_classes_FROZEN.tsv') if r['class'] == 'CAV-R6G'}
print('frozen CAV-R6G:', len(cav), sorted(cav))
print('overlap:', len(cav & set(contacts)), 'only contact:', sorted(set(contacts) - cav), 'only CAV:', sorted(cav - set(contacts)))
json.dump({'centroid': lig.mean(0).tolist(), 'contacts': {str(k): v for k, v in contacts.items()}}, open(L.OUT + '/r6g_geometry_7P05.json', 'w'))

import sys, os, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, gemmi
import p3lib as L
named = {1329: 'L1329 (ergC)', 1478: 'Y1478 (ergC)', 1365: 'M1365 (ergH)', 1383: 'W1383 (ergH)'}
ENT = ['9S2C', '9S2H', '9S4S', '9S4T', '9S4U']
out = []
for e in ENT:
    st = gemmi.read_structure(f'{L.IN}/rcsb_cdr1/{e}.cif')
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    model = st[0]
    ergs = [r for c in model for r in c if r.name == 'ERG']
    ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
    contact = collections.defaultdict(set)
    for k, erg in enumerate(ergs):
        for a in erg:
            for m in ns.find_atoms(a.pos, '\x00', radius=4.0):
                cra = m.to_cra(model)
                if cra.residue.entity_type == gemmi.EntityType.Polymer:
                    contact[cra.residue.seqid.num].add(k)
    a = json.load(open(f'{L.OUT}/align/{e}__to__7P04.json'))
    dd = {}
    for dn in ('D_NBD1', 'D_TMD1', 'D_NBD2', 'D_TMD2'):
        for p in json.load(open(f'{L.OUT}/align_dom/{e}__{dn}.json'))['pairs']:
            dd[p[1]] = (p[0], p[2])
    g = {p[1]: (p[0], p[2]) for p in a['pairs']}
    for rn, lab in named.items():
        gp, dp = (g.get(rn, (None, None)), dd.get(rn, (None, None)))
        out.append((e, lab, gp[0], gp[1], dp[0], dp[1], len(contact.get(dp[0] or gp[0], ()))))
with open(f'{L.OUT}/c4b_named_residue_check.tsv', 'w') as f:
    f.write('entry\tP33302_residue\tCg_partner_global\tCA_global_A\tCg_partner_domainwise\tCA_domainwise_A\tn_ERG_within_4A_of_partner(domainwise partner)\n')
    for r in out:
        f.write('\t'.join(('' if x is None else f'{x:.2f}' if isinstance(x, float) else str(x) for x in r)) + '\n')
        print(r)

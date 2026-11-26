import sys, os, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, gemmi
import p3lib as L
ENT = ['9S2C', '9S2H', '9S4S', '9S4T', '9S4U']
info = json.load(open(os.path.join(L.OUT, 'structure_info.json')))
q6 = ''.join((l.strip() for l in open(L.IN + '/uniprot/Q6FK23.fasta') if not l.startswith('>')))
th = [(r['name'], int(r['start']), int(r['end'])) for r in L.load_tsv(L.IN + '/harris_th_helices.tsv')]

def th_of(n):
    for name, a, b in th:
        if a <= n <= b:
            return name
    return None
classes = L.load_tsv(L.IN + '/site_classes_FROZEN.tsv')
cls_of = collections.defaultdict(set)
for c in classes:
    cls_of[int(c['residue'])].add(c['class'])
p5seq = info['7P04']['seq']
p5nums = info['7P04']['nums']
p5aa = dict(zip(p5nums, p5seq))
ca04 = {}
st04, _ = L.read_clean(L.OUT + '/pdb/7P04.pdb')
for r in L.residues(st04):
    if r['ca'] is not None:
        ca04[r['num']] = r['ca']
n04 = np.array(sorted(ca04))
x04 = np.array([ca04[n] for n in n04])
per_entry_contacts, per_entry_nerg, numbering_check = ({}, {}, {})
for e in ENT:
    st = gemmi.read_structure(f'{L.IN}/rcsb_cdr1/{e}.cif')
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    model = st[0]
    ch = max(model, key=lambda c: sum((1 for r in c if r.entity_type == gemmi.EntityType.Polymer)))
    prot = {r.seqid.num: r for r in ch if r.entity_type == gemmi.EntityType.Polymer}
    mism = [n for n, r in prot.items() if L.THREE2ONE.get(r.name, 'X') != (q6[n - 1] if 0 < n <= len(q6) else '?')]
    numbering_check[e] = (len(prot), len(mism))
    ergs = [r for c in model for r in c if r.name == 'ERG']
    per_entry_nerg[e] = len(ergs)
    ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
    contact = collections.defaultdict(set)
    for k, erg in enumerate(ergs):
        for a in erg:
            for m in ns.find_atoms(a.pos, '\x00', radius=4.0):
                cra = m.to_cra(model)
                if cra.chain.name == ch.name and cra.residue.entity_type == gemmi.EntityType.Polymer:
                    contact[cra.residue.seqid.num].add(k)
    per_entry_contacts[e] = {n: len(v) for n, v in contact.items()}
    print(e, 'ERG', len(ergs), 'contact residues', len(contact), 'numbering (n, mismatches)', numbering_check[e], flush=True)
cnt = collections.Counter()
for e in ENT:
    for n in per_entry_contacts[e]:
        cnt[n] += 1
K = sorted((n for n, c in cnt.items() if c >= 3))
print('Cdr1 residues with ERG contact in >= 3 of 5 entries:', len(K))

def generic_map(e, variant='global'):
    st = gemmi.read_structure(f'{L.IN}/rcsb_cdr1/{e}.cif')
    st.setup_entities()
    st.remove_alternative_conformations()
    model = st[0]
    ch = max(model, key=lambda c: sum((1 for r in c if r.entity_type == gemmi.EntityType.Polymer)))
    cacg = {r.seqid.num: np.array(r.find_atom('CA', '*').pos.tolist()) for r in ch if r.entity_type == gemmi.EntityType.Polymer and r.find_atom('CA', '*')}
    hel = [(h.start.res_id.seqid.num, h.end.res_id.seqid.num) for h in st.helices if h.start.chain_name == ch.name]
    a = json.load(open(f'{L.OUT}/align/{e}__to__7P04.json'))
    R, t = (np.array(a['R']), np.array(a['t']))
    pair = {p[0]: (p[1], p[2]) for p in a['pairs']}
    rt = {n: (R, t) for n in pair}
    if variant == 'domain':
        pair, rt = ({}, {})
        for dn in ('D_NBD1', 'D_TMD1', 'D_NBD2', 'D_TMD2'):
            dj = json.load(open(f'{L.OUT}/align_dom/{e}__{dn}.json'))
            if not dj['pairs']:
                continue
            for p in dj['pairs']:
                pair[p[0]] = (p[1], p[2])
                rt[p[0]] = (np.array(dj['R']), np.array(dj['t']))
    hel_th = {}
    for s, t_ in hel:
        votes = collections.Counter()
        for n in range(s, t_ + 1):
            if n in pair and pair[n][1] <= 3.0:
                votes[th_of(pair[n][0])] += 1
        if votes:
            top, c = votes.most_common(1)[0]
            if top is not None and c / (t_ - s + 1) >= 0.5:
                hel_th[s, t_] = top

    def helix_of(n):
        for (s, t_), name in hel_th.items():
            if s <= n <= t_:
                return name
        return None
    out = {}
    for n in cacg:
        why, tgt, dist = (None, None, None)
        if n not in pair:
            why = 'no_aligned_residue'
        else:
            tgt, dist = pair[n]
            Rn, tn = rt[n]
            near = n04[np.linalg.norm(x04 - (Rn @ cacg[n] + tn), axis=1) <= 3.0]
            if dist > 3.0:
                why = 'CA_gt3A'
            elif len(near) >= 2:
                why = 'ambiguous_2_candidates_within_3A'
            elif th_of(tgt) is None:
                why = 'P33302_partner_not_in_Harris_TM_helix'
            elif helix_of(n) is None:
                why = 'Cdr1_residue_not_in_a_helix_segment_matching_a_TH'
            elif helix_of(n) != th_of(tgt):
                why = 'different_helix'
        out[n] = {'accepted': why is None, 'reason': why, 'P33302': tgt, 'dist': dist, 'cdr1_helix': helix_of(n), 'p5_helix': th_of(tgt) if tgt else None}
    return out
gm = {e: generic_map(e) for e in ENT}
gmd = {e: generic_map(e, 'domain') for e in ENT} if os.path.exists(f'{L.OUT}/align_dom/9S4T__D_TMD2.json') else None
kept, disc = ([], [])
for n in K:
    acc = {e: gm[e][n] for e in ENT if n in gm[e] and gm[e][n]['accepted']}
    targets = {v['P33302'] for v in acc.values()}
    cname = q6[n - 1]
    accd = {e: gmd[e][n] for e in ENT if n in gmd[e] and gmd[e][n]['accepted']} if gmd else {}
    tgtd = {v['P33302'] for v in accd.values()}
    row = {'n_accepted_domainwise': len(accd), 'domainwise_target': ','.join((str(x) for x in sorted(tgtd))), 'Cdr1_Cg_residue': n, 'Cdr1_aa': cname, 'n_entries_ERG_contact': cnt[n], 'entries_ERG_contact': ','.join((e for e in ENT if n in per_entry_contacts[e])), 'n_entries_mapping_accepted': len(acc)}
    if len(targets) == 1:
        t = next(iter(targets))
        row.update({'P33302_residue': t, 'P33302_aa': p5aa[t], 'Harris_helix': th_of(t), 'median_CA_dist_A': float(np.median([v['dist'] for v in acc.values()])), 'aa_same': 'Y' if p5aa[t] == cname else 'N', 'prereg_classes_of_P33302_residue': ';'.join(sorted(cls_of.get(t, []))) or '-'})
        kept.append(row)
    else:
        reasons = collections.Counter((gm[e][n]['reason'] if n in gm[e] else 'absent_from_entry' for e in ENT))
        row['discard_reason'] = 'conflicting_accepted_mappings:' + ','.join((str(x) for x in sorted(targets))) if len(targets) > 1 else 'no_accepted_mapping(' + ';'.join((f'{k}:{v}' for k, v in reasons.items())) + ')'
        disc.append(row)
hdr = ['Cdr1_Cg_residue', 'Cdr1_aa', 'n_entries_ERG_contact', 'entries_ERG_contact', 'n_entries_mapping_accepted', 'P33302_residue', 'P33302_aa', 'Harris_helix', 'median_CA_dist_A', 'aa_same', 'prereg_classes_of_P33302_residue', 'n_accepted_domainwise', 'domainwise_target', 'tag']
with open(L.OUT + '/c4b_proxy_LH-erg-CgCdr1proxy.tsv', 'w') as f:
    f.write('\t'.join(hdr) + '\n')
    for r in sorted(kept, key=lambda r: r['P33302_residue']):
        r['tag'] = 'exploratory proxy, PDB-to-paper mapping inferred'
        f.write('\t'.join(('%.2f' % r[h] if isinstance(r.get(h), float) else str(r.get(h, '')) for h in hdr)) + '\n')
with open(L.OUT + '/c4b_discarded.tsv', 'w') as f:
    f.write('Cdr1_Cg_residue\tCdr1_aa\tn_entries_ERG_contact\tentries_ERG_contact\tn_entries_mapping_accepted\tdiscard_reason\n')
    for r in disc:
        f.write(f"{r['Cdr1_Cg_residue']}\t{r['Cdr1_aa']}\t{r['n_entries_ERG_contact']}\t{r['entries_ERG_contact']}\t{r['n_entries_mapping_accepted']}\t{r['discard_reason']}\n")
with open(L.OUT + '/c4b_per_entry_contacts.tsv', 'w') as f:
    f.write('entry\tn_ERG\tn_contact_residues\tCdr1_numbering_mismatches_vs_Q6FK23\tn_polymer_residues\n')
    for e in ENT:
        f.write(f'{e}\t{per_entry_nerg[e]}\t{len(per_entry_contacts[e])}\t{numbering_check[e][1]}\t{numbering_check[e][0]}\n')
names = {'L1329', 'Y1478', 'M1365', 'W1383'}
kept_names = {f"{r['P33302_aa']}{r['P33302_residue']}" for r in kept}
print('kept', len(kept), 'discarded', len(disc), '| named ergosterol residues present among kept:', sorted(names & kept_names), '| accepted in >=3 entries:', sum((1 for r in kept if r['n_entries_mapping_accepted'] >= 3)))
print(collections.Counter((r['discard_reason'].split('(')[0] for r in disc)))

import sys, csv, collections
draft_path, out_path, sum_path, fasta_path = sys.argv[1:5]
P33302 = ''.join((l.strip() for l in open(fasta_path, encoding='utf-8') if not l.startswith('>')))
assert len(P33302) == 1511
TOKEN_MAP = [('E structure (computed here)', 'COMP-contact-rule(this work; 4.0A heavy atoms on 7P0x)'), ('E structure', 'EXP-structure(literature list/figure)'), ('E mutagenesis', 'EXP-mutagenesis'), ('UNVERIFIED', 'UNVERIFIED')]

def map_evidence(ev):
    toks = [t.strip() for t in ev.split(';') if t.strip()]
    out = []
    for t in toks:
        for k, v in TOKEN_MAP:
            if t == k:
                out.append(v)
                break
        else:
            out.append(t)
    return out
draft = list(csv.DictReader(open(draft_path, encoding='utf-8'), delimiter='\t'))
aa = {i + 1: ch for i, ch in enumerate(P33302)}
by_class_rows = collections.defaultdict(list)
for r in draft:
    assert aa[int(r['residue'])] == r['aa'], r['residue']
    by_class_rows[r['class']].append(r)

def loc_of(r):
    s = r['TM_helix_or_domain']
    return s.split(' [')[0].strip()
F = {}

def add(res, cls, status, source=None, evidence=(), note=None):
    d = F.setdefault((res, cls), dict(status=status, source=[], evidence=[], notes=[]))
    assert d['status'] == status, (res, cls, d['status'], status)
    if source and source not in d['source']:
        d['source'].append(source)
    for e in evidence:
        if e not in d['evidence']:
            d['evidence'].append(e)
    if note and note not in d['notes']:
        d['notes'].append(note)
PRE = 'Nentwig et al. 2026 bioRxiv v1 (doi:10.64898/2026.10.29.762749; unreviewed)'

def rows_of(cls, pred=lambda r: True):
    return [r for r in by_class_rows[cls] if pred(r)]
TH12_MIN = 1473
C7_NOTE = 'C7: TMH12 region absent from S288c Pdr18 (Dias 2025 as relayed; to be re-verified in Phase 1) -> "not alignable" in sequences lacking it; use an intact Pdr18 allele/ortholog for structure mapping'
LH_SRC = {'LH-MIPC-sub': PRE + ' Fig. 2 legend: residues within 4 A of the M(IP)2C ligand in the substrate-binding site (TMH2, TMH4, TMH5a, TMH8, TMH11a); PRIMARY list', 'LH-MIPC-exit': PRE + ' Fig. 2 legend: residues within 4 A of the M(IP)2C OL (exit site) on TMH10, TMH11a-c, TMH12', 'LH-ergC': PRE + ' Results (cholesterol clash argument) and Discussion: ergosterol-central site; only L1329 and Y1478 are named; PDB UNRELEASED unreleased', 'LH-ergH': PRE + ' Results/Discussion: erg-horizontal site; only M1365 and W1383 are named; PDB UNRELEASED unreleased'}
LH_EVID = {'LH-MIPC-sub': ['EXP-structure(literature list/figure)'], 'LH-MIPC-exit': ['EXP-structure(literature list/figure)'], 'LH-ergC': ['EXP-structure(literature list/figure)', 'COMP-docking(authors; cholesterol clash)'], 'LH-ergH': ['EXP-structure(literature list/figure)']}
LH_STATUS = {'LH-MIPC-sub': 'exploratory', 'LH-MIPC-exit': 'exploratory', 'LH-ergC': 'per-residue report only (no class test)', 'LH-ergH': 'per-residue report only (no class test)'}
lh_members = {}
for cls in ('LH-MIPC-sub', 'LH-MIPC-exit', 'LH-ergC', 'LH-ergH'):
    lh_members[cls] = sorted({int(r['residue']) for r in by_class_rows[cls]})
    for r in rows_of(cls):
        n = int(r['residue'])
        note = 'loc: ' + loc_of(r)
        add(n, cls, LH_STATUS[cls], LH_SRC[cls], LH_EVID[cls], note)
        if n >= TH12_MIN:
            add(n, cls, LH_STATUS[cls], None, (), C7_NOTE)
lh_conf = sorted(set().union(*[set(v) for v in lh_members.values()]))
assert len(lh_conf) == 40, len(lh_conf)
for n in lh_conf:
    comps = [c for c in lh_members if n in lh_members[c]]
    ev = []
    for c in comps:
        for e in LH_EVID[c]:
            if e not in ev:
                ev.append(e)
    add(n, 'LH-CONF', 'confirmatory', 'union of literature lists (C4): ' + ' + '.join(comps) + '; ' + PRE, ev, 'component(s): ' + ','.join(comps) + '; used in H1a, H5a')
    if n >= TH12_MIN:
        add(n, 'LH-CONF', 'confirmatory', None, (), C7_NOTE)
S6D_EXTRA = [530, 569, 573, 574, 614, 1357, 1361]
cav_draft = {int(r['residue']) for r in by_class_rows['CAV-R6G']}
for n in S6D_EXTRA:
    assert n in cav_draft and n in aa
    add(n, 'LH-S6d-extra', 'sensitivity set (derived)', PRE + ' Fig. S6d labelled 4 A residues for the substrate-site M(IP)2C that are NOT in the Fig. 2 legend list (extraction: nentwig_preprint_extraction.md section 6); reconciliation with Fig. 2 UNVERIFIED', ['EXP-structure(literature list/figure)'], 'also in CAV-R6G; used only in the union sensitivity analysis')
    add(n, 'LH-CONF+S6d', 'sensitivity set (derived)', 'LH-CONF U Fig. S6d extras', ['EXP-structure(literature list/figure)'], None)
for n in lh_conf:
    add(n, 'LH-CONF+S6d', 'sensitivity set (derived)', 'LH-CONF U Fig. S6d extras', [], None)
for r in rows_of('CAV-R6G'):
    n = int(r['residue'])
    ev = map_evidence(r['evidence_type'])
    src = 'Harris 2021 (doi:10.1038/s41467-021-25574-8) Supp Fig 10c LigPlot list (read from the rendered figure) and/or this work: heavy atoms within 4.0 A of rhodamine 6G (RHQ1603) in 7P05'
    add(n, 'CAV-R6G', 'confirmatory', src, ev, 'loc: ' + loc_of(r) + '; used in H1b; contact-rule part is computational')
ALH = 'Alhumaidi et al. 2022 JBC 298:102689 (doi:10.1016/j.jbc.2022.102689)'
GATE_CORE = [682, 683, 1371, 1373]
for n in GATE_CORE:
    rr = [r for r in by_class_rows['GATE'] if int(r['residue']) == n]
    ev = []
    for r in rr:
        for e in map_evidence(r['evidence_type']):
            if e not in ev:
                ev.append(e)
    add(n, 'GATE-core', 'confirmatory', ALH + ' "The gate proper comprises Gly-682, Phe-683, Gly-1371, and Met-1373" (Results); G0 memo C3', ev, 'loc: ' + loc_of(rr[0]) + '; used in H5b (confirmatory), H4 (exploratory)')
region = sorted((set(range(680, 686)) | {1372, 1368, 1369}) - set(GATE_CORE))
for n in region:
    rr = [r for r in by_class_rows['GATE'] if int(r['residue']) == n]
    if rr:
        ev = []
        for r in rr:
            for e in map_evidence(r['evidence_type']):
                if e not in ev:
                    ev.append(e)
        note = 'loc: ' + loc_of(rr[0])
        src = ALH + ' (gating region / allosteric site: I685, V1372; "molecular diode" S1368; highly conserved gating-region residue F1369); G0 memo C3 + coordinator review'
    else:
        ev = []
        note = 'residue is inside the C3 range Y680-I685 but is not individually named in the sources read; included because C3 defines the region as a range'
        src = 'G0 memo C3 range Y680-I685'
    if n in (680, 681):
        ev = ['COMP-range(G0 memo C3; not individually named in the literature)']
        src = 'G0 memo C3 range Y680-I685 (stretch described in ' + ALH + ' and Harris 2021' + ('; Y680 is named only as the start of the stretch)' if n == 680 else '; T681 is not individually named)')
    aa_n = aa.get(n)
    add(n, 'GATE-region', 'exploratory', src, ev, note + '; disjoint from GATE-core (G682/F683 are in GATE-core only); used in H4')
GATE_ALL = sorted(set(GATE_CORE) | set(region))
for n in GATE_ALL:
    add(n, 'GATE-ALL', 'sensitivity set (derived)', 'GATE-core U GATE-region', [], None)
MOTIF_PREFIX = ('Walker A', 'Walker B', 'Q-loop', 'D-loop', 'H-loop', 'Signature')
H21 = 'Harris 2021 (doi:10.1038/s41467-021-25574-8) text, Fig 4, Supp Fig 10 LigPlot; Egner 1998 Table 4'
nbs_pool_excl = {}
for cls, tag in (('NBS1', 'NBS1'), ('NBS2', 'NBS2')):
    motif_res = set()
    for r in by_class_rows[cls]:
        sub = r['subclass']
        if sub.startswith(MOTIF_PREFIX):
            motif_res.add(int(r['residue']))
    for r in by_class_rows[cls]:
        n = int(r['residue'])
        sub = r['subclass']
        ev = map_evidence(r['evidence_type'])
        if sub.startswith('A-loop candidate'):
            nbs_pool_excl[n] = 'A-loop candidate (position-based inference; no paper names A-loop residues): UNVERIFIED; excluded from classes (C3-style rule), excluded from control pools'
            continue
        if sub.startswith(MOTIF_PREFIX):
            add(n, tag + '-motif', 'exploratory', H21 + '; sequence motif: ' + sub, ev, 'loc: ' + loc_of(r) + '; motif membership is literature/sequence-defined; ligand-shell evidence is computational')
        elif n not in motif_res:
            add(n, tag + '-contact', 'exploratory', 'this work: heavy atoms within 4.0 A of the bound nucleotide/Mg2+ in 7P04/7P05/7P06 (' + sub + ')', ev, 'loc: ' + loc_of(r) + '; contact-rule-defined (computational)')
motif_all = sorted({n for n, c in F if c in ('NBS1-motif', 'NBS2-motif')})
for n in motif_all:
    comps = [c for c in ('NBS1-motif', 'NBS2-motif') if (n, c) in F]
    ev = []
    for c in comps:
        for e in F[n, c]['evidence']:
            if e not in ev:
                ev.append(e)
    add(n, 'NBS-MOTIF', 'confirmatory', 'union of ' + ' + '.join(comps), ev, 'used in H0b (positive control, lower tail), H5c')
for r in by_class_rows['CLAMP/linker']:
    n = int(r['residue'])
    add(n, 'CLAMP', 'exploratory', 'Harris 2021 (nucleotide-sensing linker "clamping switch", Fig 3 labels, Supp Fig 7/8): ' + r['subclass'], map_evidence(r['evidence_type']), 'loc: ' + loc_of(r))
for r in by_class_rows['TRANSMISSION']:
    n = int(r['residue'])
    ev = map_evidence(r['evidence_type'])
    if ev == ['EXP-mutagenesis']:
        src = 'literature mutant of the signalling interface (' + r['source_and_location'][:90] + ')'
    else:
        src = 'this work: cis NBD-TMD heavy-atom contact <= 4.0 A in 7P04/7P05 (helix-level definition of Harris 2021: ' + r['subclass'] + ')'
    add(n, 'TRANSMISSION', 'exploratory', src, ev, 'loc: ' + loc_of(r) + '; no literature residue list exists; rule-based class')
for r in by_class_rows['FK506-ECL']:
    n = int(r['residue'])
    add(n, 'FK506-ECL', 'exploratory', 'Tanabe 2019 (FK506-resistance mutants; residue placement on a pre-7P03 homology model): ' + r['subclass'], map_evidence(r['evidence_type']) + ['COMP-homology-model(residue placement)'], 'loc: ' + loc_of(r))
for r in by_class_rows['MUTANT-only']:
    n = int(r['residue'])
    add(n, 'MUTANT-only', 'bookkeeping (not tested)', r['source_and_location'][:140], map_evidence(r['evidence_type']), 'published mutant not in any structural class; excluded from control pools; loc: ' + loc_of(r))
placed = {n for n, c in F if c not in ('LH-CONF+S6d', 'LH-S6d-extra', 'GATE-ALL')}
for r in draft:
    n = int(r['residue'])
    if n in placed:
        continue
    if r['class'] == 'UNVERIFIED-relay':
        why = 'L1239 relayed via Shin 2026 ("corresponds to Cdr1 F1230"); not found in Alhumaidi 2022 or Harris 2021: UNVERIFIED; excluded from GATE (C3)'
        ev = ['UNVERIFIED']
    elif r['class'] == 'GATE':
        why = 'named in Alhumaidi 2022 as an interaction partner in the network around the gate (S1366 interacts with F1369); not itself listed as a gate or gating-region residue; no frozen class'
        ev = map_evidence(r['evidence_type'])
    elif n in nbs_pool_excl:
        why = nbs_pool_excl[n]
        ev = ['UNVERIFIED']
    else:
        raise SystemExit(f"unplaced draft residue {n} class {r['class']}")
    add(n, 'POOL-EXCLUDED', 'excluded (not tested)', why, ev, 'excluded from every control pool')
for n, why in nbs_pool_excl.items():
    if (n, 'POOL-EXCLUDED') not in F:
        add(n, 'POOL-EXCLUDED', 'excluded (not tested)', why, ['UNVERIFIED'], 'excluded from every control pool')
EXTRA = {1365: 'preprint helix-label inconsistency: Fig. S13 calls M1365 TMH11b, the Fig. 2 legend lists L1367 (after M1365) on TMH11a; deposited-structure helix is TH11a (Harris naming)', 1383: 'preprint states one residue is "99.6 % conserved among the PDR family"; whether that is W1383 or M1365 is ambiguous in the text (W1383 is only inferred); FK506-resistance mutant W1383L (Tanabe 2019)', 1329: 'preprint quotes 78 % conservation (500 PDR-subfamily sequences; alignment source not stated)', 1478: 'preprint quotes 32 % conservation (500 PDR-subfamily sequences; alignment source not stated)'}
for (n, c), d in F.items():
    if n in EXTRA and c in ('LH-CONF', 'LH-ergC', 'LH-ergH'):
        d['notes'].append(EXTRA[n])
PRIMARY = ['LH-CONF', 'CAV-R6G', 'GATE-core', 'NBS-MOTIF', 'LH-MIPC-sub', 'LH-MIPC-exit', 'GATE-region', 'NBS1-motif', 'NBS2-motif', 'NBS1-contact', 'NBS2-contact', 'CLAMP', 'TRANSMISSION', 'FK506-ECL']
members = collections.defaultdict(set)
for n, c in F:
    members[c].add(n)
cav = members['CAV-R6G']
lh = members['LH-CONF']
for n in sorted(lh - cav):
    add(n, 'LH-ONLY', 'sensitivity set (derived)', 'LH-CONF minus CAV-R6G (C5)', [], None)
for n in sorted(cav - lh):
    add(n, 'CAV-ONLY', 'sensitivity set (derived)', 'CAV-R6G minus LH-CONF (C5)', [], None)
members = collections.defaultdict(set)
for n, c in F:
    members[c].add(n)
MULTI_SCOPE = PRIMARY + ['LH-ergC', 'LH-ergH']
for (n, c), d in F.items():
    if c in MULTI_SCOPE:
        others = [o for o in MULTI_SCOPE if o != c and n in members[o] and (not (c == 'LH-CONF' and o in ('LH-MIPC-sub', 'LH-MIPC-exit', 'LH-ergC', 'LH-ergH') or (o == 'LH-CONF' and c in ('LH-MIPC-sub', 'LH-MIPC-exit', 'LH-ergC', 'LH-ergH')) or (c == 'NBS-MOTIF' and o in ('NBS1-motif', 'NBS2-motif')) or (o == 'NBS-MOTIF' and c in ('NBS1-motif', 'NBS2-motif'))))]
        if others:
            d['notes'].append('multi-class (C5): also in ' + ','.join(others))
ORDER = ['LH-CONF', 'LH-MIPC-sub', 'LH-MIPC-exit', 'LH-ergC', 'LH-ergH', 'CAV-R6G', 'GATE-core', 'GATE-region', 'NBS-MOTIF', 'NBS1-motif', 'NBS2-motif', 'NBS1-contact', 'NBS2-contact', 'CLAMP', 'TRANSMISSION', 'FK506-ECL', 'LH-S6d-extra', 'LH-CONF+S6d', 'LH-ONLY', 'CAV-ONLY', 'GATE-ALL', 'MUTANT-only', 'POOL-EXCLUDED']
assert {c for _, c in F} <= set(ORDER), {c for _, c in F} - set(ORDER)
keys = sorted(F, key=lambda k: (ORDER.index(k[1]), k[0]))
clean = lambda s: s.replace('\t', ' ').replace('\n', ' ')
with open(out_path, 'w', encoding='utf-8', newline='') as f:
    f.write('residue\taa\tclass\tconfirmatory/exploratory\tdefinition source\tevidence type\tnotes\n')
    for n, c in keys:
        d = F[n, c]
        f.write('\t'.join([str(n), aa[n], c, d['status'], clean(' || '.join(d['source'])), clean('; '.join(d['evidence'])) or 'derived (inherits)', clean('; '.join(d['notes']))]) + '\n')
lines = ['site_classes_FROZEN summary (generated by freeze_classes.py)', '']
lines.append('class\tstatus\tn_residues\tn_TMH12(1473-1500)\tn_EXP_only\tn_with_COMP_contact_rule\tn_COMP_range_model_or_docking')
for c in ORDER:
    ms = sorted(members[c]) if c in members else []
    if not ms:
        continue
    st = F[ms[0], c]['status']
    ev = {n: ' '.join(F[n, c]['evidence']) for n in ms}
    n_cc = sum((1 for n in ms if 'COMP-contact-rule' in ev[n]))
    n_oc = sum((1 for n in ms if any((t in ev[n] for t in ('COMP-range', 'COMP-homology', 'COMP-docking')))))
    n_exp = sum((1 for n in ms if 'COMP' not in ev[n] and ev[n] != ''))
    lines.append(f'{c}\t{st}\t{len(ms)}\t{sum((1 for n in ms if n >= TH12_MIN and n <= 1500))}\t{n_exp}\t{n_cc}\t{n_oc}')
lines.append('')
lines.append('residues in more than one primary/tested class (C5):')
mc = collections.defaultdict(list)
for c in MULTI_SCOPE:
    if c in ('LH-CONF', 'NBS-MOTIF'):
        continue
    for n in members[c]:
        mc[n].append(c)
for n in sorted(mc):
    if len(mc[n]) > 1:
        lines.append(f'{n}{aa[n]}\t' + ','.join(mc[n]))
lines.append('')
lines.append('pairwise overlaps between confirmatory classes: ' + '; '.join((f'{a}&{b}={len(members[a] & members[b])}' for a in ('LH-CONF', 'CAV-R6G', 'GATE-core', 'NBS-MOTIF') for b in ('LH-CONF', 'CAV-R6G', 'GATE-core', 'NBS-MOTIF') if a < b)))
lines.append(f'distinct residues in the frozen table: {len({n for n, _ in F})}')
lines.append(f'pool-excluded (named, no frozen class): ' + ','.join((f'{n}{aa[n]}' for n in sorted(members['POOL-EXCLUDED']))))
open(sum_path, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
print('\n'.join(lines))

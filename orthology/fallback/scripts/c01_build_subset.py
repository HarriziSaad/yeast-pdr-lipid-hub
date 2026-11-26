import csv, collections, hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
P1C = os.path.dirname(HERE)
PROJ = os.path.dirname(os.path.dirname(P1C))
P1A = os.path.join(PROJ, 'phase1_orthology', 'phase1a')
P1B = os.path.join(PROJ, 'phase1_orthology', 'phase1b')
OUT = os.path.join(P1C, 'outputs')
os.makedirs(OUT, exist_ok=True)
MANIFEST = os.path.join(P1B, 'outputs', 'set', 'analysis_set_manifest.tsv')
MANIFEST_SHA = '735b5231a4b9cdde5372dd16f679d96f1664c07269897e8c6e7d166f920c4f85'
YGOB_MAP = os.path.join(P1B, 'outputs', 'set', 'ygob', 'ygob_map.tsv')
YGOB_MAP_SHA = 'f69b7c40639b031486d6f386e669b268f27bcd73fdd6d6467003f5cc8a2c7157'
YGOB_MATCH = os.path.join(P1B, 'outputs', 'set', 'sampling', 'ygob_match.tsv')
DIAS_MATCH = os.path.join(P1B, 'outputs', 'set', 'sampling', 'dias_match.tsv')
DIAS_A1 = os.path.join(P1A, 'inputs', 'dias2025_tableA1_strains.tsv')
SPECIES_TABLE = os.path.join(P1A, 'outputs', 'qc_extra', 'species_table.tsv')
G1_LABELS = {'mafft': os.path.join(P1B, 'outputs', 'g1', 'fast_mafft', 'g1_tree_clade_labels.tsv'), 'muscle': os.path.join(P1B, 'outputs', 'g1', 'fast_muscle', 'g1_tree_clade_labels.tsv')}
GROENEWALD_SACCHAROMYCETACEAE = {'Cyniclomyces', 'Eremothecium', 'Ashbya', 'Grigorovia', 'Hagleromyces', 'Kazachstania', 'Kluyveromyces', 'Lachancea', 'Nakaseomyces', 'Naumovozyma', 'Saccharomyces', 'Savitreea', 'Tetrapisispora', 'Torulaspora', 'Vanderwaltozyma', 'Yueomyces', 'Zygosaccharomyces', 'Zygotorulaspora'}
SUBSET_ORDER = 'Saccharomycetales'
S288C = {'PDR5': 'Saccharomyces_cerevisiae|g002312.m1', 'PDR15': 'Saccharomyces_cerevisiae|g004097.m1', 'PDR10': 'Saccharomyces_cerevisiae|g002467.m1', 'SNQ2': 'Saccharomyces_cerevisiae|g003736.m1', 'PDR18': 'Saccharomyces_cerevisiae|g002018.m1', 'PDR12': 'Saccharomyces_cerevisiae|g002739.m1', 'PDR11': 'Saccharomyces_cerevisiae|g005516.m1', 'AUS1': 'Saccharomyces_cerevisiae|g002189.m1'}
YOL075C = 'Saccharomyces_cerevisiae|g002110.m1'

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def read_tsv(path):
    with open(path, encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path, rows, fields=None):
    fields = fields or list(rows[0].keys())
    with open(path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
checks = []

def check(ok, msg):
    checks.append(('OK  ' if ok else 'FAIL') + ' ' + msg)
    print(checks[-1])
    if not ok:
        with open(os.path.join(OUT, 'subset_checks.txt'), 'w', encoding='utf-8') as f:
            f.write('\n'.join(checks) + '\n')
        sys.exit('check failed: ' + msg)
check(sha256(MANIFEST) == MANIFEST_SHA, 'analysis_set_manifest.tsv SHA-256 equals the A5 value')
check(sha256(YGOB_MAP) == YGOB_MAP_SHA, 'ygob_map.tsv SHA-256 equals the A5 value')
man = read_tsv(MANIFEST)
check(len(man) == 1318, 'manifest has 1,318 rows (1,253 ingroup + 35 OG658 + 30 Dikarya)')
roles = collections.Counter((r['role'] for r in man))
check(roles == {'INGROUP': 1253, 'OUT_OG658': 35, 'OUT_DIKARYA': 30}, 'role counts %s' % dict(roles))
sp_tab = read_tsv(SPECIES_TABLE)
genus_of = {r['speciesCode']: r['genus'] for r in sp_tab if r['kind'] == 'y1000_genome'}
fam2018_of = {r['speciesCode']: r['family_2018'] for r in sp_tab if r['kind'] == 'y1000_genome'}
order_of_tab = {r['speciesCode']: r['order'] for r in sp_tab if r['kind'] == 'y1000_genome'}
sub = [r for r in man if r['role'] == 'INGROUP' and r['order'] == SUBSET_ORDER or r['role'] == 'OUT_OG658']
ing = [r for r in sub if r['role'] == 'INGROUP']
og = [r for r in sub if r['role'] == 'OUT_OG658']
drop = [r for r in man if r not in sub]
check(all((order_of_tab[r['species_code']] == r['order'] for r in man if r['role'] != 'OUT_DIKARYA')), 'manifest order label equals species-table order label for all Y1000+ rows')
genera = collections.Counter((genus_of[r['species_code']] for r in ing))
off_list = sorted((g for g in genera if g not in GROENEWALD_SACCHAROMYCETACEAE))
check(not off_list, 'every genus of the retained ingroup genomes is on the Groenewald 2023 Saccharomycetaceae list (off-list: %s)' % off_list)
fam18 = collections.Counter((fam2018_of[r['species_code']] for r in {r['species_code']: r for r in ing}.values()))
check(set(fam18) <= {'', 'Saccharomycetaceae'}, 'family_2018 of retained genomes is blank or Saccharomycetaceae only: %s' % dict(fam18))
dropped_in_list = sorted({(r['species_code'], r['order']) for r in drop if r['role'] == 'INGROUP' and genus_of[r['species_code']] in GROENEWALD_SACCHAROMYCETACEAE})
check(not dropped_in_list, 'no dropped ingroup genome has a genus on the Saccharomycetaceae list: %s' % dropped_in_list)
ids = {r['seq_id'] for r in sub}
check(all((v in ids for v in S288C.values())), 'all eight S288C PDR sequences are in the subset')
check(all((next((r for r in man if r['seq_id'] == v))['role'] == 'INGROUP' for v in S288C.values())), 'the eight S288C sequences are INGROUP rows')
check(YOL075C in ids, 'S. cerevisiae YOL075C (OG0000658 anchor) is in the subset outgroup')
ygm = read_tsv(YGOB_MAP)
mapped = [r['mapped_seq'] for r in ygm if r['mapped_seq']]
check(all((m in ids for m in mapped)), 'all %d mapped YGOB sequences (ygob_map.tsv) are in the subset' % len(mapped))
yg_match = read_tsv(YGOB_MATCH)
yg_genomes = {r['y1000_genome'] for r in yg_match if r['matched'] == 'Y'}
sub_genomes = {r['species_code'] for r in ing}
yg_in = sorted((g for g in yg_genomes if g in sub_genomes))
check(len(yg_genomes) == 20 and len(yg_in) == 20, 'YGOB: %d of %d matched YGOB genomes have sequences in the subset' % (len(yg_in), len(yg_genomes)))

def is_wgd(r):
    return 'WGD' in r['layer'].split('+')

def species_stats(rows):
    return len({r['species_code'] for r in rows})
C = {}
C['rule'] = 'INGROUP & order==Saccharomycetales  +  all OUT_OG658  (Dikarya and other orders dropped)'
C['manifest_sha256'] = MANIFEST_SHA
C['n_total'] = len(sub)
C['n_ingroup'] = len(ing)
C['n_ingroup_genomes'] = species_stats(ing)
C['n_ingroup_species_keys'] = len({r['species_key'] for r in ing})
C['n_outgroup_og658'] = len(og)
C['n_dropped'] = len(drop)
C['dropped_by_role_order'] = {'%s|%s' % k: v for k, v in sorted(collections.Counter(((r['role'], r['order'] if r['role'] == 'INGROUP' else '-') for r in drop)).items())}
C['ingroup_by_family_seq'] = dict(collections.Counter((r['family'] for r in ing)))
C['ingroup_by_family_species'] = {f: len({r['species_code'] for r in ing if r['family'] == f}) for f in sorted({r['family'] for r in ing})}
C['ingroup_by_besthit_seq'] = dict(sorted(collections.Counter((r['best_hit'] for r in ing)).items()))
C['ingroup_genera'] = {g: dict(seq=n, genomes=len({r['species_code'] for r in ing if genus_of[r['species_code']] == g})) for g, n in sorted(genera.items())}
C['ingroup_by_layer_seq'] = dict(sorted(collections.Counter((r['layer'] for r in ing)).items()))
wgd_rows = [r for r in ing if is_wgd(r)]
non_rows = [r for r in ing if not is_wgd(r)]
C['wgd'] = dict(seq=len(wgd_rows), species=species_stats(wgd_rows))
C['non_wgd'] = dict(seq=len(non_rows), species=species_stats(non_rows))
C['wgd_by_family_seq'] = dict(collections.Counter((r['family'] for r in wgd_rows)))
C['non_wgd_by_family_seq'] = dict(collections.Counter((r['family'] for r in non_rows)))
C['wgd_by_besthit_seq'] = dict(sorted(collections.Counter((r['best_hit'] for r in wgd_rows)).items()))
C['non_wgd_by_besthit_seq'] = dict(sorted(collections.Counter((r['best_hit'] for r in non_rows)).items()))
C['seq_per_genome_summary'] = {'min': min(collections.Counter((r['species_code'] for r in ing)).values()), 'max': max(collections.Counter((r['species_code'] for r in ing)).values())}
C['outgroup_by_order'] = dict(sorted(collections.Counter((r['order'] for r in og)).items()))
C['outgroup_from_saccharomycetales_genomes'] = sorted((r['seq_id'] for r in og if r['order'] == SUBSET_ORDER))
full_by_fam = collections.Counter((r['family'] for r in man if r['role'] == 'INGROUP'))
C['full_ingroup_by_family_seq'] = dict(full_by_fam)
C['fraction_of_full_ingroup_kept'] = round(len(ing) / 1253, 4)
C['fraction_by_family_kept'] = {f: round(C['ingroup_by_family_seq'].get(f, 0) / n, 3) for f, n in full_by_fam.items()}
plan = read_tsv(os.path.join(P1B, 'outputs', 'set', 'sampling', 'sampling_plan.tsv'))
sel_sacch = [r for r in plan if r['selected'] == 'Y' and r['order'] == SUBSET_ORDER]
C['plan_selected_genomes_in_order'] = len(sel_sacch)
C['plan_selected_without_sequence'] = sorted((r['speciesCode'] for r in sel_sacch if r['speciesCode'] not in sub_genomes))
dias_rows = read_tsv(DIAS_MATCH)
dias_a1 = read_tsv(DIAS_A1)
dias_family = collections.defaultdict(set)
for r in dias_a1:
    dias_family[r['species'].replace(' (HYBRID)', '')].add(r['family'])
dias_out = []
for r in dias_rows:
    key, genome = (r['y1000_species_key'], r['y1000_genome'])
    g_rows = [m for m in man if m['role'] == 'INGROUP' and m['species_key'] == key] if key else []
    if r['matched'] != 'Y':
        status, order = ('NOT_MATCHED(' + r['note'] + ')', '')
    elif not g_rows:
        status, order = ('MATCHED_NO_PASS_SEQUENCE', order_of_tab.get(genome, ''))
    else:
        order = g_rows[0]['order']
        status = 'KEPT_IN_SUBSET' if order == SUBSET_ORDER else 'LOST_NOT_SACCHAROMYCETACEAE'
    dias_out.append(dict(dias_species=r['dias_species'], y1000_species_key=key, y1000_order=order, dias_family=';'.join(sorted(dias_family.get(r['dias_species'].replace(' (HYBRID)', ''), []))), n_seq_in_full_set=len(g_rows), n_seq_in_subset=len(g_rows) if status == 'KEPT_IN_SUBSET' else 0, status=status))
write_tsv(os.path.join(OUT, 'subset_dias_species.tsv'), dias_out)
C['dias_species_listed'] = len(dias_out)
C['dias_status'] = dict(collections.Counter((d['status'] for d in dias_out)))
C['dias_kept_species'] = [d['dias_species'] for d in dias_out if d['status'] == 'KEPT_IN_SUBSET']
C['dias_lost_species_orders'] = dict(collections.Counter((d['y1000_order'] for d in dias_out if d['status'] == 'LOST_NOT_SACCHAROMYCETACEAE')))
C['dias_kept_sequences'] = sum((d['n_seq_in_subset'] for d in dias_out))
C['dias_family_label_consistency'] = {d['dias_species']: (d['dias_family'], d['y1000_order']) for d in dias_out if (d['dias_family'] == 'Saccharomycetaceae') != (d['y1000_order'] == SUBSET_ORDER) and d['status'] != 'NOT_MATCHED(hybrid, not matched)'}
C['s288c_retained'] = {k: v in ids for k, v in S288C.items()}
C['ygob_genomes_retained'] = len(yg_in)
C['ygob_post_wgd_genomes'] = sorted((r['y1000_genome'] for r in yg_match if r['group'] == 'post-WGD'))
C['ygob_mapped_sequences'] = len(mapped)
desc = {}
for aln, path in G1_LABELS.items():
    lab = read_tsv(path)
    d = {}
    for x in ('PDR5', 'PDR15', 'PDR10', 'SNQ2', 'PDR18', 'PDR12', 'PDR11', 'AUS1', 'unassigned'):
        mem = [l for l in lab if l['tree_clade'] == x and l['role'] == 'INGROUP']
        mem_sub = [l for l in mem if l['seq_id'] in ids]
        d[x] = dict(ingroup_in_full_clade=len(mem), of_which_in_subset=len(mem_sub), species_in_subset=len({l['species'] for l in mem_sub}))
    desc[aln] = d
C['descriptive_full_fast_tree_clade_members_in_subset(not_subset_trees)'] = desc
with open(os.path.join(OUT, 'subset_counts.json'), 'w', encoding='utf-8') as f:
    json.dump(C, f, indent=1)
with open(os.path.join(OUT, 'subset_ids.txt'), 'w', encoding='utf-8', newline='') as f:
    for r in sub:
        f.write(r['seq_id'] + '\n')
write_tsv(os.path.join(OUT, 'subset_manifest.tsv'), sub, fields=list(man[0].keys()))
sp_rows = []
for code in sorted(sub_genomes):
    rs = [r for r in ing if r['species_code'] == code]
    sp_rows.append(dict(species_code=code, species_key=rs[0]['species_key'], genus=genus_of[code], family_2018=fam2018_of[code], wgd='Y' if is_wgd(rs[0]) else 'N', layer=rs[0]['layer'], n_seq=len(rs), n_PDR5like=sum((r['family'] == 'PDR5-like' for r in rs)), n_SNQ2like=sum((r['family'] == 'SNQ2-like' for r in rs)), n_PDR12like=sum((r['family'] == 'PDR12-like' for r in rs)), n_PDR11AUS1like=sum((r['family'] == 'PDR11/AUS1-like' for r in rs))))
write_tsv(os.path.join(OUT, 'subset_species.tsv'), sp_rows)
with open(os.path.join(OUT, 'subset_checks.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(checks) + '\n')
sums = []
for fn in ('subset_ids.txt', 'subset_manifest.tsv', 'subset_species.tsv', 'subset_dias_species.tsv', 'subset_counts.json', 'subset_checks.txt'):
    sums.append('%s  %s' % (sha256(os.path.join(OUT, fn)), fn))
sums.append('%s  scripts/c01_build_subset.py' % sha256(os.path.abspath(__file__)))
with open(os.path.join(OUT, 'SHA256SUMS.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(sums) + '\n')
print(json.dumps({k: C[k] for k in ('n_total', 'n_ingroup', 'n_ingroup_genomes', 'n_outgroup_og658', 'n_dropped', 'ingroup_by_family_seq', 'wgd', 'non_wgd', 'dias_status')}, indent=1))

import argparse, collections, csv, hashlib, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_tsv, write_tsv, log, SEED
ap = argparse.ArgumentParser()
ap.add_argument('--run', default='/home/pdr/work/phase1a/run_full')
ap.add_argument('--plan', required=True)
ap.add_argument('--outdir', required=True)
ap.add_argument('--k-per-order', type=int, default=3)
ap.add_argument('--seed', type=int, default=SEED)
ap.add_argument('--ingroup-rule', choices=['genome_dedup', 'rep_only'], default='genome_dedup')
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

def H(code):
    return int(hashlib.sha256(f'{a.seed}|{code}'.encode()).hexdigest(), 16)

def read_fasta(path):
    seqs, name, buf = ({}, None, [])
    for l in open(path, encoding='utf-8'):
        l = l.rstrip('\n')
        if l.startswith('>'):
            if name is not None:
                seqs[name] = ''.join(buf)
            name = l[1:].split()[0]
            buf = []
        else:
            buf.append(l.strip())
    if name is not None:
        seqs[name] = ''.join(buf)
    return seqs
plan = {r['speciesCode']: r for r in read_tsv(a.plan)}
sel = {c for c, r in plan.items() if r['selected'] == 'Y'}
ret = read_tsv(os.path.join(a.run, 'qc/retrieval_table.tsv'))
final = read_fasta(os.path.join(a.run, 'selected_final.faa'))
log('selected genomes', len(sel), 'final sequences', len(final))
rows = []
clus = {}
for c in csv.reader(open(os.path.join(a.run, 'qc/cdhit99_clusters.tsv'), encoding='utf-8'), delimiter='\t'):
    if len(c) >= 10 and '|' in c[0]:
        clus[c[0]] = (c[7], c[8])

def keep_ingroup(r):
    if r['cdhit99_rep_within_species'] == 'Y':
        return True
    if a.ingroup_rule == 'rep_only':
        return False
    rep, flag = clus[r['seq_id']]
    return rep.split('|')[0] != r['seq_id'].split('|')[0]
n_kept_other_strain = 0
for r in ret:
    if r['genome_kind'] != 'y1000_genome' or r['species_code'] not in sel:
        continue
    if r['passes_rule_final_seq'] == 'Y' and r['retrieval_class'] == 'PASS' and keep_ingroup(r):
        if r['cdhit99_rep_within_species'] != 'Y':
            n_kept_other_strain += 1
        rows.append(dict(seq_id=r['seq_id'], role='INGROUP', species_code=r['species_code'], species_key=r['species_key'], order=r['order'], layer=plan[r['species_code']]['layer'], source='Y1000+', length=r['length_final'], ctg_call=r['ctg_call'], cug_code=r['CUG_code'], retrans_status=r['retrans_status'], og=r['OG_id'], qc_flags=r['qc_flags'], pf14510=r['PF14510_present'], seq_sha256_phase1a=r['seq_sha256']))
n_in = len(rows)
cand = collections.defaultdict(list)
for r in ret:
    if r['genome_kind'] == 'y1000_genome' and r['species_code'] in sel and (r['OG_id'] == 'OG0000658') and (r['final_class'] == 'NEAR_MISS:NO_PF06422') and (int(r['n_PF00005']) >= 2) and (int(r['n_PF01061']) >= 2) and (1100 <= int(r['length_final']) <= 1800):
        cand[r['species_code']].append(r)
by_order = collections.defaultdict(list)
for sp in cand:
    by_order[cand[sp][0]['order']].append(sp)
chosen_g = {}
for o, sps in sorted(by_order.items()):
    for sp in sorted(sps, key=H)[:a.k_per_order]:
        chosen_g[sp] = 'ORDER_RULE'
anchor_g = 'Saccharomyces_cerevisiae'
if anchor_g in cand and anchor_g not in chosen_g:
    chosen_g[anchor_g] = 'ANCHOR(YOL075C)'
elif anchor_g in cand:
    chosen_g[anchor_g] = 'ORDER_RULE+ANCHOR(YOL075C)'
for sp, why in sorted(chosen_g.items()):
    r = sorted(cand[sp], key=lambda x: (-int(x['length_final']), x['seq_id']))[0]
    rows.append(dict(seq_id=r['seq_id'], role='OUT_OG658', species_code=sp, species_key=r['species_key'], order=r['order'], layer=why, source='Y1000+', length=r['length_final'], ctg_call=r['ctg_call'], cug_code=r['CUG_code'], retrans_status=r['retrans_status'], og=r['OG_id'], qc_flags=r['qc_flags'], pf14510=r['PF14510_present'], seq_sha256_phase1a=r['seq_sha256']))
n_og = len(rows) - n_in
dseq = read_fasta(os.path.join(a.run, 'outgroups/outgroups.faa'))
dik = []
for r in read_tsv(os.path.join(a.run, 'outgroups/outgroup_hits.tsv')):
    if r['source'] == 'UniProt' and r['passes_rule'] == 'Y' and (not r['organism'].startswith('Candida_albicans')):
        sid = f"{r['organism']}_{r['proteome']}|{r['accession']}"
        assert sid in dseq, sid
        dik.append(dict(seq_id=sid, role='OUT_DIKARYA', species_code=r['organism'], species_key=r['organism'], order='Dikarya_outgroup', layer='UNIPROT_PROTEOME', source='UniProt ' + r['proteome'], length=str(len(dseq[sid])), ctg_call='', cug_code='', retrans_status='', og='', qc_flags=r['rule_flags'], pf14510='Y' if int(r['n_PF14510'] or 0) else 'N', seq_sha256_phase1a=''))
rows.extend(dik)
log('ingroup', n_in, '(of which kept although the within-species representative is in another strain:', n_kept_other_strain, ') OG658 outgroup', n_og, 'Dikarya', len(dik))
seqs = {}
for r in rows:
    s = final[r['seq_id']] if r['role'] != 'OUT_DIKARYA' else dseq[r['seq_id']]
    assert len(s) == int(r['length']), (r['seq_id'], len(s), r['length'])
    seqs[r['seq_id']] = s
    r['seq_sha256'] = hashlib.sha256(s.encode()).hexdigest()
    if r['seq_sha256_phase1a']:
        assert r['seq_sha256'] == r['seq_sha256_phase1a'], ('sha256 differs from the Phase 1a retrieval table', r['seq_id'])
ids = [r['seq_id'] for r in rows]
assert len(set(ids)) == len(ids)
bad = [i for i in ids if re.search('[\\s():;,\\[\\]\\\'"]', i)]
assert not bad, bad
with open(os.path.join(a.outdir, 'analysis_set.faa'), 'w') as f:
    for r in sorted(rows, key=lambda x: ({'INGROUP': 0, 'OUT_OG658': 1, 'OUT_DIKARYA': 2}[x['role']], x['seq_id'])):
        f.write('>%s\n' % r['seq_id'])
        s = seqs[r['seq_id']]
        for i in range(0, len(s), 80):
            f.write(s[i:i + 80] + '\n')
write_tsv(os.path.join(a.outdir, 'analysis_set_manifest.tsv'), sorted(rows, key=lambda x: (x['role'], x['order'], x['seq_id'])))
cnt = collections.Counter(((r['role'], r['order']) for r in rows))
tot = collections.Counter((r['role'] for r in rows))
gen = {role: len({r['species_code'] for r in rows if r['role'] == role}) for role in tot}
lay = collections.Counter()
for r in rows:
    if r['role'] == 'INGROUP':
        for l in (r['layer'] or 'NONE').split('+'):
            lay[l] += 1
sp_layers = collections.defaultdict(set)
for r in rows:
    if r['role'] == 'INGROUP':
        sp_layers[r['species_code']].update((r['layer'] or 'NONE').split('+'))
lay_sp = collections.Counter((l for v in sp_layers.values() for l in v))
summ = dict(ingroup_rule=a.ingroup_rule, ingroup_kept_other_strain_collapsed=n_kept_other_strain, seed=a.seed, k_per_order=a.k_per_order, n_total=len(rows), by_role=dict(tot), genomes_by_role=gen, ingroup_by_order={o: cnt['INGROUP', o] for o in sorted({r['order'] for r in rows if r['role'] == 'INGROUP'})}, out_og658_by_order={o: cnt['OUT_OG658', o] for o in sorted({r['order'] for r in rows if r['role'] == 'OUT_OG658'})}, ingroup_sequences_by_layer=dict(lay), ingroup_species_by_layer=dict(lay_sp), n_dikarya_by_organism=dict(collections.Counter((r['species_code'] for r in dik))))
json.dump(summ, open(os.path.join(a.outdir, 'analysis_set_summary.json'), 'w'), indent=1)
print(json.dumps(summ, indent=1))

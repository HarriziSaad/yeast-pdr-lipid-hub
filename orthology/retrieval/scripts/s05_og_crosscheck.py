import argparse, collections, json, os
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from common import read_tsv, write_tsv, log
ap = argparse.ArgumentParser()
ap.add_argument('--orthogroups', required=True)
ap.add_argument('--species-table', required=True)
ap.add_argument('--summary', required=True)
ap.add_argument('--fasta', required=True)
ap.add_argument('--outdir', required=True)
ap.add_argument('--og', default='OG0000012')
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
st = read_tsv(a.species_table)
col2code = {}
for r in st:
    if r['kind'] == 'y1000_genome':
        col2code[r['pep_stem'] + '.final'] = r['speciesCode']
    elif r['pep_stem'] in ('candida_albicans.cgd', 'saccharomyces_cerevisiae.sgd'):
        col2code[r['pep_stem']] = r['speciesCode']
code2order = {r['speciesCode']: r['order'] for r in st}
summ = read_tsv(a.summary)
interest = {}
for r in summ:
    if r['status'] == 'PASS' or r['status'].startswith('NEAR_MISS'):
        interest[r['seq_id']] = r
members = []
og_of = {}
n_cols_ignored = 0
with open(a.orthogroups, encoding='utf-8') as f:
    header = f.readline().rstrip('\n').split('\t')
    cols = {i: col2code[h] for i, h in enumerate(header) if h in col2code}
    n_cols_ignored = len(header) - 1 - len(cols)
    missing_cols = sorted(set(col2code) - set(header))
    log('columns mapped', len(cols), 'ignored (outgroups / genomes outside this input)', n_cols_ignored, 'species-table columns missing in OG file', len(missing_cols))
    for line in f:
        tab = line.rstrip('\n').split('\t')
        og = tab[0]
        is_target = og == a.og
        for i, code in cols.items():
            if i >= len(tab) or not tab[i]:
                continue
            for g in tab[i].split(', '):
                sid = code + '|' + g
                if is_target:
                    members.append(sid)
                if sid in interest:
                    og_of[sid] = og
log('OG', a.og, 'members in mapped columns', len(members), 'interest sequences assigned to an OG', len(og_of), 'of', len(interest))
member_set = set(members)
need = member_set
length = {}
with open(a.fasta, encoding='utf-8') as f:
    sid = None
    for line in f:
        if line.startswith('>'):
            sid = line[1:].strip()
            if sid not in need:
                sid = None
        elif sid is not None:
            length[sid] = len(line.strip())
            sid = None
mem_rows = []
for s in members:
    code, _, g = s.partition('|')
    mem_rows.append(dict(seq_id=s, species_code=code, gene_id=g, in_search_fasta='Y' if s in length else 'N', seq_len=length.get(s, ''), summary_status=next((r['status'] for r in [interest.get(s)] if r), None) or ''))
write_tsv(os.path.join(a.outdir, 'og0000012_members.tsv'), mem_rows)
write_tsv(os.path.join(a.outdir, 'og_of_interest.tsv'), [dict(seq_id=k, OG=v) for k, v in sorted(og_of.items())])
by_id = {r['seq_id']: r for r in summ}
out = []
for r in summ:
    if r['status'] == 'PASS' or r['status'].startswith('NEAR_MISS'):
        s = r['seq_id']
        inog = 'Y' if s in member_set else 'N'
        r2 = dict(r)
        r2['in_OG0000012'] = inog
        r2['OG_id'] = og_of.get(s, '')
        r2['retrieval_class'] = r['status']
        out.append(r2)
n_og_nodomain = n_og_other = 0
for s in members:
    r = by_id.get(s)
    if r is not None and (r['status'] == 'PASS' or r['status'].startswith('NEAR_MISS')):
        continue
    code, _, g = s.partition('|')
    if r is None:
        n_og_nodomain += 1
        row = dict(seq_id=s, species_code=code, gene_id=g, seq_len=length.get(s, ''), n_PF00005=0, n_PF01061=0, n_PF06422=0, n_PF14510=0, n_PF19055=0, arch='', arch_AG='', PF14510_present='N', status='OG_FAIL:NO_DOMAIN_HIT', near_miss_reasons='no included Pfam domain', repair_candidate='')
    else:
        n_og_other += 1
        row = dict(r)
        row['status'] = 'OG_FAIL:OTHER_ABC'
    row['in_OG0000012'] = 'Y'
    row['OG_id'] = a.og
    row['retrieval_class'] = row['status']
    out.append(row)
write_tsv(os.path.join(a.outdir, 'classification_full.tsv'), out)
cnt = collections.Counter()
by_order = collections.defaultdict(collections.Counter)
for r in out:
    cls = 'PASS' if r['retrieval_class'] == 'PASS' else 'NEAR_MISS' if r['retrieval_class'].startswith('NEAR_MISS') else 'OG_FAIL'
    key = (cls, r['in_OG0000012'])
    cnt[key] += 1
    by_order[code2order.get(r['species_code'], '')][key] += 1
conc = dict(OG=a.og, n_members_mapped=len(members), n_members_unique=len(member_set), PASS_in_OG=cnt['PASS', 'Y'], PASS_not_in_OG=cnt['PASS', 'N'], OG_not_PASS=len(members) - cnt['PASS', 'Y'], OG_not_PASS_breakdown=dict(NEAR_MISS=cnt['NEAR_MISS', 'Y'], OTHER_ABC_domains_only=n_og_other, no_included_domain=n_og_nodomain), NEAR_MISS_not_in_OG=cnt['NEAR_MISS', 'N'], PASS_not_in_OG_by_OG=dict(collections.Counter((r['OG_id'] or 'unassigned/none' for r in out if r['retrieval_class'] == 'PASS' and r['in_OG0000012'] == 'N'))), by_order={o: {f'{k[0]}|inOG={k[1]}': v for k, v in c.items()} for o, c in by_order.items()}, columns_mapped=len(cols), columns_ignored=n_cols_ignored)
json.dump(conc, open(os.path.join(a.outdir, 'og_concordance.json'), 'w'), indent=1)
log(json.dumps({k: v for k, v in conc.items() if k != 'by_order'}))

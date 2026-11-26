import csv, json, collections, statistics, sys, os
root = sys.argv[1]
O = lambda f: os.path.join(root, 'outputs', f)
rd = lambda f: list(csv.DictReader(open(O(f), encoding='utf-8'), delimiter='\t'))
mp = rd('site_class_map.tsv')
ctx = rd('site_class_map_context_outgroups.tsv')
exp = rd('site_class_map_Cdr1_experimental.tsv')
stats = {}

def distinct(rows, name):
    d = {}
    for r in rows:
        if r['paralog'] == name:
            d[int(r['P33302_residue'])] = r
    return d
p5 = distinct(mp, 'Pdr5')
n = len(p5)
cnt = collections.Counter()
for k, r in p5.items():
    for col in ('paralog_residue', 'residue_bestref', 'residue_domainwise', 'residue_seqaln', 'consensus_residue'):
        if r[col] != '' and int(r[col]) == k:
            cnt[col] += 1
    votes = collections.Counter((x for x in (r['residue_domainwise'], r['paralog_residue'], r['residue_bestref']) if x != ''))
    top = max(votes.values())
    pick = next((c for c in (r['residue_domainwise'], r['paralog_residue'], r['residue_bestref']) if c != '' and votes[c] == top))
    if int(pick) == k:
        cnt['structure_only_consensus'] += 1
stats['pdr5_selftest_n'] = n
stats['pdr5_selftest_correct'] = dict(cnt)
cdr = distinct(ctx, 'CaCdr1_AF')
e9 = distinct(exp, '9VT1')
agree = collections.Counter()
tot = 0
for k in e9:
    if e9[k]['paralog_residue'] == '' or cdr[k]['paralog_residue'] == '':
        continue
    tot += 1
    ref = e9[k]['paralog_residue']
    for col in ('paralog_residue', 'residue_bestref', 'residue_domainwise', 'residue_seqaln', 'consensus_residue'):
        if cdr[k][col] == ref:
            agree[col] += 1
    votes = collections.Counter((x for x in (cdr[k]['residue_domainwise'], cdr[k]['paralog_residue'], cdr[k]['residue_bestref']) if x != ''))
    if votes:
        top = max(votes.values())
        pick = next((c for c in (cdr[k]['residue_domainwise'], cdr[k]['paralog_residue'], cdr[k]['residue_bestref']) if c != '' and votes[c] == top))
        if pick == ref:
            agree['structure_only_consensus'] += 1
stats['cdr1_AF_vs_9VT1_n'] = tot
stats['cdr1_AF_vs_9VT1_agree'] = dict(agree)
dw = {}
for name in ['Pdr5', 'Snq2', 'Pdr10', 'Pdr11', 'Pdr12', 'Pdr15', 'Pdr18', 'Aus1']:
    d = distinct(mp, name)
    vals = [float(r['CA_dist_domainwise_A']) for r in d.values() if r['CA_dist_domainwise_A'] != '']
    prim = [float(r['CA_distance_A']) for r in d.values() if r['CA_distance_A'] != '']
    dw[name] = {'n_domainwise_pairs': len(vals), 'le3_domainwise': sum((v <= 3 for v in vals)), 'le3_primary': sum((v <= 3 for v in prim)), 'n_primary': len(prim)}
stats['CA_le3A'] = dw
fl = {}
for name in ['Pdr5', 'Snq2', 'Pdr10', 'Pdr11', 'Pdr12', 'Pdr15', 'Pdr18', 'Aus1']:
    d = distinct(mp, name)
    c = collections.Counter()
    for r in d.values():
        pl = r['pLDDT']
        if pl != '':
            c['pLDDT<70' if float(pl) < 70 else 'pLDDT70-80' if float(pl) < 80 else 'pLDDT>=80'] += 1
        if 'MAPPING_DIFFERS' in r['flags']:
            c['mapping_differs'] += 1
        if 'MAPPING_UNCERTAIN' in r['flags']:
            c['mapping_uncertain'] += 1
        if 'NOT_ALIGNABLE' in r['flags']:
            c['not_alignable'] += 1
        if r['sidechain_ok'] == 'Y':
            c['sidechain_ok'] += 1
    fl[name] = dict(c)
stats['distinct_residue_flags'] = fl
json.dump(stats, open(O('summary_stats.json'), 'w'), indent=1)
print(json.dumps(stats, indent=1))

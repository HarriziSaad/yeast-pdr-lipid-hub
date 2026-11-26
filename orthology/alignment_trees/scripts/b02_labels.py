import argparse, collections, json, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_tsv, write_tsv, log
ap = argparse.ArgumentParser()
ap.add_argument('--setdir', required=True)
ap.add_argument('--seeds', required=True)
ap.add_argument('--threads', type=int, default=8)
a = ap.parse_args()
FAM = {'PDR5': 'PDR5-like', 'PDR15': 'PDR5-like', 'PDR10': 'PDR5-like', 'SNQ2': 'SNQ2-like', 'PDR18': 'SNQ2-like', 'PDR12': 'PDR12-like', 'PDR11': 'PDR11/AUS1-like', 'AUS1': 'PDR11/AUS1-like'}
w = os.path.join(a.setdir, 'labels')
os.makedirs(w, exist_ok=True)
seeds = {}
name = None
for l in open(a.seeds):
    if l.startswith('>'):
        gn = [t for t in l.split() if t.startswith('GN=')][0][3:]
        name = gn
        seeds[name] = []
    else:
        seeds[name].append(l.strip())
with open(os.path.join(w, 'seeds.faa'), 'w') as f:
    for k, v in seeds.items():
        f.write('>%s\n%s\n' % (k, ''.join(v)))
subprocess.run(['diamond', 'makedb', '--in', os.path.join(w, 'seeds.faa'), '-d', os.path.join(w, 'seeds'), '--quiet'], check=True)
subprocess.run(['diamond', 'blastp', '-q', os.path.join(a.setdir, 'analysis_set.faa'), '-d', os.path.join(w, 'seeds'), '--very-sensitive', '-e', '1e-3', '--max-target-seqs', '8', '--max-hsps', '1', '-p', str(a.threads), '--outfmt', '6', 'qseqid', 'sseqid', 'bitscore', 'qlen', 'length', '-o', os.path.join(w, 'hits.tsv'), '--quiet'], check=True)
best = collections.defaultdict(dict)
for l in open(os.path.join(w, 'hits.tsv')):
    q, s, b, ql, ln = l.rstrip('\n').split('\t')
    best[q][s] = max(best[q].get(s, 0), float(b))
man = read_tsv(os.path.join(a.setdir, 'analysis_set_manifest.tsv'))
for r in man:
    h = sorted(best.get(r['seq_id'], {}).items(), key=lambda x: -x[1])
    r['best_hit'] = h[0][0] if h else 'NO_HIT'
    r['best_bits'] = '%.0f' % h[0][1] if h else ''
    r['second_hit'] = h[1][0] if len(h) > 1 else ''
    r['second_bits'] = '%.0f' % h[1][1] if len(h) > 1 else ''
    r['family'] = FAM.get(r['best_hit'], 'NO_HIT')
write_tsv(os.path.join(a.setdir, 'analysis_set_manifest.tsv'), man)
ing = [r for r in man if r['role'] == 'INGROUP']
tab = collections.Counter(((r['order'], r['family']) for r in ing))
orders = sorted({r['order'] for r in ing})
fams = ['PDR5-like', 'SNQ2-like', 'PDR12-like', 'PDR11/AUS1-like', 'NO_HIT']
rows = [dict(order=o, **{f: tab[o, f] for f in fams}, total=sum((tab[o, f] for f in fams))) for o in orders]
rows.append(dict(order='TOTAL', **{f: sum((tab[o, f] for o in orders)) for f in fams}, total=len(ing)))
write_tsv(os.path.join(a.setdir, 'counts_order_x_family.tsv'), rows)
bh = collections.Counter((r['best_hit'] for r in ing))
write_tsv(os.path.join(a.setdir, 'counts_besthit_label.tsv'), [dict(best_hit=k, ingroup_sequences=v) for k, v in sorted(bh.items())])
spf = collections.defaultdict(set)
for r in ing:
    spf[r['best_hit']].add(r['species_code'])
write_tsv(os.path.join(a.setdir, 'counts_besthit_species.tsv'), [dict(best_hit=k, sequences=bh[k], species=len(spf[k])) for k in sorted(bh)])
for r in rows:
    log(r)
log(dict(bh))
for role in ('OUT_OG658', 'OUT_DIKARYA'):
    log(role, dict(collections.Counter((r['best_hit'] for r in man if r['role'] == role))))

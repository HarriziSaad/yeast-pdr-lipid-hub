import argparse, collections, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_tsv, write_tsv, log
ap = argparse.ArgumentParser()
ap.add_argument('--ygob-dir', default='/mnt/c/Users/MSI/Documents/Advisor/project/phase0_verify/data_access/downloads/ygob')
ap.add_argument('--setdir', default='/home/pdr/work/phase1b/set')
ap.add_argument('--ygob-match', default='/home/pdr/work/phase1b/set/sampling/ygob_match.tsv')
ap.add_argument('--out', default='/home/pdr/work/phase1b/set/ygob')
ap.add_argument('--threads', type=int, default=4)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
POST1 = ['Vanderwaltozyma polyspora', 'Tetrapisispora phaffii', 'Tetrapisispora blattae', 'Naumovozyma dairenensis', 'Naumovozyma castellii', 'Kazachstania naganishii', 'Kazachstania africana', 'Candida glabrata (= Nakaseomyces glabratus)', 'Saccharomyces bayanus var. uvarum (= S. uvarum)', 'Saccharomyces kudriavzevii', 'Saccharomyces mikatae', 'Saccharomyces cerevisiae']
PRE = ['Zygosaccharomyces rouxii', 'Torulaspora delbrueckii', 'Kluyveromyces lactis', 'Eremothecium gossypii', 'Eremothecium cymbalariae', 'Lachancea kluyveri', 'Lachancea thermotolerans', 'Lachancea waltii']
POST2 = list(reversed(POST1))
cols = [(s, 'post1') for s in POST1] + [('ANC', 'anc')] + [(s, 'pre') for s in PRE] + [(s, 'post2') for s in POST2]
assert len(cols) == 33
PILLARS = {'Anc_5.500': 'PDR5/PDR15', 'Anc_7.123a': 'PDR11/AUS1', 'Anc_3.220': 'SNQ2', 'Anc_8.517': 'PDR12'}
gm = {r['ygob_species']: r for r in read_tsv(a.ygob_match)}
rows = []
for l in open(os.path.join(a.ygob_dir, 'Pillars.tab'), encoding='utf-8'):
    f = l.rstrip('\n').split('\t')
    if len(f) != 33 or f[12] not in PILLARS:
        continue
    for (sp, pos), g in zip(cols, f):
        if g == '---' or pos == 'anc':
            continue
        rows.append(dict(pillar=f[12], family=PILLARS[f[12]], ygob_species=sp, position=pos, ygob_gene=g, group=gm[sp]['group'], y1000_genome=gm[sp]['y1000_genome']))
log('pillar members', len(rows))
aa = {}
nm = None
for l in open(os.path.join(a.ygob_dir, 'AA.fsa'), encoding='utf-8'):
    if l.startswith('>'):
        nm = l[1:].split()[0]
        aa[nm] = []
    else:
        aa[nm].append(l.strip().rstrip('*'))
need = {r['ygob_gene'] for r in rows}
missing = [g for g in need if g not in aa]
log('genes needed', len(need), 'missing in AA.fsa', len(missing))
with open(os.path.join(a.out, 'ygob_pdr_genes.faa'), 'w') as f:
    for g in sorted(need - set(missing)):
        f.write('>%s\n%s\n' % (g, ''.join(aa[g])))
man = read_tsv(os.path.join(a.setdir, 'analysis_set_manifest.tsv'))
ygenomes = {r['y1000_genome'] for r in gm.values()}
seqs = {}
nm = None
for l in open(os.path.join(a.setdir, 'analysis_set.faa')):
    if l.startswith('>'):
        nm = l[1:].strip()
        seqs[nm] = []
    else:
        seqs[nm].append(l.strip())
sub = [r['seq_id'] for r in man if r['species_code'] in ygenomes]
with open(os.path.join(a.out, 'ygob_genomes_set.faa'), 'w') as f:
    for s in sub:
        f.write('>%s\n%s\n' % (s, ''.join(seqs[s])))
log('subject sequences (YGOB genomes in the analysis set)', len(sub))
subprocess.run(['diamond', 'makedb', '--in', os.path.join(a.out, 'ygob_genomes_set.faa'), '-d', os.path.join(a.out, 'ygobdb'), '--quiet'], check=True)
subprocess.run(['diamond', 'blastp', '-q', os.path.join(a.out, 'ygob_pdr_genes.faa'), '-d', os.path.join(a.out, 'ygobdb'), '--very-sensitive', '-e', '1e-10', '-k', '10', '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qstart', 'qend', 'bitscore', '-p', str(a.threads), '-o', os.path.join(a.out, 'ygob_hits.tsv'), '--quiet'], check=True)
hits = collections.defaultdict(list)
for l in open(os.path.join(a.out, 'ygob_hits.tsv')):
    q, s, pid, ln, ql, sl, qs, qe, bs = l.rstrip('\n').split('\t')
    hits[q].append((float(bs), s, float(pid), int(ln), int(ql), int(sl), int(qs), int(qe)))
for r in rows:
    g = r['ygob_gene']
    best = None
    for bs, s, pid, ln, ql, sl, qs, qe in sorted(hits.get(g, []), reverse=True):
        if s.split('|')[0] == r['y1000_genome'] or s.startswith(r['y1000_genome'] + '|'):
            best = (s, pid, (qe - qs + 1) / ql, ql, sl, bs)
            break
    if best is None:
        r.update(mapped_seq='', pident='', qcov='', ygob_len='', set_len='', status='NO_HIT_IN_SET')
    else:
        ok = best[1] >= 85.0 and best[2] >= 0.8
        r.update(mapped_seq=best[0] if ok else '', pident=best[1], qcov=round(best[2], 3), ygob_len=best[3], set_len=best[4], status='MAPPED' if ok else 'WEAK_HIT')
write_tsv(os.path.join(a.out, 'ygob_map.tsv'), rows)
cnt = collections.Counter(((r['family'], r['group'], r['status']) for r in rows))
for k, v in sorted(cnt.items()):
    log(k, v)
tgt = collections.Counter((r['mapped_seq'] for r in rows if r['mapped_seq']))
log('set sequences hit by more than one YGOB gene:', {k: v for k, v in tgt.items() if v > 1})

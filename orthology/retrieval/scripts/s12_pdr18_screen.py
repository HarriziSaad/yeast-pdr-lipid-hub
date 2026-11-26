import argparse, collections, os, subprocess, tempfile
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from common import read_tsv, write_tsv, fasta_iter, write_fasta, log
ap = argparse.ArgumentParser()
ap.add_argument('--retrieval', required=True)
ap.add_argument('--final-faa', required=True)
ap.add_argument('--seeds', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--genera', default='Saccharomyces')
ap.add_argument('--threads', type=int, default=4)
a = ap.parse_args()
gen = set(a.genera.split(','))
rows = [r for r in read_tsv(a.retrieval) if r['retrieval_class'] == 'PASS' and r['genome_kind'] == 'y1000_genome' and (r['species_key'].split('_')[0] in gen)]
ids = {r['seq_id'] for r in rows}
seqs = {sid: s for sid, _, s in fasta_iter(a.final_faa) if sid in ids}
seed_name = {'P33302': 'PDR5', 'P32568': 'SNQ2', 'P51533': 'PDR10', 'P40550': 'PDR11', 'Q02785': 'PDR12', 'Q04182': 'PDR15', 'P53756': 'PDR18', 'Q08409': 'AUS1'}
tmp = tempfile.mkdtemp(prefix='pdr18_')
with open(os.path.join(tmp, 'seeds.faa'), 'w') as fh:
    for sid, d, s in fasta_iter(a.seeds):
        acc = sid.split('|')[1] if '|' in sid else sid
        write_fasta(fh, acc, s)
with open(os.path.join(tmp, 'q.faa'), 'w') as fh:
    for sid, s in seqs.items():
        write_fasta(fh, sid, s)
subprocess.run(['diamond', 'makedb', '--in', os.path.join(tmp, 'seeds.faa'), '-d', os.path.join(tmp, 'seeds')], check=True, capture_output=True)
p = subprocess.run(['diamond', 'blastp', '-d', os.path.join(tmp, 'seeds'), '-q', os.path.join(tmp, 'q.faa'), '--very-sensitive', '--max-target-seqs', '8', '--outfmt', '6', 'qseqid', 'sseqid', 'bitscore', '-p', str(a.threads), '--evalue', '1e-5'], check=True, capture_output=True, text=True)
hits = collections.defaultdict(list)
for line in p.stdout.splitlines():
    q, s, b = line.split('\t')
    hits[q].append((float(b), s))
out = []
for r in sorted(rows, key=lambda r: (r['species_key'], r['gene_id'])):
    h = sorted(hits.get(r['seq_id'], []), reverse=True)
    lab = seed_name.get(h[0][1], h[0][1]) if h else 'no_hit'
    tie = 'Y' if len(h) > 1 and h[0][0] == h[1][0] else 'N'
    out.append(dict(species_code=r['species_code'], species_key=r['species_key'], gene_id=r['gene_id'], length=r['length_final'], tail_after_last_PF01061=r['tail_after_last_PF01061'], PF14510_present=r['PF14510_present'], nearest_S288C_paralog_label=lab, label_tie=tie, candidate_intact_Pdr18='Y' if lab == 'PDR18' and r['tail_after_last_PF01061'] and (int(r['tail_after_last_PF01061']) >= 100) else 'N' if lab == 'PDR18' else '', qc_flags=r['qc_flags']))
write_tsv(a.out, out)
c = collections.Counter((r['nearest_S288C_paralog_label'] for r in out))
pd = [r for r in out if r['nearest_S288C_paralog_label'] == 'PDR18']
log('Saccharomyces PASS sequences', len(out), dict(c), '| labelled PDR18:', len(pd), 'with tail>=100:', sum((1 for r in pd if r['candidate_intact_Pdr18'] == 'Y')))

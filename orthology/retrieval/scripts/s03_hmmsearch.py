import argparse, json, os, sys, time
import pyhmmer
from pyhmmer import easel, plan7
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from common import log

def _s(x):
    return x.decode() if isinstance(x, (bytes, bytearray)) else '' if x is None else str(x)
HMMS = ['PF00005', 'PF01061', 'PF06422', 'PF14510', 'PF19055']
ap = argparse.ArgumentParser()
ap.add_argument('--fasta', required=True)
ap.add_argument('--hmm-dir', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--threads', type=int, default=10)
ap.add_argument('--chunk', type=int, default=100000, help='sequences per block (memory control only; Z is fixed to the total)')
a = ap.parse_args()
hmms = []
meta = {}
for acc in HMMS:
    with plan7.HMMFile(os.path.join(a.hmm_dir, acc + '.hmm')) as hf:
        h = hf.read()
    hmms.append(h)
    meta[acc] = dict(name=_s(h.name), accession=_s(h.accession), length=h.M, GA=list(h.cutoffs.gathering) if h.cutoffs.gathering_available() else None)
    log('HMM', acc, meta[acc])
    assert meta[acc]['GA'] is not None, acc
t0 = time.time()
with open(a.fasta, 'rb') as f:
    n_total = sum((1 for line in f if line.startswith(b'>')))
log('sequences in', a.fasta, n_total)
rows_out = 0
n_seen = 0
hit_seqs = set()
hdr = ['seq_id', 'species_code', 'gene_id', 'seq_len', 'hmm_acc', 'hmm_name', 'dom_index', 'n_dom_in_hit', 'env_from', 'env_to', 'ali_from', 'ali_to', 'hmm_from', 'hmm_to', 'hmm_len', 'dom_score', 'dom_cEvalue', 'dom_iEvalue', 'seq_score']
with open(a.out, 'w', encoding='utf-8', newline='\n') as out, easel.SequenceFile(a.fasta, digital=True, alphabet=easel.Alphabet.amino()) as sf:
    out.write('\t'.join(hdr) + '\n')
    while True:
        block = sf.read_block(sequences=a.chunk)
        if len(block) == 0:
            break
        n_seen += len(block)
        for tophits in pyhmmer.hmmsearch(hmms, block, cpus=a.threads, bit_cutoffs='gathering', Z=n_total):
            acc = _s(tophits.query.accession)
            acc = acc.split('.')[0]
            qname = _s(tophits.query.name)
            for hit in tophits:
                if not hit.included:
                    continue
                sid = _s(hit.name)
                code, _, gid = sid.partition('|')
                doms = [d for d in hit.domains if d.included]
                for i, d in enumerate(doms, 1):
                    al = d.alignment
                    out.write('\t'.join(map(str, [sid, code, gid, hit.length, acc, qname, i, len(doms), d.env_from, d.env_to, al.target_from, al.target_to, al.hmm_from, al.hmm_to, al.hmm_length, round(d.score, 2), '%.3g' % d.c_evalue, '%.3g' % d.i_evalue, round(hit.score, 2)])) + '\n')
                    rows_out += 1
                if doms:
                    hit_seqs.add(sid)
        log('searched', n_seen, 'of', n_total, 'domain rows', rows_out, 'seqs with hits', len(hit_seqs), '%.0fs' % (time.time() - t0))
run = dict(tool='pyhmmer', pyhmmer_version=pyhmmer.__version__, hmmer_version_of_pyhmmer_build='see tools/software_versions.tsv', options=dict(bit_cutoffs='gathering', Z=n_total), hmms=meta, n_sequences=n_total, n_sequences_with_included_hit=len(hit_seqs), n_domain_rows=rows_out, fasta=a.fasta, threads=a.threads, chunk=a.chunk, wall_seconds=round(time.time() - t0, 1))
json.dump(run, open(a.out + '.run.json', 'w'), indent=1)
log('done', run['n_sequences_with_included_hit'], 'sequences with >=1 included domain')

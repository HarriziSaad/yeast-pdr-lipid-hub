import argparse, collections, os
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from common import read_tsv, write_tsv, log
ap = argparse.ArgumentParser()
ap.add_argument('--merged', required=True)
ap.add_argument('--extras-dir', default='')
ap.add_argument('--species-table', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--summary', required=True)
a = ap.parse_args()
st = read_tsv(a.species_table)
genomes = {r['speciesCode'] for r in st if r['kind'] == 'y1000_genome'}
extras = {r['pep_file']: r['speciesCode'] for r in st if r['kind'] == 'reference_annotation'}
stats = collections.defaultdict(lambda: dict(n=0, res=0, n_internal_stop=0, n_X=0, n_empty=0))
seen = set()
n_skip_empty_code = n_bad = 0
n_unknown_code = collections.Counter()

def add(code, gid, seq, out):
    s = stats[code]
    s['n'] += 1
    s['res'] += len(seq)
    body = seq[:-1] if seq.endswith('*') else seq
    if '*' in body:
        s['n_internal_stop'] += 1
    if 'X' in body.upper():
        s['n_X'] += 1
    if not seq:
        s['n_empty'] += 1
    out.write('>' + code + '|' + gid + '\n' + body + '\n')
with open(a.out, 'w', encoding='utf-8', newline='\n') as out:
    with open(a.merged, encoding='utf-8') as f:
        hid, buf = (None, [])

        def flush():
            global n_skip_empty_code, n_bad
            if hid is None:
                return
            code, _, gid = hid.partition('|')
            seq = ''.join(buf)
            if code == '':
                n_skip_empty_code += 1
                return
            if code not in genomes:
                n_unknown_code[code] += 1
                return
            full = code + '|' + gid
            if full in seen or '|' in gid or (not seq):
                n_bad += 1
                return
            seen.add(full)
            add(code, gid, seq, out)
        for line in f:
            if line.startswith('>'):
                flush()
                hid = line[1:].split(None, 1)[0].rstrip('\r\n') if len(line) > 1 else ''
                buf = []
            else:
                buf.append(line.strip())
        flush()
    n_extra = 0
    for fn, code in extras.items():
        p = os.path.join(a.extras_dir, fn) if a.extras_dir else ''
        if not p or not os.path.exists(p):
            log('extra file not found, skipped:', fn)
            continue
        with open(p, encoding='utf-8') as f:
            hid, buf = (None, [])
            recs = []
            for line in f:
                if line.startswith('>'):
                    if hid is not None:
                        recs.append((hid, ''.join(buf)))
                    hid = line[1:].split(None, 1)[0].rstrip('\r\n')
                    buf = []
                else:
                    buf.append(line.strip())
            if hid is not None:
                recs.append((hid, ''.join(buf)))
        for gid, seq in recs:
            if code + '|' + gid in seen or not seq or '|' in gid:
                n_bad += 1
                continue
            seen.add(code + '|' + gid)
            add(code, gid, seq, out)
            n_extra += 1
        log('extra added', fn, '->', code, len(recs))
log('records kept', len(seen), 'skipped empty-code records', n_skip_empty_code, 'bad/dup', n_bad, 'unknown codes', dict(n_unknown_code))
rows = [dict(speciesCode=k, **v) for k, v in sorted(stats.items())]
write_tsv(a.summary, rows, ['speciesCode', 'n', 'res', 'n_internal_stop', 'n_X', 'n_empty'])
missing = sorted(genomes - set(stats))
log('genomes in species table without proteins in the merged FASTA:', len(missing), missing[:5])

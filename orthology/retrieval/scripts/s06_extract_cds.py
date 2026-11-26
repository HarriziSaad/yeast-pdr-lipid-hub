import argparse, collections, gzip, io, os, re, tarfile, zlib
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from common import read_tsv, write_tsv, log
ap = argparse.ArgumentParser()
ap.add_argument('--ids', required=True, help='text file, one seq_id (speciesCode|gene) per line')
ap.add_argument('--species-table', required=True)
ap.add_argument('--cds-tar', default='')
ap.add_argument('--cds-merged', default='')
ap.add_argument('--out', required=True)
ap.add_argument('--log', required=True)
ap.add_argument('--allow-truncated', action='store_true')
a = ap.parse_args()
want = collections.defaultdict(set)
for line in open(a.ids, encoding='utf-8'):
    s = line.strip()
    if s:
        c, _, g = s.partition('|')
        want[c].add(g)
st = {r['speciesCode']: r for r in read_tsv(a.species_table)}
stem2code = {r['pep_stem']: c for c, r in st.items() if r['kind'] == 'y1000_genome'}
_norm = lambda s: re.sub('[^A-Za-z0-9]', '', s).lower()
_ncount = collections.Counter((_norm(s) for s in stem2code))
norm2code = {_norm(s): c for s, c in stem2code.items() if _ncount[_norm(s)] == 1}
fuzzy = {}
found = collections.Counter()
members_seen = set()
CDS_ID_RE = re.compile('^(.*)\\.(g\\d+\\.m\\d+)$')
MAX_BYTES = 400000000

def parse_cds(data):
    sid, buf = (None, [])
    for line in data.decode('ascii', errors='replace').split('\n'):
        if line.startswith('>'):
            if sid is not None:
                yield (sid, ''.join(buf))
            sid = line[1:].split(None, 1)[0].strip()
            buf = []
        else:
            buf.append(line.strip())
    if sid is not None:
        yield (sid, ''.join(buf))
with open(a.out, 'w', encoding='ascii', newline='\n') as out:
    if a.cds_tar:
        try:
            with tarfile.open(a.cds_tar, 'r|gz') as tf:
                for m in tf:
                    if not m.isfile():
                        continue
                    base = os.path.basename(m.name)
                    if not base.endswith('.final.cds.tar.gz'):
                        continue
                    stem = base[:-len('.final.cds.tar.gz')]
                    code = stem2code.get(stem)
                    if code is None and _norm(stem) in norm2code:
                        code = norm2code[_norm(stem)]
                        fuzzy[code] = stem
                    if code is None or code not in want:
                        continue
                    try:
                        raw = tf.extractfile(m).read()
                        with tarfile.open(fileobj=io.BytesIO(raw), mode='r:gz') as inner:
                            fl = [x for x in inner if x.isfile()]
                            if len(fl) != 1 or fl[0].size > MAX_BYTES:
                                raise ValueError('unexpected nested layout ' + m.name)
                            data = inner.extractfile(fl[0]).read()
                    except (EOFError, tarfile.ReadError, OSError, zlib.error) as e:
                        if not a.allow_truncated:
                            raise
                        log('TEST MODE: truncated member', m.name, repr(e))
                        break
                    members_seen.add(code)
                    for tok, seq in parse_cds(data):
                        mm = CDS_ID_RE.match(tok)
                        gid = mm.group(2) if mm else tok
                        if gid in want[code]:
                            out.write('>' + code + '|' + gid + '\n' + seq.upper() + '\n')
                            found[code] += 1
        except (EOFError, tarfile.ReadError, OSError, zlib.error) as e:
            if not a.allow_truncated:
                raise
            log('TEST MODE: archive ended early:', repr(e))
    elif a.cds_merged:
        with gzip.open(a.cds_merged, 'rt', encoding='ascii') as f:
            sid, buf, keep = (None, [], False)
            for line in f:
                if line.startswith('>'):
                    if keep:
                        out.write('>' + sid + '\n' + ''.join(buf) + '\n')
                        found[sid.split('|')[0]] += 1
                    sid = line[1:].split(None, 1)[0].strip()
                    c, _, g = sid.partition('|')
                    keep = g in want.get(c, ())
                    members_seen.add(c)
                    buf = []
                elif keep:
                    buf.append(line.strip())
            if keep:
                out.write('>' + sid + '\n' + ''.join(buf) + '\n')
                found[sid.split('|')[0]] += 1
    else:
        raise SystemExit('give --cds-tar or --cds-merged')
rows = []
for code in sorted(want):
    r = st.get(code, {})
    kind = r.get('kind', '')
    if kind != 'y1000_genome':
        status = 'NO_CDS_FOR_REFERENCE_ANNOTATION'
    elif code in members_seen:
        status = 'OK' if found[code] == len(want[code]) else 'PARTIAL'
    else:
        status = 'CDS_NOT_AVAILABLE'
    rows.append(dict(speciesCode=code, kind=kind, ctg_class=r.get('ctg_class', ''), n_requested=len(want[code]), n_found=found[code], status=status, cds_member_stem_if_not_exact=fuzzy.get(code, '')))
write_tsv(a.log, rows)
log('genomes requested', len(rows), 'with CDS found', sum((1 for r in rows if r['status'] in ('OK', 'PARTIAL'))), 'CDS records written', sum(found.values()), dict(collections.Counter((r['status'] for r in rows))))

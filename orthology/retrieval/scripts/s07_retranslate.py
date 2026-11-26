import argparse, collections
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from common import read_tsv, write_tsv, write_fasta, fasta_iter, retranslate_one, log
ap = argparse.ArgumentParser()
ap.add_argument('--classification', required=True)
ap.add_argument('--pep', required=True)
ap.add_argument('--cds', required=True)
ap.add_argument('--species-table', required=True)
ap.add_argument('--out-faa', required=True)
ap.add_argument('--out-log', required=True)
a = ap.parse_args()
sel = [r['seq_id'] for r in read_tsv(a.classification)]
selset = set(sel)
st = {r['speciesCode']: r for r in read_tsv(a.species_table)}
pep = {}
with open(a.pep, encoding='utf-8') as f:
    sid = None
    for line in f:
        if line.startswith('>'):
            s = line[1:].strip()
            sid = s if s in selset else None
        elif sid is not None:
            pep[sid] = line.strip()
            sid = None
cds = {}
for sid, _, seq in fasta_iter(a.cds):
    cds[sid] = seq
log('selected', len(sel), 'pep found', len(pep), 'cds found', len(cds))
rows, cnt = ([], collections.Counter())
with open(a.out_faa, 'w', encoding='utf-8', newline='\n') as out:
    for sid in sel:
        code = sid.split('|')[0]
        r = st.get(code, {})
        call = r.get('ctg_call', '')
        kind = r.get('kind', '')
        p = pep.get(sid)
        if p is None:
            rows.append(dict(seq_id=sid, status='PEP_NOT_FOUND', retranslated='N', ctg_call=call))
            cnt['PEP_NOT_FOUND'] += 1
            continue
        if kind != 'y1000_genome':
            res = dict(status='NA_REFERENCE_ANNOTATION', seq=p[:-1] if p.endswith('*') else p, retranslated='N', table='', n_ctg='', n_changed=0, n_mismatch='', note='')
        elif sid not in cds:
            res = dict(status='CDS_NOT_AVAILABLE', seq=p[:-1] if p.endswith('*') else p, retranslated='N', table='', n_ctg='', n_changed=0, n_mismatch='', note='')
        else:
            res = retranslate_one(cds[sid], p, call)
        cug_pending = 'Y' if res['status'] in ('CDS_NOT_AVAILABLE', 'MISMATCH_PEP_CDS', 'CDS_NOT_MULTIPLE_OF_3_LEGACY') and call in ('S', 'A', '?') and (kind == 'y1000_genome') else 'N'
        cnt[res['status']] += 1
        write_fasta(out, sid, res['seq'], 'retrans=%s status=%s ctg_call=%s' % (res['retranslated'], res['status'], call))
        rows.append(dict(seq_id=sid, ctg_call=call, ctg_class=r.get('ctg_class', ''), status=res['status'], retranslated=res['retranslated'], table=res['table'], n_inframe_CTG=res['n_ctg'], n_changed=res['n_changed'], n_mismatch=res['n_mismatch'], cug_pending=cug_pending, len_pep=len(p[:-1] if p.endswith('*') else p), len_final=len(res['seq']), note=res['note']))
write_tsv(a.out_log, rows)
log('status counts', dict(cnt))

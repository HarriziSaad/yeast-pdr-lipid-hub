import argparse, collections, glob, gzip, os, re
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from common import read_tsv, write_tsv, log
ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest='cmd', required=True)
p1 = sub.add_parser('prepare')
p1.add_argument('--uniprot-dir', required=True)
p1.add_argument('--extra-fasta', default='', help='comma-separated additional UniProt proteome FASTA.gz (e.g. C. albicans UP000000559)')
p1.add_argument('--out-faa', required=True)
p1.add_argument('--out-meta', required=True)
p2 = sub.add_parser('table')
p2.add_argument('--meta', required=True)
p2.add_argument('--summary', required=True)
p2.add_argument('--out-hits', required=True)
p2.add_argument('--out-summary', required=True)
p2.add_argument('--y1000-table', default='')
a = ap.parse_args()
if a.cmd == 'prepare':
    files = sorted(glob.glob(os.path.join(a.uniprot_dir, 'UP*.fasta.gz')))
    if a.extra_fasta:
        files += [f for f in a.extra_fasta.split(',') if f]
    meta, n = ([], 0)
    with open(a.out_faa, 'w', encoding='utf-8', newline='\n') as out:
        for fn in files:
            m = re.match('(UP\\d+)_(.+?)\\.fasta\\.gz$', os.path.basename(fn))
            upid, org = (m.group(1), m.group(2))
            code = f'{org}_{upid}'
            cnt = 0
            with gzip.open(fn, 'rt', encoding='utf-8') as f:
                hid, buf = (None, [])
                recs = []
                for line in f:
                    if line.startswith('>'):
                        if hid is not None:
                            recs.append((hid, ''.join(buf)))
                        hid, buf = (line[1:].rstrip('\n'), [])
                    else:
                        buf.append(line.strip())
                if hid is not None:
                    recs.append((hid, ''.join(buf)))
            for hid, seq in recs:
                mm = re.match('(\\w+)\\|([A-Za-z0-9_\\-]+)\\|(\\S+)\\s+(.*?)\\s+OS=(.*?)\\s+OX=(\\d+)(?:\\s+GN=(\\S+))?', hid)
                if not mm:
                    raise SystemExit('unparsable UniProt header: ' + hid[:120])
                db, acc, entry, desc, os_, ox, gn = mm.groups()
                out.write(f'>{code}|{acc}\n{seq}\n')
                meta.append(dict(seq_id=f'{code}|{acc}', proteome=upid, organism=org, accession=acc, entry_name=entry, protein_name=desc, gene_name=gn or '', OS=os_, OX=ox, db=db, length=len(seq)))
                cnt += 1
            log(code, cnt)
            n += cnt
    write_tsv(a.out_meta, meta)
    log('outgroup proteins', n, 'from', len(files), 'proteomes')
else:
    meta = {r['seq_id']: r for r in read_tsv(a.meta)}
    summ = read_tsv(a.summary)
    searched = collections.Counter((r['proteome'] + '|' + r['organism'] for r in meta.values()))
    rows = []
    cnt = collections.defaultdict(collections.Counter)
    for r in summ:
        m = meta[r['seq_id']]
        key = m['proteome'] + '|' + m['organism']
        if r['status'] == 'PASS' or r['status'].startswith('NEAR_MISS'):
            cnt[key]['PASS' if r['status'] == 'PASS' else r['status']] += 1
            rows.append(dict(source='UniProt', organism=m['organism'], proteome=m['proteome'], accession=m['accession'], entry_name=m['entry_name'], gene_name=m['gene_name'], protein_name=m['protein_name'], length=r['seq_len'], n_PF00005=r['n_PF00005'], n_PF01061=r['n_PF01061'], n_PF06422=r['n_PF06422'], n_PF14510=r['n_PF14510'], n_PF19055=r['n_PF19055'], coords_PF00005=r['coords_PF00005'], coords_PF01061=r['coords_PF01061'], coords_PF06422=r['coords_PF06422'], coords_PF14510=r['coords_PF14510'], coords_PF19055=r['coords_PF19055'], arch=r['arch'], passes_rule='Y' if r['status'] == 'PASS' else 'N', retrieval_class=r['status'], near_miss_reasons=r['near_miss_reasons'], rule_flags=r['rule_flags']))
        elif r['status'] == 'OTHER_ABC':
            cnt[key]['OTHER_ABC_NBD_only'] += 1
    if a.y1000_table:
        for r in read_tsv(a.y1000_table):
            if r['species_key'] == 'Candida_albicans' and (r['retrieval_class'] == 'PASS' or r['retrieval_class'].startswith('NEAR_MISS')):
                rows.append(dict(source='Y1000+' if r['genome_kind'] == 'y1000_genome' else 'Y1000+ pep archive (CGD reference annotation)', organism='Candida_albicans', proteome=r['species_code'], accession=r['gene_id'], entry_name='', gene_name='', protein_name='', length=r['length_final'], n_PF00005=r['n_PF00005'], n_PF01061=r['n_PF01061'], n_PF06422=r['n_PF06422'], n_PF14510=r['n_PF14510'], n_PF19055=r['n_PF19055'], coords_PF00005=r['coords_PF00005'], coords_PF01061=r['coords_PF01061'], coords_PF06422=r['coords_PF06422'], coords_PF14510=r['coords_PF14510'], coords_PF19055=r['coords_PF19055'], arch=r['arch'], passes_rule=r['passes_rule'], retrieval_class=r['retrieval_class'], near_miss_reasons=r['near_miss_reasons'], rule_flags=r['qc_flags']))
    rows.sort(key=lambda r: (r['organism'], r['source'], r['retrieval_class'] != 'PASS', r['accession']))
    write_tsv(a.out_hits, rows)
    srows = []
    for key, n in sorted(searched.items()):
        c = cnt[key]
        upid, org = key.split('|')
        srows.append(dict(organism=org, proteome=upid, proteins_searched=n, PASS=c['PASS'], **{k: v for k, v in sorted(c.items()) if k != 'PASS'}))
    write_tsv(a.out_summary, srows)
    log('outgroup hits rows', len(rows), 'PASS', sum((1 for r in rows if r['passes_rule'] == 'Y')))

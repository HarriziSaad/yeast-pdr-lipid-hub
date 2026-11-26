import json, csv, sys, os
root = sys.argv[1]
meta = {m['uniprotAccession']: m for m in json.load(open(os.path.join(root, 'alphafold_meta.json')))}
man = {}
for r in csv.DictReader(open(os.path.join(root, 'download_manifest.tsv')), delimiter='\t'):
    man[os.path.basename(r['local_file'])] = r
out = open(sys.argv[2], 'w')
out.write('evidence\tlabel\tid\tsource_version\tmodel_or_release_date\tAFDB_global_pLDDT\tfile\tbytes\tsha256\turl\tretrieved_utc\n')
for acc, m in meta.items():
    for suf in ('model.cif', 'pae.json', 'confidence.json', 'api.json'):
        f = f'{acc}_{suf}'
        r = man[f]
        out.write('\t'.join(['C-AlphaFold', m['label'], m['entryId'], 'AFDB v%s (%s)' % (m['latestVersion'], m['toolUsed']), m['modelCreatedDate'][:10], str(m['globalMetricValue']), 'alphafold/' + f, r['bytes'], r['sha256'], r['url'], r['retrieved_utc']]) + '\n')
for f, r in man.items():
    if r['local_file'].startswith('rcsb_cdr1') or r['local_file'].startswith('pdr5_copy'):
        pid = f.split('.')[0].split('_')[0]
        out.write('\t'.join(['E-cryoEM', r['note'], pid, 'RCSB mmCIF' if f.endswith('.cif') else 'RCSB core entry JSON', 'see pdb_template_table.md', '', r['local_file'].replace('\\', '/'), r['bytes'], r['sha256'], r['url'], r['retrieved_utc']]) + '\n')

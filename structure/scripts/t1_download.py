import sys, os, json, hashlib, shutil, datetime, urllib.request, urllib.error
outroot, resdir = (sys.argv[1], sys.argv[2])
AF = {'P33302': 'Pdr5', 'P32568': 'Snq2', 'P51533': 'Pdr10', 'P40550': 'Pdr11', 'Q02785': 'Pdr12', 'Q04182': 'Pdr15', 'P53756': 'Pdr18', 'Q08409': 'Aus1', 'Q08234': 'YOL075C', 'P25371': 'Adp1', 'Q5ANA3': 'CaCdr1'}
CDR1 = {'9IUK': 'Ca Cdr1 apo (Peng 2024)', '9IUL': 'Ca Cdr1 fluconazole (Peng 2024)', '9IUM': 'Ca Cdr1 milbemycin oxime (Peng 2024)', '22FR': 'Ca Cdr1 apo (Wang 2026)', '9VT1': 'Ca Cdr1 ATP/ADP (Wang 2026)', '9VTM': 'Ca Cdr1 AMPPNP (Wang 2026)', '9VTO': 'Ca Cdr1 AMPPNP/AMPPNP (Wang 2026)', '9VSU': 'Ca Cdr1 FK506 (Wang 2026)', '9VSV': 'Ca Cdr1 curcumin (Wang 2026)', '9VSW': 'Ca Cdr1 beauvericin (Wang 2026)', '9VSZ': 'Ca Cdr1 fluconazole near (Wang 2026)', '9VT2': 'Ca Cdr1 fluconazole middle (Wang 2026)', '9VTN': 'Ca Cdr1 AMPPNP+fluconazole far (Wang 2026)', '9S2C': 'Cg Cdr1 ATP,ADP-VO4 (Pata 2026 inferred)', '9S2H': 'Cg Cdr1 ATP+ADP (Pata 2026 inferred)', '9S4S': 'Cg Cdr1 itraconazole+ATP+ADP-Mg-VO4 (Pata 2026 inferred)', '9S4T': 'Cg Cdr1 apo (Pata 2026 inferred)', '9S4U': 'Cg Cdr1 itraconazole+ATP+ADP-Mg (Pata 2026 inferred)'}
PDR5 = ['7P03', '7P04', '7P05', '7P06']
manifest = []

def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def fetch(url, dest, note=''):
    req = urllib.request.Request(url, headers={'User-Agent': 'pdr-phase3-structure-retrieval/1.0'})
    status, lm = ('', '')
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
            status = str(r.status)
            lm = r.headers.get('Last-Modified', '')
    except urllib.error.HTTPError as e:
        status = f'HTTP{e.code}'
        data = b''
    except Exception as e:
        status = 'ERR:' + type(e).__name__
        data = b''
    ts = now()
    if data:
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        with open(dest, 'wb') as f:
            f.write(data)
        manifest.append((url, ts, len(data), sha256(dest), status, lm, os.path.relpath(dest, outroot), note))
    else:
        manifest.append((url, ts, 0, '', status, lm, '', note + ' FAILED'))
    return data
afdir = os.path.join(outroot, 'alphafold')
af_meta = []
for acc, lab in AF.items():
    api = f'https://alphafold.ebi.ac.uk/api/prediction/{acc}'
    d = fetch(api, os.path.join(afdir, f'{acc}_api.json'), f'AF DB API record {lab}')
    if not d:
        continue
    rec = json.loads(d)
    if len(rec) != 1:
        print('WARNING', acc, 'API returned', len(rec), 'records')
    r = rec[0]
    af_meta.append({k: r.get(k) for k in ('uniprotAccession', 'entryId', 'latestVersion', 'modelCreatedDate', 'sequenceVersionDate', 'globalMetricValue', 'toolUsed', 'sequenceStart', 'sequenceEnd', 'sequenceChecksum', 'cifUrl', 'paeDocUrl', 'plddtDocUrl', 'organismScientificName', 'uniprotId')})
    af_meta[-1]['label'] = lab
    af_meta[-1]['seq_len'] = len(r['sequence'])
    for key, suffix in (('cifUrl', 'model.cif'), ('paeDocUrl', 'pae.json'), ('plddtDocUrl', 'confidence.json')):
        u = r.get(key)
        if u:
            fetch(u, os.path.join(afdir, f'{acc}_{suffix}'), f'AF DB {suffix} {lab}')
        else:
            manifest.append((f'(none in API record: {key})', now(), 0, '', 'ABSENT', '', '', f'{acc} {lab}'))
with open(os.path.join(outroot, 'alphafold_meta.json'), 'w') as f:
    json.dump(af_meta, f, indent=1)
for pid, lab in CDR1.items():
    fetch(f'https://files.rcsb.org/download/{pid}.cif', os.path.join(outroot, 'rcsb_cdr1', f'{pid}.cif'), lab)
    fetch(f'https://data.rcsb.org/rest/v1/core/entry/{pid}', os.path.join(outroot, 'rcsb_cdr1', f'{pid}_entry.json'), 'RCSB core entry record ' + lab)
for pid in PDR5:
    src = os.path.join(resdir, f'{pid}.cif')
    dst = os.path.join(outroot, 'pdr5_copy', f'{pid}.cif')
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(src, dst)
    manifest.append((f'file copy of {src}', now(), os.path.getsize(dst), sha256(dst), 'COPY', '', os.path.relpath(dst, outroot), f'Pdr5 {pid} (phase0 download)'))
with open(os.path.join(outroot, 'download_manifest.tsv'), 'w') as f:
    f.write('url\tretrieved_utc\tbytes\tsha256\thttp_status\tlast_modified\tlocal_file\tnote\n')
    for row in manifest:
        f.write('\t'.join((str(x) for x in row)) + '\n')
print(len(manifest), 'manifest rows;', sum((1 for m in manifest if 'FAILED' in m[7])), 'failed')

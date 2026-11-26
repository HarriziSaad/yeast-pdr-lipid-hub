import csv, gzip, io, json, os, re, sys, time, itertools
SEED = 20261118
BASES = 'TCAG'
CODONS = [''.join(c) for c in itertools.product(BASES, repeat=3)]
CODON_INDEX = {c: i for i, c in enumerate(CODONS)}
NCBI_AA = {1: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', 12: 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', 26: 'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'}
CTG_CALL_TO_TABLE = {'L': 1, 'S': 12, 'A': 26}

def log(*a):
    print(time.strftime('%H:%M:%S'), *a, file=sys.stderr, flush=True)

def open_text(path, mode='rt'):
    if str(path).endswith('.gz'):
        return gzip.open(path, mode, encoding='utf-8', newline='')
    return open(path, mode, encoding='utf-8', newline='')

def read_tsv(path):
    with open_text(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path, rows, fields=None):
    if fields is None:
        fields = []
        for r in rows:
            for k in r:
                if k not in fields:
                    fields.append(k)
    with open(path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

def fasta_iter(path):
    with open_text(path) as f:
        sid, desc, buf = (None, '', [])
        for line in f:
            if line.startswith('>'):
                if sid is not None:
                    yield (sid, desc, ''.join(buf))
                head = line[1:].rstrip('\r\n')
                parts = head.split(None, 1)
                sid = parts[0] if parts else ''
                desc = parts[1] if len(parts) > 1 else ''
                buf = []
            else:
                buf.append(line.strip())
        if sid is not None:
            yield (sid, desc, ''.join(buf))

def write_fasta(fh, sid, seq, desc='', width=0):
    fh.write('>' + sid + (' ' + desc if desc else '') + '\n')
    if width:
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + '\n')
    else:
        fh.write(seq + '\n')
IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
_codon_cache = {}

def _codon_aa(c, table_id, ctg):
    key = (c, table_id, ctg)
    v = _codon_cache.get(key)
    if v is not None:
        return v
    aa = NCBI_AA[table_id]
    try:
        sets = [IUPAC[x] for x in c]
    except KeyError:
        _codon_cache[key] = 'X'
        return 'X'
    res = set()
    for x in sets[0]:
        for y in sets[1]:
            for z in sets[2]:
                cc = x + y + z
                res.add('X' if ctg == 'X' and cc == 'CTG' else aa[CODON_INDEX[cc]])
    v = res.pop() if len(res) == 1 else 'X'
    _codon_cache[key] = v
    return v

def translate(cds, table_id=1, ctg='table'):
    out = []
    for i in range(0, len(cds) - len(cds) % 3, 3):
        out.append(_codon_aa(cds[i:i + 3], table_id, ctg))
    return ''.join(out)

def sanitize(s):
    return re.sub('[^A-Za-z0-9_.\\-]', '_', s)

def species_key_of(final_tip_id):
    t = final_tip_id.split('_')
    if len(t) < 2:
        return (final_tip_id, '')
    if t[1] in ('sp.', 'sp', 'aff.', 'cf.', 'nov.') or t[1].startswith('sp.'):
        return (final_tip_id, '')
    return (t[0] + '_' + t[1], '_'.join(t[2:]))
STOPS = ('TAA', 'TAG', 'TGA')

def ctg_compatible(codon):
    return len(codon) == 3 and all((b in IUPAC and a in IUPAC[b] for a, b in zip('CTG', codon)))

def retranslate_one(cds, pep, ctg_call):
    pep_c = pep[:-1] if pep.endswith('*') else pep
    c = cds
    trailing = len(c) % 3
    body = c[:len(c) - trailing]
    if len(body) >= 3 and body[-3:] in STOPS:
        body = body[:-3]
    std = translate(body, 1)
    n_ctg = sum((1 for i in range(0, len(body), 3) if body[i:i + 3] == 'CTG'))
    note = f'trailing_nt={trailing}' if trailing else ''
    if std != pep_c:
        n_mis = sum((1 for x, y in zip(std, pep_c) if x != y)) + abs(len(std) - len(pep_c))
        return dict(status='MISMATCH_PEP_CDS', seq=pep_c, retranslated='N', table='', n_ctg=n_ctg, n_changed=0, n_mismatch=n_mis, note=(note + ';' if note else '') + f'std_len={len(std)} pep_len={len(pep_c)}')
    if ctg_call == 'L':
        return dict(status='STANDARD_VERIFIED', seq=pep_c, retranslated='N', table=1, n_ctg=n_ctg, n_changed=0, n_mismatch=0, note=note)
    if ctg_call in ('S', 'A'):
        tid = CTG_CALL_TO_TABLE[ctg_call]
        new = translate(body, tid)
        status, tab = ('RETRANSLATED_T%d' % tid, tid)
    elif ctg_call == '?':
        new = translate(body, 1, ctg='X')
        status, tab = ('RETRANSLATED_MASKED_X', '1+CTG->X')
    else:
        return dict(status='UNKNOWN_CALL', seq=pep_c, retranslated='N', table='', n_ctg=n_ctg, n_changed=0, n_mismatch=0, note='ctg_call=' + repr(ctg_call))
    changed = [i for i, (x, y) in enumerate(zip(new, std)) if x != y]
    assert all((ctg_compatible(body[3 * i:3 * i + 3]) for i in changed)), 'a non-CTG codon changed'
    assert len(new) == len(std)
    return dict(status=status, seq=new, retranslated='Y', table=tab, n_ctg=n_ctg, n_changed=len(changed), n_mismatch=0, note=note)

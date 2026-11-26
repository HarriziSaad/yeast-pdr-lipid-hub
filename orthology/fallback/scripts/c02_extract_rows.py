import argparse, hashlib, json, os, sys
ap = argparse.ArgumentParser()
ap.add_argument('--aln', required=True)
ap.add_argument('--ids', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--report', required=True)
ap.add_argument('--expect-cols', type=int, default=None, help='stop if the input alignment does not have this many columns (MAFFT trim50: 1444; MUSCLE trim50: 1477)')
ap.add_argument('--expect-rows', type=int, default=1318)
a = ap.parse_args()
GAPAMB = set('-?*XBZJUOxbzjuo')

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()
seqs, order = ({}, [])
name = None
with open(a.aln, encoding='utf-8') as f:
    for line in f:
        line = line.rstrip('\r\n')
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in seqs:
                sys.exit('duplicate id in alignment: ' + name)
            seqs[name] = []
            order.append(name)
        elif name is not None and line:
            seqs[name].append(line.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}
L = {len(s) for s in seqs.values()}
if len(L) != 1:
    sys.exit('input is not an alignment (row lengths differ): %s' % sorted(L)[:5])
ncol = L.pop()
if len(seqs) != a.expect_rows:
    sys.exit('expected %d rows in the input alignment, found %d' % (a.expect_rows, len(seqs)))
if a.expect_cols is not None and ncol != a.expect_cols:
    sys.exit('expected %d columns, found %d' % (a.expect_cols, ncol))
ids = [l.strip() for l in open(a.ids, encoding='utf-8') if l.strip()]
if len(set(ids)) != len(ids):
    sys.exit('duplicate ids in the id list')
missing = [i for i in ids if i not in seqs]
if missing:
    sys.exit('ids not in the alignment: %s' % missing[:5])
idset = set(ids)
keep_rows = [n for n in order if n in idset]
assert len(keep_rows) == len(ids)
informative = [any((seqs[n][j] not in GAPAMB for n in keep_rows)) for j in range(ncol)]
keep_cols = [j for j, ok in enumerate(informative) if ok]
with open(a.out, 'w', encoding='utf-8', newline='') as f:
    for n in keep_rows:
        f.write('>' + n + '\n' + ''.join((seqs[n][j] for j in keep_cols)) + '\n')
occ = []
for n in keep_rows:
    s = ''.join((seqs[n][j] for j in keep_cols))
    occ.append(sum((c not in GAPAMB for c in s)) / len(s))
rep = dict(input=os.path.basename(a.aln), input_sha256=sha256(a.aln), input_rows=len(seqs), input_cols=ncol, ids_file=os.path.basename(a.ids), ids_sha256=sha256(a.ids), output=os.path.basename(a.out), output_sha256=sha256(a.out), output_rows=len(keep_rows), output_cols=len(keep_cols), columns_removed_all_gap_or_ambiguous=ncol - len(keep_cols), removed_column_indices_0based=[j for j, ok in enumerate(informative) if not ok][:200], min_row_occupancy=round(min(occ), 4), median_row_occupancy=round(sorted(occ)[len(occ) // 2], 4), n_rows_below_0_5=sum((o < 0.5 for o in occ)), pdr5_row_present='Saccharomyces_cerevisiae|g002312.m1' in idset)
json.dump(rep, open(a.report, 'w', encoding='utf-8'), indent=1)
print(json.dumps(rep, indent=1))

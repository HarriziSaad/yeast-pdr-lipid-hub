import argparse, collections, hashlib, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_tsv, write_tsv, log
ap = argparse.ArgumentParser()
ap.add_argument('--aln', required=True)
ap.add_argument('--name', required=True)
ap.add_argument('--outdir', required=True)
ap.add_argument('--manifest', default='/home/pdr/work/phase1b/set/analysis_set_manifest.tsv')
ap.add_argument('--classes', default='/mnt/c/Users/MSI/Documents/Advisor/project/preregistration/site_classes_FROZEN.tsv')
ap.add_argument('--cov', default='/mnt/c/Users/MSI/Documents/Advisor/project/preregistration/outputs/structure_covariates_7P04.tsv')
ap.add_argument('--p33302', default='/mnt/c/Users/MSI/Documents/Advisor/project/phase0_verify/residues/downloads/P33302.fasta')
ap.add_argument('--trim', type=float, default=0.5)
ap.add_argument('--pdr5-id', default='Saccharomyces_cerevisiae|g002312.m1')
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
GAPCH = set('-XBZJUO*?xbzjuo')

def read_fasta(path):
    seqs, order, name, buf = ({}, [], None, [])
    for l in open(path, encoding='utf-8'):
        l = l.rstrip('\n')
        if l.startswith('>'):
            if name is not None:
                seqs[name] = ''.join(buf)
            name = l[1:].split()[0]
            order.append(name)
            buf = []
        else:
            buf.append(l.strip())
    if name is not None:
        seqs[name] = ''.join(buf)
    return (seqs, order)
seqs, order = read_fasta(a.aln)
L = {len(s) for s in seqs.values()}
assert len(L) == 1, 'alignment rows differ in length'
L = L.pop()
N = len(seqs)
man = {r['seq_id']: r for r in read_tsv(a.manifest)}
assert set(man) == set(seqs), (len(set(man) - set(seqs)), len(set(seqs) - set(man)))
ing = [i for i in order if man[i]['role'] == 'INGROUP']
g_all = [0] * L
g_ing = [0] * L
for i in order:
    s = seqs[i]
    isi = man[i]['role'] == 'INGROUP'
    for c, ch in enumerate(s):
        if ch in GAPCH:
            g_all[c] += 1
            if isi:
                g_ing[c] += 1
gf_all = [x / N for x in g_all]
gf_ing = [x / len(ing) for x in g_ing]
keep = [c for c in range(L) if gf_all[c] <= a.trim]
log(a.name, 'N', N, 'L', L, 'kept', len(keep))
with open(os.path.join(a.outdir, f'{a.name}.trim{int(a.trim * 100)}.faa'), 'w') as f:
    for i in order:
        s = seqs[i]
        t = ''.join((s[c] for c in keep))
        f.write('>%s\n' % i)
        for k in range(0, len(t), 80):
            f.write(t[k:k + 80] + '\n')
with open(os.path.join(a.outdir, f'{a.name}.kept_columns.txt'), 'w') as f:
    f.write('\n'.join((str(c + 1) for c in keep)) + '\n')

def occ(i, cols=None):
    s = seqs[i]
    cols = range(L) if cols is None else cols
    n = sum((1 for c in cols if s[c] not in GAPCH))
    return n / len(cols)
occ_all = {i: occ(i) for i in order}
occ_trim = {i: occ(i, keep) for i in order}
rows = []
for i in order:
    nong = sum((1 for ch in seqs[i] if ch not in GAPCH))
    rows.append(dict(seq_id=i, role=man[i]['role'], order=man[i]['order'], residues=nong, occ_full=round(occ_all[i], 4), occ_trimmed=round(occ_trim[i], 4)))
write_tsv(os.path.join(a.outdir, f'{a.name}.seq_occupancy.tsv'), rows)

def q(v, p):
    v = sorted(v)
    return v[int(round(p * (len(v) - 1)))]
p5 = seqs[a.pdr5_id]
p33302 = ''.join((l.strip() for l in open(a.p33302) if not l.startswith('>')))
p5_ungapped = p5.replace('-', '')
identical = p5_ungapped == p33302
col_of = []
for c, ch in enumerate(p5):
    if ch != '-':
        col_of.append(c)
mapping_note = 'identical to P33302' if identical else 'DIFFERS from P33302'
if not identical:
    from Bio import Align
    al = Align.PairwiseAligner()
    al.mode = 'global'
    al.match_score, al.mismatch_score, al.open_gap_score, al.extend_gap_score = (2, -1, -5, -0.5)
    best = al.align(p33302, p5_ungapped)[0]
    mp = {}
    for (s1, e1), (s2, e2) in zip(*best.aligned):
        for k in range(e1 - s1):
            mp[s1 + k + 1] = s2 + k + 1
    ident = sum((1 for k, v in mp.items() if p33302[k - 1] == p5_ungapped[v - 1]))
    mapping_note = f'DIFFERS: {len(p33302)} vs {len(p5_ungapped)} aa; pairwise-mapped {len(mp)}, identical at {ident}'
else:
    mp = {k: k for k in range(1, len(p33302) + 1)}
cls = read_tsv(a.classes)
by_res = collections.defaultdict(lambda: dict(aa='', classes=[]))
for r in cls:
    d = by_res[int(r['residue'])]
    d['aa'] = r['aa']
    d['classes'].append(r['class'])
keepset = set(keep)
resrows = []
for resn in sorted(by_res):
    d = by_res[resn]
    k = mp.get(resn)
    if k is None:
        resrows.append(dict(residue=resn, aa=d['aa'], classes=';'.join(sorted(set(d['classes']))), column='', aa_in_pdr5_row='', aa_match='', gapfrac_all='', gapfrac_ingroup='', retained_in_trimmed='N', NA_for_scoring_ingroup_gt50='', note='no counterpart'))
        continue
    c = col_of[k - 1]
    resrows.append(dict(residue=resn, aa=d['aa'], classes=';'.join(sorted(set(d['classes']))), column=c + 1, aa_in_pdr5_row=p5[c], aa_match='Y' if p5[c] == d['aa'] else 'N', gapfrac_all=round(gf_all[c], 4), gapfrac_ingroup=round(gf_ing[c], 4), retained_in_trimmed='Y' if c in keepset else 'N', NA_for_scoring_ingroup_gt50='Y' if gf_ing[c] > 0.5 else 'N', note=''))
write_tsv(os.path.join(a.outdir, f'{a.name}.class_residue_columns.tsv'), resrows)
rmap = {int(r['residue']): r for r in resrows}
crow = []
classes_all = sorted({c for d in by_res.values() for c in d['classes']})
for cname in classes_all:
    rs = sorted((r for r, d in by_res.items() if cname in d['classes']))
    ret = sum((1 for r in rs if rmap[r]['retained_in_trimmed'] == 'Y'))
    na = sum((1 for r in rs if rmap[r]['NA_for_scoring_ingroup_gt50'] == 'Y'))
    mism = sum((1 for r in rs if rmap[r]['aa_match'] == 'N'))
    crow.append(dict(cls=cname, n_residues=len(rs), retained_trimmed=ret, lost_trimmed=len(rs) - ret, NA_scoring_ingroup=na, aa_mismatch_vs_P33302=mism, min_gapfrac_ingroup=min((float(rmap[r]['gapfrac_ingroup']) for r in rs)), max_gapfrac_ingroup=max((float(rmap[r]['gapfrac_ingroup']) for r in rs))))
n_ret_all = sum((1 for r in rmap.values() if r['retained_in_trimmed'] == 'Y'))
crow.append(dict(cls='ALL_273_DISTINCT', n_residues=len(rmap), retained_trimmed=n_ret_all, lost_trimmed=len(rmap) - n_ret_all, NA_scoring_ingroup=sum((1 for r in rmap.values() if r['NA_for_scoring_ingroup_gt50'] == 'Y')), aa_mismatch_vs_P33302=sum((1 for r in rmap.values() if r['aa_match'] == 'N')), min_gapfrac_ingroup='', max_gapfrac_ingroup=''))
write_tsv(os.path.join(a.outdir, f'{a.name}.class_survival.tsv'), crow)
cov = read_tsv(a.cov)
tm = [int(r['residue']) for r in cov if r['TM_helix'] == '1']
runs, st, pv = ([], None, None)
for x in tm:
    if st is None:
        st = pv = x
    elif x == pv + 1:
        pv = x
    else:
        runs.append((st, pv))
        st = pv = x
runs.append((st, pv))
hrows = []
for hi, (s, e) in enumerate(runs, 1):
    cols = [col_of[mp[x] - 1] for x in range(s, e + 1) if x in mp]
    ret = sum((1 for c in cols if c in keepset))
    hrows.append(dict(segment=hi, start=s, end=e, n_res=e - s + 1, retained_trimmed=ret, frac=round(ret / max(1, len(cols)), 3), median_gapfrac_ingroup=round(sorted((gf_ing[c] for c in cols))[len(cols) // 2], 4), monotone_columns='Y' if cols == sorted(cols) else 'N'))
write_tsv(os.path.join(a.outdir, f'{a.name}.tm_segments.tsv'), hrows)
summ = dict(name=a.name, aln=a.aln, aln_sha256=hashlib.sha256(open(a.aln, 'rb').read()).hexdigest(), n_seq=N, n_ingroup=len(ing), alignment_length=L, trim_gap_threshold=a.trim, retained_columns=len(keep), retained_fraction=round(len(keep) / L, 4), removed_columns=L - len(keep), gapfrac_all_mean_full=round(sum(gf_all) / L, 4), gapfrac_all_mean_trimmed=round(sum((gf_all[c] for c in keep)) / len(keep), 4), gapfrac_ingroup_mean_full=round(sum(gf_ing) / L, 4), gapfrac_ingroup_mean_trimmed=round(sum((gf_ing[c] for c in keep)) / len(keep), 4), columns_gapfrac_all_gt_0_9=sum((1 for x in gf_all if x > 0.9)), columns_gapfrac_all_gt_0_5=sum((1 for x in gf_all if x > 0.5)), columns_gapfrac_all_le_0_05=sum((1 for x in gf_all if x <= 0.05)), columns_ingroup_gapfrac_gt_0_5=sum((1 for x in gf_ing if x > 0.5)), seq_occupancy_trimmed_min=round(min(occ_trim.values()), 4), seq_occupancy_trimmed_median=round(q(list(occ_trim.values()), 0.5), 4), n_seq_trimmed_occupancy_lt_0_5=sum((1 for v in occ_trim.values() if v < 0.5)), n_seq_trimmed_occupancy_lt_0_3=sum((1 for v in occ_trim.values() if v < 0.3)), pdr5_row=a.pdr5_id, pdr5_vs_P33302=mapping_note, n_class_residues=len(rmap), class_residues_retained=n_ret_all, tm_segments=len(runs), tm_residues_retained=sum((h['retained_trimmed'] for h in hrows)), tm_residues=sum((h['n_res'] for h in hrows)), tm_segments_monotone=sum((1 for h in hrows if h['monotone_columns'] == 'Y')))
json.dump(summ, open(os.path.join(a.outdir, f'{a.name}.aln_qc.json'), 'w'), indent=1)
print(json.dumps(summ, indent=1))
for r in crow:
    print(r)

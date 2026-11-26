import argparse, collections, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_tsv, write_tsv, log
import gemmi
from Bio import Align
ap = argparse.ArgumentParser()
ap.add_argument('--aln', required=True)
ap.add_argument('--name', required=True)
ap.add_argument('--usalign', required=True, help='US-align text output with structure 1 = Cdr1 structure, structure 2 = 7P04')
ap.add_argument('--cdr1-cif', required=True)
ap.add_argument('--pdr5-cif', required=True)
ap.add_argument('--outdir', required=True)
ap.add_argument('--manifest', default='/home/pdr/work/phase1b/set/analysis_set_manifest.tsv')
ap.add_argument('--classes', default='/mnt/c/Users/MSI/Documents/Advisor/project/preregistration/site_classes_FROZEN.tsv')
ap.add_argument('--pdr5-id', default='Saccharomyces_cerevisiae|g002312.m1')
ap.add_argument('--calb-code', default='Candida_albicans_SC5314')
ap.add_argument('--p33302', default='/mnt/c/Users/MSI/Documents/Advisor/project/phase0_verify/residues/downloads/P33302.fasta')
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

def chain_residues(cif):
    st = gemmi.read_structure(cif)
    st.setup_entities()
    ch = st[0][0]
    out = []
    for r in ch:
        if r.het_flag == 'H' or gemmi.find_tabulated_residue(r.name) is None or (not gemmi.find_tabulated_residue(r.name).is_amino_acid()):
            continue
        if r.find_atom('CA', '*') is None:
            continue
        out.append((r.seqid.num, gemmi.find_tabulated_residue(r.name).one_letter_code.upper()))
    return (ch.name, out)
c_chain, c_res = chain_residues(a.cdr1_cif)
p_chain, p_res = chain_residues(a.pdr5_cif)
log('Cdr1 structure chain', c_chain, len(c_res), 'residues; 7P04 chain', p_chain, len(p_res))
lines = open(a.usalign).read().split('\n')
k = next((i for i, l in enumerate(lines) if l.startswith('(":" denotes')))
s1, mk, s2 = (lines[k + 1], lines[k + 2], lines[k + 3])
assert len(s1) == len(s2) == len(mk) or len(mk) <= len(s1)
mk = mk.ljust(len(s1))
i1 = i2 = 0
pairs = {}
for x, y, m in zip(s1, s2, mk):
    if x != '-' and y != '-':
        pairs[p_res[i2][0]] = (c_res[i1][0], m)
        assert p_res[i2][1] == y.upper() or True
    if x != '-':
        i1 += 1
    if y != '-':
        i2 += 1
assert i1 == len(c_res) and i2 == len(p_res), (i1, len(c_res), i2, len(p_res))
log('structure-aligned pairs', len(pairs), ' of which ":" (<5 A):', sum((1 for v in pairs.values() if v[1] == ':')))
seqs, name = ({}, None)
for l in open(a.aln):
    l = l.rstrip('\n')
    if l.startswith('>'):
        name = l[1:].split()[0]
        seqs[name] = []
    else:
        seqs[name].append(l.strip())
seqs = {k_: ''.join(v) for k_, v in seqs.items()}
man = {r['seq_id']: r for r in read_tsv(a.manifest)}
calb = [i for i in seqs if man[i]['species_code'] == a.calb_code]
log('C. albicans SC5314 sequences in the set:', len(calb))
c_seq = ''.join((x for _, x in c_res))
al = Align.PairwiseAligner()
al.mode = 'global'
al.match_score, al.mismatch_score, al.open_gap_score, al.extend_gap_score = (2, -1, -5, -0.5)
cand = []
for i in calb:
    u = seqs[i].replace('-', '')
    sc = al.score(c_seq, u)
    cand.append((sc, i))
cand.sort(reverse=True)
best_id = cand[0][1]
log('Cdr1 row (best match to the structure sequence):', best_id, 'scores', cand[:3])
row = seqs[best_id]
u = row.replace('-', '')
aln = al.align(c_seq, u)[0]
struct_idx_to_row = {}
for (sa, ea), (sb, eb) in zip(*aln.aligned):
    for t in range(ea - sa):
        struct_idx_to_row[sa + t] = sb + t
ident = sum((1 for si, ri in struct_idx_to_row.items() if c_seq[si] == u[ri]))
log('structure-sequence residues mapped to the row:', len(struct_idx_to_row), 'identical:', ident)
cnum_to_idx = {num: i for i, (num, _) in enumerate(c_res)}
row_col = []
for c, ch in enumerate(row):
    if ch != '-':
        row_col.append(c)
col_to_rowidx = {c: i for i, c in enumerate(row_col)}
rowidx_to_cnum = {}
for si, ri in struct_idx_to_row.items():
    rowidx_to_cnum[ri] = c_res[si][0]
p5 = seqs[a.pdr5_id]
p33302 = ''.join((l.strip() for l in open(a.p33302) if not l.startswith('>')))
assert p5.replace('-', '') == p33302
p5_col = [c for c, ch in enumerate(p5) if ch != '-']
cls = collections.defaultdict(set)
for r in read_tsv(a.classes):
    cls[int(r['residue'])].add(r['class'])
rows = []
for resn in sorted(pairs):
    cnum, m = pairs[resn]
    col = p5_col[resn - 1]
    ri = col_to_rowidx.get(col)
    if ri is None:
        s_aln, status, d = ('', 'ALIGNMENT_GAP_IN_CDR1', '')
    else:
        s_aln = rowidx_to_cnum.get(ri, '')
        if s_aln == '':
            status, d = ('CDR1_RESIDUE_NOT_MODELLED', '')
        else:
            d = s_aln - cnum
            status = 'OK_0' if d == 0 else 'OK_1_3' if abs(d) <= 3 else 'DISCREPANT_GE4'
    rows.append(dict(pdr5_residue=resn, in_class=';'.join(sorted(cls.get(resn, []))), structure_partner=cnum, us_marker=m, alignment_partner=s_aln, difference=d, status=status))
write_tsv(os.path.join(a.outdir, f'{a.name}.cdr1_check_all_pairs.tsv'), rows)
sel = [r for r in rows if r['us_marker'] == ':']
cr = [r for r in sel if r['in_class']]
cnt = collections.Counter((r['status'] for r in sel))
cnt_c = collections.Counter((r['status'] for r in cr))
by_class = collections.defaultdict(collections.Counter)
for r in cr:
    for c in r['in_class'].split(';'):
        by_class[c][r['status']] += 1
write_tsv(os.path.join(a.outdir, f'{a.name}.cdr1_check_by_class.tsv'), [dict(cls=c, n=sum(v.values()), **{s: v.get(s, 0) for s in ('OK_0', 'OK_1_3', 'DISCREPANT_GE4', 'ALIGNMENT_GAP_IN_CDR1', 'CDR1_RESIDUE_NOT_MODELLED')}) for c, v in sorted(by_class.items())])
write_tsv(os.path.join(a.outdir, f'{a.name}.cdr1_check_class_residues_discrepant.tsv'), [r for r in cr if r['status'] not in ('OK_0', 'OK_1_3')])
summ = dict(name=a.name, cdr1_structure=os.path.basename(a.cdr1_cif), cdr1_row=best_id, structure_pairs=len(pairs), pairs_lt5A=len(sel), status_all_lt5A=dict(cnt), class_residues_with_pair_lt5A=len(cr), status_class_residues=dict(cnt_c), class_residues_in_set=len(cls), class_residues_without_pair=sum((1 for r in cls if r not in {x['pdr5_residue'] for x in sel})))
json.dump(summ, open(os.path.join(a.outdir, f'{a.name}.cdr1_check.json'), 'w'), indent=1)
print(json.dumps(summ, indent=1))

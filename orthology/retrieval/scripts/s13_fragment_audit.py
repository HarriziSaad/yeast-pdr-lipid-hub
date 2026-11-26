import argparse, collections, json, os, subprocess, tempfile
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import pandas as pd
from scipy.stats import chi2_contingency
from common import read_tsv, write_fasta, fasta_iter, log
ap = argparse.ArgumentParser()
ap.add_argument('--retrieval', required=True)
ap.add_argument('--members', required=True)
ap.add_argument('--final-faa', required=True)
ap.add_argument('--seeds', required=True)
ap.add_argument('--outdir', required=True)
ap.add_argument('--threads', type=int, default=8)
ap.add_argument('--low-margin-frac', type=float, default=0.05, help='LOW_MARGIN if (best - second) < this fraction of the best bit score (screening threshold)')
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
FAMILY = {'PDR5': 'PDR5-like', 'PDR15': 'PDR5-like', 'PDR10': 'PDR5-like', 'SNQ2': 'SNQ2-like', 'PDR18': 'SNQ2-like', 'PDR12': 'PDR12-like', 'PDR11': 'PDR11/AUS1-like', 'AUS1': 'PDR11/AUS1-like', 'no_hit': 'no_hit'}
FAMREF = {'PDR5-like': 'P33302', 'SNQ2-like': 'P32568', 'PDR12-like': 'Q02785', 'PDR11/AUS1-like': 'P40550'}
NAME = {'P33302': 'PDR5', 'P32568': 'SNQ2', 'P51533': 'PDR10', 'P40550': 'PDR11', 'Q02785': 'PDR12', 'Q04182': 'PDR15', 'P53756': 'PDR18', 'Q08409': 'AUS1'}
ret = {r['seq_id']: r for r in read_tsv(a.retrieval)}
mem = [r['seq_id'] for r in read_tsv(a.members)]
mem = [m for m in mem if ret.get(m, {}).get('genome_kind') == 'y1000_genome']
log('OG0000012 members in Y1000+ genomes', len(mem))
seq = {}
want = set(mem)
for sid, _, s in fasta_iter(a.final_faa):
    if sid in want:
        seq[sid] = s
tmp = tempfile.mkdtemp(prefix='fragaudit_')
seedlen = {}
with open(os.path.join(tmp, 'seeds.faa'), 'w') as fh:
    for sid, d, s in fasta_iter(a.seeds):
        acc = sid.split('|')[1]
        write_fasta(fh, acc, s)
        seedlen[acc] = len(s)
with open(os.path.join(tmp, 'q.faa'), 'w') as fh:
    for m in mem:
        write_fasta(fh, m, seq[m].replace('*', ''))
subprocess.run(['diamond', 'makedb', '--in', os.path.join(tmp, 'seeds.faa'), '-d', os.path.join(tmp, 'seeds')], check=True, capture_output=True)
p = subprocess.run(['diamond', 'blastp', '-d', os.path.join(tmp, 'seeds'), '-q', os.path.join(tmp, 'q.faa'), '--very-sensitive', '-k', '8', '--evalue', '1e-3', '-p', str(a.threads), '--outfmt', '6', 'qseqid', 'sseqid', 'bitscore', 'evalue', 'qlen', 'qstart', 'qend', 'sstart', 'send'], check=True, capture_output=True, text=True)
best = collections.defaultdict(dict)
for line in p.stdout.splitlines():
    q, s, b, e, ql, qs, qe, ss, se = line.split('\t')
    b = float(b)
    if s not in best[q] or b > best[q][s][0]:
        best[q][s] = (b, int(qs), int(qe), min(int(ss), int(se)), max(int(ss), int(se)))
rows = []
for m in mem:
    r = ret[m]
    h = sorted(best.get(m, {}).items(), key=lambda kv: -kv[1][0])
    cls = r['retrieval_class']
    if cls == 'PASS':
        reason = 'PASS'
    elif cls.startswith('NEAR_MISS'):
        reason = r['near_miss_class']
    else:
        reason = cls.split(':', 1)[1]
    d = dict(seq_id=m, species_code=r['species_code'], order=r['order'], CUG_code=r['CUG_code'], length=r['length_final'], retrieval_class=cls, reason=reason, pass_flag='PASS' if cls == 'PASS' else 'nonPASS')
    if h:
        (s1, v1), v2 = (h[0], h[1] if len(h) > 1 else None)
        marg = v1[0] - (v2[1][0] if v2 else 0.0)
        d.update(best_paralog=NAME[s1], best_bits=round(v1[0], 1), second_paralog=NAME[v2[0]] if v2 else '', second_bits=round(v2[1][0], 1) if v2 else '', margin_bits=round(marg, 1), low_margin='Y' if marg < a.low_margin_frac * v1[0] else 'N', q_start=v1[1], q_end=v1[2], seed_start=v1[3], seed_end=v1[4], seed_cov=round((v1[4] - v1[3] + 1) / seedlen[s1], 3), seed_mid_frac=round((v1[3] + v1[4]) / 2 / seedlen[s1], 3))
    else:
        d.update(best_paralog='no_hit', best_bits='', second_paralog='', second_bits='', margin_bits='', low_margin='', q_start='', q_end='', seed_start='', seed_end='', seed_cov='', seed_mid_frac='')
    rows.append(d)
df = pd.DataFrame(rows)
df['family'] = df.best_paralog.map(FAMILY)
refseq = {sid.split('|')[1]: s for sid, d, s in fasta_iter(a.seeds)}
with open(os.path.join(tmp, 'refs.faa'), 'w') as fh:
    for acc in FAMREF.values():
        write_fasta(fh, acc, refseq[acc])
subprocess.run(['diamond', 'makedb', '--in', os.path.join(tmp, 'refs.faa'), '-d', os.path.join(tmp, 'refs')], check=True, capture_output=True)
p2 = subprocess.run(['diamond', 'blastp', '-d', os.path.join(tmp, 'refs'), '-q', os.path.join(tmp, 'q.faa'), '--very-sensitive', '-k', '4', '--evalue', '1e-3', '-p', str(a.threads), '--outfmt', '6', 'qseqid', 'sseqid', 'bitscore', 'sstart', 'send'], check=True, capture_output=True, text=True)
refspan = {}
for line in p2.stdout.splitlines():
    q, s, b, ss, se = line.split('\t')
    b = float(b)
    if (q, s) not in refspan or b > refspan[q, s][0]:
        refspan[q, s] = (b, min(int(ss), int(se)), max(int(ss), int(se)))
df['famref_start'] = [refspan.get((m, FAMREF.get(f, '')), (0, -1, -1))[1] for m, f in zip(df.seq_id, df.family)]
df['famref_end'] = [refspan.get((m, FAMREF.get(f, '')), (0, -1, -1))[2] for m, f in zip(df.seq_id, df.family)]
df.to_csv(os.path.join(a.outdir, 'og_member_paralog_assignment.tsv'), sep='\t', index=False)
reasons = ['PASS', 'HALF_SIZE_ABCG', 'PARTIAL_DOMAINS', 'LEN_LONG', 'LEN_SHORT', 'NO_PF06422', 'OTHER_ABC', 'NO_DOMAIN_HIT']

def summarize(g, key):
    t = g.groupby(key).apply(lambda x: pd.Series({**{'n_members': len(x), 'n_PASS': int((x.pass_flag == 'PASS').sum()), 'n_nonPASS': int((x.pass_flag != 'PASS').sum())}, **{'n_' + r: int((x.reason == r).sum()) for r in reasons[1:]}, 'n_low_margin': int((x.low_margin == 'Y').sum())}), include_groups=False)
    t['frac_nonPASS'] = (t.n_nonPASS / t.n_members).round(3)
    for r in reasons[1:]:
        t['frac_nonPASS_' + r] = (t['n_' + r] / t.n_nonPASS.where(t.n_nonPASS > 0)).round(3)
    return t.reset_index()
by_par = summarize(df, 'best_paralog')
by_par.to_csv(os.path.join(a.outdir, 'audit_by_paralog.tsv'), sep='\t', index=False)
by_ord = summarize(df, 'order')
by_ord.to_csv(os.path.join(a.outdir, 'audit_by_order.tsv'), sep='\t', index=False)
by_ctg = summarize(df, 'CUG_code')
by_ctg.to_csv(os.path.join(a.outdir, 'audit_by_ctg_code.tsv'), sep='\t', index=False)
xo = df.groupby(['order', 'best_paralog']).apply(lambda x: pd.Series({'n_members': len(x), 'n_nonPASS': int((x.pass_flag != 'PASS').sum()), 'frac_nonPASS': round((x.pass_flag != 'PASS').mean(), 3)}), include_groups=False).reset_index()
xo.to_csv(os.path.join(a.outdir, 'audit_order_x_paralog.tsv'), sep='\t', index=False)
nz = df[(df.pass_flag != 'PASS') & (df.best_paralog != 'no_hit')].copy()
nz['seed_cov'] = nz.seed_cov.astype(float)
nz['seed_mid_frac'] = nz.seed_mid_frac.astype(float)
nz['part'] = pd.cut(nz.seed_mid_frac, [0, 1 / 3, 2 / 3, 1.0001], labels=['N-terminal third', 'middle third', 'C-terminal third'])
cov = nz.groupby(['best_paralog', 'part'], observed=False).size().unstack(fill_value=0)
cov['median_seed_cov'] = nz.groupby('best_paralog').seed_cov.median().round(3)
cov.reset_index().to_csv(os.path.join(a.outdir, 'audit_seed_coverage.tsv'), sep='\t', index=False)
pairs = []
for (sp, par), g in df[df.best_paralog != 'no_hit'].groupby(['species_code', 'best_paralog']):
    has_pass = bool((g.pass_flag == 'PASS').any())
    iv = sorted(((int(x.seed_start), int(x.seed_end)) for x in g[g.pass_flag != 'PASS'].itertuples()))
    cov_len, cur_e = (0, 0)
    for s, e in iv:
        s = max(s, cur_e + 1)
        if e >= s:
            cov_len += e - s + 1
            cur_e = e
    ucov = cov_len / [v for k, v in seedlen.items() if NAME[k] == par][0] if iv else 0.0
    pairs.append(dict(species_code=sp, order=g['order'].iloc[0], paralog=par, n_members=len(g), n_PASS=int((g.pass_flag == 'PASS').sum()), n_nonPASS=int((g.pass_flag != 'PASS').sum()), has_PASS='Y' if has_pass else 'N', union_seed_coverage_of_nonPASS=round(ucov, 3)))
pr = pd.DataFrame(pairs)
pr.to_csv(os.path.join(a.outdir, 'audit_pairs_genome_paralog.tsv'), sep='\t', index=False)
nopass = pr[pr.has_PASS == 'N']

def pair_sum(g, key):
    t = g.groupby(key).apply(lambda x: pd.Series({'n_pairs_genome_x_paralog': len(x), 'n_pairs_with_PASS': int((x.has_PASS == 'Y').sum()), 'n_pairs_only_fragments': int((x.has_PASS == 'N').sum()), 'frac_pairs_only_fragments': round((x.has_PASS == 'N').mean(), 3), 'only_fragments_union_cov_ge_0.8': int(((x.has_PASS == 'N') & (x.union_seed_coverage_of_nonPASS >= 0.8)).sum()), 'only_fragments_union_cov_lt_0.3': int(((x.has_PASS == 'N') & (x.union_seed_coverage_of_nonPASS < 0.3)).sum())}), include_groups=False)
    return t.reset_index()
ps_par = pair_sum(pr, 'paralog')
ps_ord = pair_sum(pr, 'order')
ps_par.insert(0, 'by', 'paralog')
ps_par = ps_par.rename(columns={'paralog': 'group'})
ps_ord.insert(0, 'by', 'order')
ps_ord = ps_ord.rename(columns={'order': 'group'})
pd.concat([ps_par, ps_ord]).to_csv(os.path.join(a.outdir, 'audit_pairs_summary.tsv'), sep='\t', index=False)
fam = summarize(df, 'family')
fam.to_csv(os.path.join(a.outdir, 'audit_by_family.tsv'), sep='\t', index=False)
xf = df.groupby(['order', 'family']).apply(lambda x: pd.Series({'n_members': len(x), 'n_nonPASS': int((x.pass_flag != 'PASS').sum()), 'frac_nonPASS': round((x.pass_flag != 'PASS').mean(), 3)}), include_groups=False).reset_index()
xf.to_csv(os.path.join(a.outdir, 'audit_order_x_family.tsv'), sep='\t', index=False)
reflen = {f: seedlen[acc] for f, acc in FAMREF.items()}
fpairs = []
for (sp, f), g in df[df.family != 'no_hit'].groupby(['species_code', 'family']):
    has_pass = bool((g.pass_flag == 'PASS').any())
    iv = sorted(((int(x.famref_start), int(x.famref_end)) for x in g[g.pass_flag != 'PASS'].itertuples() if int(x.famref_start) > 0))
    cov_len, cur_e = (0, 0)
    for s_, e_ in iv:
        s_ = max(s_, cur_e + 1)
        if e_ >= s_:
            cov_len += e_ - s_ + 1
            cur_e = e_
    fpairs.append(dict(species_code=sp, order=g['order'].iloc[0], family=f, n_members=len(g), n_PASS=int((g.pass_flag == 'PASS').sum()), n_nonPASS=int((g.pass_flag != 'PASS').sum()), has_PASS='Y' if has_pass else 'N', union_ref_coverage_of_nonPASS=round(cov_len / reflen[f], 3) if iv else 0.0))
fp = pd.DataFrame(fpairs)
fp.to_csv(os.path.join(a.outdir, 'audit_pairs_genome_family.tsv'), sep='\t', index=False)

def fpair_sum(g, key):
    t_ = g.groupby(key).apply(lambda x: pd.Series({'n_pairs_genome_x_family': len(x), 'n_pairs_with_PASS': int((x.has_PASS == 'Y').sum()), 'n_pairs_only_fragments': int((x.has_PASS == 'N').sum()), 'frac_pairs_only_fragments': round((x.has_PASS == 'N').mean(), 3), 'only_fragments_union_cov_ge_0.8': int(((x.has_PASS == 'N') & (x.union_ref_coverage_of_nonPASS >= 0.8)).sum()), 'only_fragments_union_cov_lt_0.3': int(((x.has_PASS == 'N') & (x.union_ref_coverage_of_nonPASS < 0.3)).sum())}), include_groups=False)
    return t_.reset_index()
f1 = fpair_sum(fp, 'family')
f1.insert(0, 'by', 'family')
f1 = f1.rename(columns={'family': 'group'})
f2 = fpair_sum(fp, 'order')
f2.insert(0, 'by', 'order')
f2 = f2.rename(columns={'order': 'group'})
pd.concat([f1, f2]).to_csv(os.path.join(a.outdir, 'audit_pairs_family_summary.tsv'), sep='\t', index=False)
summ = dict(n_members=len(df), n_PASS=int((df.pass_flag == 'PASS').sum()), n_nonPASS=int((df.pass_flag != 'PASS').sum()), n_no_hit=int((df.best_paralog == 'no_hit').sum()), no_hit_reasons=df[df.best_paralog == 'no_hit'].reason.value_counts().to_dict(), n_low_margin=int((df.low_margin == 'Y').sum()))
for name, key in (('paralog', 'best_paralog'), ('family', 'family'), ('order', 'order'), ('ctg_code', 'CUG_code')):
    sub = df[df.best_paralog != 'no_hit'] if key in ('best_paralog', 'family') else df
    ct = pd.crosstab(sub[key], sub.pass_flag)
    chi, pv, dof, _ = chi2_contingency(ct)
    n = ct.values.sum()
    summ['heterogeneity_' + name] = dict(chi2=round(float(chi), 1), dof=int(dof), p=float(pv), cramers_v=round(float((chi / (n * (min(ct.shape) - 1))) ** 0.5), 3))
summ['non_pass_fraction_overall'] = round(summ['n_nonPASS'] / summ['n_members'], 3)
summ['repair_candidate_members'] = int(sum((1 for m in mem if ret[m].get('repair_candidate') == 'Y')))
summ['family_pairs_total'] = len(fp)
summ['family_pairs_only_fragments'] = int((fp.has_PASS == 'N').sum())
summ['family_pairs_only_fragments_cov_ge_0.8'] = int(((fp.has_PASS == 'N') & (fp.union_ref_coverage_of_nonPASS >= 0.8)).sum())
summ['pairs_total'] = len(pr)
summ['pairs_only_fragments'] = int((pr.has_PASS == 'N').sum())
summ['pairs_only_fragments_cov_ge_0.8'] = int(((pr.has_PASS == 'N') & (pr.union_seed_coverage_of_nonPASS >= 0.8)).sum())
json.dump(summ, open(os.path.join(a.outdir, 'audit_summary.json'), 'w'), indent=1)
log(json.dumps(summ))
print(by_par.to_string())
print(fam.to_string())
print(pd.concat([f1, f2]).to_string())
print(by_ord[['order', 'n_members', 'n_PASS', 'n_nonPASS', 'frac_nonPASS']].to_string())

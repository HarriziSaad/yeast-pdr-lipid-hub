import argparse, json, os
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import pandas as pd
from common import log
ap = argparse.ArgumentParser()
ap.add_argument('--audit-dir', required=True)
ap.add_argument('--plan', required=True)
ap.add_argument('--out', required=True)
a = ap.parse_args()
m = pd.read_csv(os.path.join(a.audit_dir, 'og_member_paralog_assignment.tsv'), sep='\t')
pairs = pd.read_csv(os.path.join(a.audit_dir, 'audit_pairs_genome_family.tsv'), sep='\t')
plan = pd.read_csv(a.plan, sep='\t').fillna('')
lay = plan.set_index('speciesCode')[['selected', 'layer', 'representative_of_species']]
m = m.join(lay, on='species_code')
pairs = pairs.join(lay, on='species_code')

def tag(layer, sel):
    if sel != 'Y':
        return 'not selected'
    ls = set(layer.split('+'))
    if 'WGD' in ls:
        return 'selected: WGD block'
    if ls & {'DIAS', 'YGOB', 'ANCHOR'}:
        return 'selected: Dias/YGOB/anchor (non-WGD)'
    return 'selected: stratified'
m['group'] = [tag(l, s) for l, s in zip(m.layer, m.selected)]
pairs['group'] = [tag(l, s) for l, s in zip(pairs.layer, pairs.selected)]
rows = []
for g in ['selected: WGD block', 'selected: Dias/YGOB/anchor (non-WGD)', 'selected: stratified', 'not selected']:
    x, y = (m[m.group == g], pairs[pairs.group == g])
    rows.append(dict(group=g, genomes=x.species_code.nunique(), OG_members=len(x), PASS=int((x.pass_flag == 'PASS').sum()), nonPASS=int((x.pass_flag != 'PASS').sum()), frac_nonPASS=round((x.pass_flag != 'PASS').mean(), 3), pairs_genome_x_family=len(y), pairs_only_fragments=int((y.has_PASS == 'N').sum()), frac_pairs_only_fragments=round((y.has_PASS == 'N').mean(), 3), only_fragments_cov_ge_0_8=int(((y.has_PASS == 'N') & (y.union_ref_coverage_of_nonPASS >= 0.8)).sum()), pairs_with_PASS_and_fragments_covering_ge_0_8=int(((y.has_PASS == 'Y') & (y.union_ref_coverage_of_nonPASS >= 0.8)).sum())))
t = pd.DataFrame(rows)
t.to_csv(a.out, sep='\t', index=False)
print(t.to_string())
allp = pairs
print('all pairs: with PASS and non-PASS fragments covering >=0.8 of the family reference (probable additional, fragmented copy):', int(((allp.has_PASS == 'Y') & (allp.union_ref_coverage_of_nonPASS >= 0.8)).sum()), 'of', int((allp.has_PASS == 'Y').sum()))
w = m[m.group == 'selected: WGD block']
print('WGD block members per family:\n', pd.crosstab(w.family, w.pass_flag).to_string())

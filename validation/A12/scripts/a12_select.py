import csv
import re
from collections import Counter, defaultdict
UNBIASED_CORE = ('egner1998', 'tanabe2019', 'niimi2012', 'niimi2022', 'shekhar2016', 'zhao2026', 'sauna2008', 'downes2013', 'ananthaswamy2010')
UNBIASED_ABSTRACT = ('kolaczkowski2013',)
DESIGNED_RE = re.compile('site-directed|site directed|\\bdesigned\\b|\\bconstructed\\b|\\breconstruct', re.I)
SCREEN_RE = re.compile('isolat|spontaneous|recovered|arose|found in|selected|uv screen|screen', re.I)
FULLTEXT = ('full text', 'preprint full text')

def read_table(path):
    with open(path, encoding='utf-8', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def pdr5_residue(row):
    if row['protein'] == 'Pdr5':
        v = row['residue_native']
    else:
        v = row['pdr5_equivalent_residue']
    v = (v or '').strip()
    return int(v) if v.isdigit() else None

def _flags(row):
    return set(row['other_flags'].replace(';', ' ').split())

def _base_ok(row, allow_identity_secondary=False):
    fl = _flags(row)
    bad = {'SECONDARY'} | (set() if allow_identity_secondary else {'IDENTITY_FROM_SECONDARY'})
    return row['evidence_use'] == 'single' and row['expr_loc_defect_flag'].strip() != 'Y' and (not fl & bad) and (row['evidence_access'] in FULLTEXT) and row['residue_verified_vs_UniProt'].strip().upper().startswith('Y')

def eligible(row, tier):
    r = pdr5_residue(row)
    if r is None:
        return False
    pk, cat = (row['paper_key'], row['phenotype_category'])
    designed = bool(DESIGNED_RE.search(row['phenotype_summary']))
    if tier == 'U1':
        return pk in UNBIASED_CORE and _base_ok(row) and (cat == 'altered specificity') and (not designed)
    if tier == 'U2':
        if pk in UNBIASED_CORE:
            return _base_ok(row, allow_identity_secondary=True) and cat == 'altered specificity' and (not designed)
        return pk in UNBIASED_ABSTRACT and row['evidence_use'] == 'single' and (cat == 'altered specificity') and (row['expr_loc_defect_flag'].strip() != 'Y') and ('SECONDARY' not in _flags(row)) and row['residue_verified_vs_UniProt'].strip().upper().startswith('Y') and (not designed)
    if tier == 'UANY':
        return pk in UNBIASED_CORE and _base_ok(row) and (cat in ('altered specificity', 'gain of function', 'loss of function')) and bool(SCREEN_RE.search(row['phenotype_summary'])) and (not designed)
    if tier == 'TAS':
        return pk not in UNBIASED_CORE and pk not in UNBIASED_ABSTRACT and _base_ok(row) and (cat == 'altered specificity')
    if tier == 'TESTED':
        return pk not in UNBIASED_CORE and pk not in UNBIASED_ABSTRACT and _base_ok(row)
    raise ValueError(tier)

def residue_table(rows, tier, protein=None, exclude_paper=None):
    out = {}
    for row in rows:
        if exclude_paper and row['paper_key'] == exclude_paper:
            continue
        if protein and row['protein'] != protein:
            continue
        if not eligible(row, tier):
            continue
        r = pdr5_residue(row)
        d = out.setdefault(r, dict(proteins=set(), papers=set(), n_rows=0))
        d['proteins'].add(row['protein'])
        d['papers'].add(row['paper_key'])
        d['n_rows'] += 1
    return out

def hit_sets(rows, tier):
    pooled = set(residue_table(rows, tier))
    own = set(residue_table(rows, tier, protein='Pdr5'))
    cdr = set(residue_table(rows, tier, protein='Cdr1'))
    papers = sorted({r['paper_key'] for r in rows if eligible(r, tier)})
    loo = {p: set(residue_table(rows, tier, exclude_paper=p)) for p in papers}
    return dict(pooled=pooled, pdr5=own, cdr1=cdr, loo=loo, papers=papers)

def counts_report(rows):
    lines = []
    for tier in ('U1', 'U2', 'UANY', 'TAS', 'TESTED'):
        hs = hit_sets(rows, tier)
        rows_t = [r for r in rows if eligible(r, tier)]
        lines.append(f"tier {tier}: rows={len(rows_t)} residues pooled={len(hs['pooled'])} Pdr5-own={len(hs['pdr5'])} Cdr1-mapped={len(hs['cdr1'])} both={len(hs['pdr5'] & hs['cdr1'])}")
        if tier in ('U1', 'U2'):
            c = Counter((r['paper_key'] for r in rows_t))
            res_by_paper = defaultdict(set)
            for r in rows_t:
                res_by_paper[r['paper_key']].add(pdr5_residue(r))
            for p in sorted(c):
                lines.append(f'    {p}: rows={c[p]} residues={len(res_by_paper[p])}')
    dropped = Counter()
    for r in rows:
        if r['paper_key'] in UNBIASED_CORE and _base_ok(r) and (r['phenotype_category'] == 'altered specificity') and (pdr5_residue(r) is not None) and DESIGNED_RE.search(r['phenotype_summary']):
            dropped[r['paper_key']] += 1
    lines.append("rows removed by the 'designed' rule (otherwise eligible for U1): " + (str(dict(dropped)) or 'none'))
    unmapped = sum((1 for r in rows if r['protein'] == 'Cdr1' and pdr5_residue(r) is None))
    lines.append(f'Cdr1 rows without a Pdr5 equivalent (unmappable): {unmapped}')
    return '\n'.join(lines)

def write_tables(rows, outdir):
    import os
    os.makedirs(outdir, exist_ok=True)
    for tier in ('U1', 'U2', 'UANY', 'TAS', 'TESTED'):
        t = residue_table(rows, tier)
        with open(os.path.join(outdir, f'a12_hits_{tier}.tsv'), 'w', encoding='utf-8', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['residue', 'proteins', 'papers', 'n_rows'])
            for r in sorted(t):
                w.writerow([r, ';'.join(sorted(t[r]['proteins'])), ';'.join(sorted(t[r]['papers'])), t[r]['n_rows']])
if __name__ == '__main__':
    import sys
    rows = read_table(sys.argv[1])
    print(counts_report(rows))
    if len(sys.argv) > 2:
        write_tables(rows, sys.argv[2])
        print('tables written to', sys.argv[2])

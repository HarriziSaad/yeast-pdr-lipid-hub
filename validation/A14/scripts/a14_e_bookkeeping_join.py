import argparse
import csv
import os
import re
import sys
import openpyxl
TIP_RE = re.compile('[(,]\\s*([^:(),;\\[\\]]+?)\\s*(?=:)')

def read_s2(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb['Table S2']
    rows = list(ws.iter_rows(values_only=True))
    header = list(rows[0])
    ix_sp = header.index('species')
    ix_as = header.index('assembly')
    drug_cols = [(i, h) for i, h in enumerate(header[:12]) if h and h.endswith('_resistance')]
    data = []
    for r in rows[1:]:
        if r is None or r[ix_as] is None:
            continue
        a = str(r[ix_as]).strip()
        s = str(r[ix_sp]).strip() if r[ix_sp] is not None else ''
        if not a or a.startswith('=') or s.startswith('='):
            continue
        data.append((s, a, r))
    return (header, ix_sp, ix_as, drug_cols, data)

def read_names(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    h = list(rows[0])
    out = {}
    for r in rows[1:]:
        if r[0]:
            out[str(r[0]).strip()] = {h[j]: str(r[j]).strip() if r[j] is not None else '' for j in range(1, len(h))}
    return (h, out)

def read_tree_tips(path):
    s = open(path, encoding='utf-8').read()
    tips = [t for t in TIP_RE.findall(s)]
    return tips

def read_tsv(path, cols):
    out = {}
    with open(path, encoding='utf-8', newline='') as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        for r in rd:
            out[r[cols[0]]] = r
        return (rd.fieldnames, out)

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--s2', required=True)
    ap.add_argument('--names', required=True)
    ap.add_argument('--tree', required=True)
    ap.add_argument('--qc', required=True, help='phase1a qc_summary_by_genome.tsv')
    ap.add_argument('--index', required=True, help='data_staging y1000_proteome_index.tsv')
    ap.add_argument('--s3', default=None, help='Harrison 2024 S3 TSV (cycloheximide columns; counts only)')
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    log = []

    def say(s=''):
        print(s)
        log.append(s)
    header, ix_sp, ix_as, drug_cols, data = read_s2(a.s2)
    assemblies = [d[1] for d in data]
    say('S2 data rows (non-formula, with an assembly id): %d; distinct assembly ids: %d; distinct species strings: %d' % (len(data), len(set(assemblies)), len(set((d[0] for d in data)))))
    say('S2 drug columns: %s' % ', '.join((h for _, h in drug_cols)))
    for i, h in drug_cols:
        vals = [d[2][i] for d in data]
        good = [v for v in vals if str(v).strip() in ('0', '1', '0.0', '1.0')]
        ones = [v for v in good if str(v).strip() in ('1', '1.0')]
        say('  %-28s non-missing 0/1 calls %d, of which 1: %d (marginal only)' % (h, len(good), len(ones)))
    nh, names = read_names(a.names)
    say('ALL_NAMES rows: %d; columns: %s' % (len(names), ', '.join(nh)))
    tips = read_tree_tips(a.tree)
    tipset = set(tips)
    say('time tree tips: %d (distinct %d)' % (len(tips), len(tipset)))
    qc_fields, qc = read_tsv(a.qc, ['species_code'])
    idx_fields, idx = read_tsv(a.index, ['speciesCode'])
    say('Phase 1a qc_summary_by_genome rows: %d (columns include n_PASS: %s)' % (len(qc), 'n_PASS' in qc_fields))
    say('proteome index rows: %d (columns include n_proteins: %s)' % (len(idx), 'n_proteins' in idx_fields))
    cand = ['Tip_ID', 'Tip_ID_v2', 'FINAL_Tip_ID']
    in_names = [x for x in assemblies if x in names]
    say('S2 assemblies found in ALL_NAMES (assembly_fullID_updated): %d of %d' % (len(in_names), len(assemblies)))
    best = None
    for c in cand:
        tipm = [x for x in in_names if names[x].get(c, '') in tipset]
        say('  via %-13s -> tree tip: %d' % (c, len(tipm)))
        qcm = [x for x in in_names if names[x].get(c, '') in qc]
        say('  via %-13s -> Phase 1a qc_summary_by_genome species_code: %d' % (c, len(qcm)))
        idxm = [x for x in in_names if names[x].get(c, '') in idx]
        say('  via %-13s -> proteome index speciesCode: %d' % (c, len(idxm)))
    idx_by_asm = {}
    for k, r in idx.items():
        idx_by_asm[r.get('assembly_fullID_updated', '')] = r
    say('proteome index rows keyed by assembly_fullID_updated: %d; S2 assemblies found there: %d' % (len(idx_by_asm), sum((1 for x in assemblies if x in idx_by_asm))))
    rows = []
    n_tip = n_qc = n_both = 0
    for s, asm, _ in data:
        nm = names.get(asm, {})
        tip = ''
        for c in ('FINAL_Tip_ID', 'Tip_ID', 'Tip_ID_v2'):
            if nm.get(c, '') in tipset:
                tip = nm[c]
                tip_col = c
                break
        ir = idx_by_asm.get(asm) or idx_by_asm.get(nm.get('assembly_fullID_updated', ''))
        qcode = ''
        if ir:
            qcode = ir['speciesCode'] if ir['speciesCode'] in qc else ''
        rows.append((s, asm, tip, qcode))
        n_tip += 1 if tip else 0
        n_qc += 1 if qcode else 0
        n_both += 1 if tip and qcode else 0
    say('RESULT: S2 rows %d; with a time-tree tip %d; with a Phase 1a genome row %d; with both %d' % (len(rows), n_tip, n_qc, n_both))
    dup_tip = len(rows) - len(set((r[2] for r in rows if r[2])))
    say('  duplicated tree tips among the matched rows: %d' % dup_tip)
    if a.s3:
        s2_tips = set((r[2] for r in rows if r[2]))
        with open(a.s3, encoding='utf-8', newline='') as fh:
            rd = csv.DictReader(fh, delimiter='\t')
            s3 = list(rd)
        say('S3 (Harrison 2024) rows: %d' % len(s3))
        chx_cols = [c for c in s3[0].keys() if 'ycloheximide' in c]
        for c in chx_cols:
            sc = [r for r in s3 if str(r[c]).strip() in ('0', '1')]
            on_tree = [r for r in sc if r['species name'] in tipset]
            in_qc = [r for r in sc if r['species name'] in qc]
            in_s2 = [r for r in sc if r['species name'] in s2_tips]
            say('  %-24s scored %d; on time tree %d; in Phase 1a %d; also in S2 (532) %d (counts only)' % (c, len(sc), len(on_tree), len(in_qc), len(in_s2)))
        both = [r for r in s3 if all((str(r[c]).strip() in ('0', '1') for c in chx_cols))]
        say('  both CHX concentrations scored: %d (on tree %d)' % (len(both), sum((1 for r in both if r['species name'] in tipset))))
    with open(os.path.join(a.outdir, 'e_join_table.tsv'), 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['species_S2', 'assembly', 'tree_tip', 'phase1a_species_code'])
        w.writerows(rows)
    with open(os.path.join(a.outdir, 'e_bookkeeping_log.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(log) + '\n')
if __name__ == '__main__':
    main()

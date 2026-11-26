import argparse
import csv
import gzip
import hashlib
import io
import os
import re
import statistics
import sys
import tarfile
GENES = [('YOR153W', 'PDR5'), ('YDR011W', 'SNQ2'), ('YDR406W', 'PDR15'), ('YOR328W', 'PDR10'), ('YPL058C', 'PDR12'), ('YIL013C', 'PDR11'), ('YNR070W', 'PDR18'), ('YOR011W', 'AUS1')]
HDR_RE = re.compile('^>(?P<iso>.+?)_(?P<orf>Y[A-P][LR]\\d{3}[WC](?:-[A-Z])?)_(?P<rest>.*)$')

def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()

def parse_fasta_text(text):
    hdr, buf = (None, [])
    for line in text.splitlines():
        if line.startswith('>'):
            if hdr is not None:
                yield (hdr, ''.join(buf))
            hdr, buf = (line.rstrip(), [])
        else:
            buf.append(line.strip())
    if hdr is not None:
        yield (hdr, ''.join(buf))

def read_strain_table(path):
    names = []
    with open(path, encoding='utf-8', errors='replace', newline='') as fh:
        rd = csv.reader(fh, delimiter='\t')
        header = next(rd)
        for row in rd:
            if row and row[0].strip():
                names.append(row[0].strip())
    return (header, names)

def read_first_column(path):
    ids = []
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt', encoding='utf-8', errors='replace', newline='') as fh:
        rd = csv.reader(fh, delimiter='\t')
        next(rd)
        for row in rd:
            if row:
                ids.append(row[0].strip())
    return ids

def read_presence_columns(path, orfs):
    out = {}
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace', newline='') as fh:
        rd = csv.reader(fh, delimiter='\t')
        header = next(rd)
        idx = {}
        for orf in orfs:
            cols = [i for i, h in enumerate(header) if h.split('.')[-1] == orf or h.endswith('.' + orf) or '.' + orf in h]
            idx[orf] = cols
        counts = {orf: [0] * len(idx[orf]) for orf in orfs}
        n = 0
        for row in rd:
            if not row:
                continue
            n += 1
            for orf in orfs:
                for k, i in enumerate(idx[orf]):
                    try:
                        if float(row[i]) > 0:
                            counts[orf][k] += 1
                    except (ValueError, IndexError):
                        pass
        for orf in orfs:
            out[orf] = [(header[i], counts[orf][k], n) for k, i in enumerate(idx[orf])]
    return out

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--tar', required=True)
    ap.add_argument('--extracted-dir', required=True)
    ap.add_argument('--strains', required=True)
    ap.add_argument('--pheno', required=True)
    ap.add_argument('--presence', required=True)
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    log = []

    def say(s=''):
        print(s)
        log.append(s)
    want = {'Proteome_1011/%s.fasta' % g: g for g, _ in GENES}
    members_total = 0
    tar_texts = {}
    tar_sha = {}
    with tarfile.open(a.tar, 'r:gz') as tf:
        for m in tf:
            if not m.isfile():
                continue
            members_total += 1
            if m.name in want:
                data = tf.extractfile(m).read()
                tar_sha[want[m.name]] = sha256_bytes(data)
                tar_texts[want[m.name]] = data.decode('ascii', errors='replace')
    say('tar members (files): %d; PDR gene members found: %d of %d' % (members_total, len(tar_texts), len(GENES)))
    st_header, strain_names = read_strain_table(a.strains)
    strain_set = set(strain_names)
    pheno_ids = read_first_column(a.pheno)
    pheno_set = set(pheno_ids)
    say('strain table rows (non-empty Standardized name): %d (distinct %d); phenotype-matrix row ids: %d (distinct %d)' % (len(strain_names), len(strain_set), len(pheno_ids), len(pheno_set)))
    pres = read_presence_columns(a.presence, [g for g, _ in GENES])
    rows = []
    iso_by_gene = {}
    for orf, label in GENES:
        text = tar_texts.get(orf)
        if text is None:
            say('MISSING in archive: %s' % orf)
            continue
        ext_path = os.path.join(a.extracted_dir, orf + '.fasta')
        ext_ok = 'NA'
        if os.path.exists(ext_path):
            with open(ext_path, 'rb') as fh:
                ext_ok = 'identical' if sha256_bytes(fh.read()) == tar_sha[orf] else 'DIFFERENT'
        recs = list(parse_fasta_text(text))
        n = len(recs)
        parsed, isos, has_coord, lens, alpha = (0, [], 0, [], set())
        nonterm_star = False
        term_star = 0
        for hdr, seq in recs:
            mo = HDR_RE.match(hdr)
            if mo and mo.group('orf') == orf:
                parsed += 1
                isos.append(mo.group('iso'))
                if re.search('chromosome\\w+:\\d+-\\d+', mo.group('rest')):
                    has_coord += 1
            alpha.update(seq)
            if seq.endswith('*'):
                term_star += 1
            if '*' in seq[:-1]:
                nonterm_star = True
            lens.append(len(seq.rstrip('*')))
        iso_set = set(isos)
        iso_by_gene[orf] = iso_set
        modal = max(set(lens), key=lens.count) if lens else 0
        n_in_strain = len(iso_set & strain_set)
        n_in_pheno = len(iso_set & pheno_set)
        pcols = pres.get(orf, [])
        pstr = ';'.join(('%s=%d/%d' % c for c in pcols))
        rows.append({'gene': label, 'orf': orf, 'records': n, 'header_parsed': parsed, 'distinct_isolate_codes': len(iso_set), 'max_records_per_isolate': max((isos.count(i) for i in iso_set), default=0) if len(iso_set) < 5000 else 'NA', 'headers_with_coordinates': has_coord, 'isolates_in_strain_table': n_in_strain, 'isolates_in_phenotype_matrix_ids': n_in_pheno, 'len_min_noStop': min(lens) if lens else 0, 'len_median_noStop': statistics.median(lens) if lens else 0, 'len_max_noStop': max(lens) if lens else 0, 'share_at_modal_length': round(lens.count(modal) / len(lens), 4) if lens else 0, 'records_with_terminal_star': term_star, 'any_nonterminal_star': 'yes' if nonterm_star else 'no', 'alphabet': ''.join(sorted(alpha)), 'extracted_copy_vs_archive': ext_ok, 'archive_member_sha256': tar_sha[orf], 'presence_matrix_columns(name=present/rows)': pstr})
    if iso_by_gene:
        allg = set.intersection(*iso_by_gene.values())
        say('isolates with a record for ALL %d genes: %d (of which in phenotype ids: %d)' % (len(iso_by_gene), len(allg), len(allg & pheno_set)))
        anyg = set.union(*iso_by_gene.values())
        say('isolates with a record for >= 1 gene: %d; of these not in strain table: %d; not in phenotype ids: %d' % (len(anyg), len(anyg - strain_set), len(anyg - pheno_set)))
        say('phenotype ids without a PDR5 record: %d; strain-table names without a PDR5 record: %d' % (len(pheno_set - iso_by_gene.get('YOR153W', set())), len(strain_set - iso_by_gene.get('YOR153W', set()))))
        non_strain = sorted(anyg - strain_set)
        say('isolate codes not in the strain table (first 12): %s' % ', '.join(non_strain[:12]))
    out_tsv = os.path.join(a.outdir, 'p_bookkeeping_1011_counts.tsv')
    if rows:
        with open(out_tsv, 'w', encoding='utf-8', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
            w.writeheader()
            w.writerows(rows)
    with open(os.path.join(a.outdir, 'p_bookkeeping_1011_log.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(log) + '\n')
    for r in rows:
        say('%-6s records %4d distinct %4d  in-strain-table %4d  in-pheno-ids %4d  len %d/%s/%d  modal %.3f  nonterm* %s  alphabet %s' % (r['gene'], r['records'], r['distinct_isolate_codes'], r['isolates_in_strain_table'], r['isolates_in_phenotype_matrix_ids'], r['len_min_noStop'], r['len_median_noStop'], r['len_max_noStop'], r['share_at_modal_length'], r['any_nonterminal_star'], r['alphabet']))
    return rows
if __name__ == '__main__':
    main()

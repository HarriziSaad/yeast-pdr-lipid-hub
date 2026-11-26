import argparse
import json
import os
import platform
import sys
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a14_p_lib as P
A11_LIB_SHA256 = 'bb1cfcfd45b792ed0269ea1aaf8b8c30c6bbff110132c46c11a38df12f7346ee'
MAFFT_SHA256 = '77e7a5fddd2365af9854bbaaa25f21ced1a906be97eff3871bd82c4e60d63d25'
MUSCLE_SHA256 = '841f438c9f754c9c1d87bc203e50c9d5b6162cef1b0fe4678db1a295ee142ab3'
MANIFEST_SHA256 = '735b5231a4b9cdde5372dd16f679d96f1664c07269897e8c6e7d166f920c4f85'
FROZEN_SHA256 = 'e409356b6e46c1b0b10fdaa79a1f1b1c835cfbe038b939ee03ceeb8d853fb0a3'
COVARIATES_SHA256 = '5a32e3fd32c93aee4bd0093da69caea657a3fa2f0aa2345c4ede52d0c9644eba'
SEED = 20261118
RUN_TABLE = [('PRIMARY', 'distinct', 'all', 'all', 'primary'), ('S_nosingle', 'nosingle', 'all', 'all', 'primary'), ('S_af', 'af', 'all', 'all', 'primary'), ('S_nostop', 'distinct', 'all', 'all_nostop', 'primary'), ('S_major', 'distinct', 'all', 'major', 'primary'), ('S_notmajor', 'distinct', 'all', 'notmajor', 'primary'), ('S_pdr5only', 'distinct', 'PDR5', 'all', 'primary'), ('C_comp', 'distinct', 'all', 'all', 'comp'), ('C_exposed', 'distinct', 'all', 'all', 'exposed')]

def read_class_sets(frozen_rows, truth_rows):
    h5 = {}
    for r in frozen_rows:
        if r['class'] in ('LH-CONF', 'GATE-core', 'NBS-MOTIF'):
            h5.setdefault(r['class'], set()).add(int(r['residue']))
    truth = {}
    for r in truth_rows:
        truth.setdefault(r['set'], set()).add(int(r['residue']))
    return (h5, truth)

def locate_s288c_rows(ids, seqs, refs):
    out = {}
    for g, ref in refs.items():
        hits = [s for i, s in zip(ids, seqs) if i.startswith('Saccharomyces_cerevisiae|') and s.replace('-', '').replace('.', '') == ref]
        if len(hits) != 1:
            raise RuntimeError(f'S288C row of {g}: {len(hits)} exact matches (expected 1)')
        out[g] = hits[0]
    return out

def build_alignment_context(a11, ids, seqs, refs, cov):
    rows = locate_s288c_rows(ids, seqs, refs)
    colmap = P.build_column_map(rows, 'PDR5')
    gcol = a11.gap_fractions(seqs)
    na = {c for c, g in enumerate(gcol) if g > a11.GAP_NA_THRESHOLD}
    rc = a11.reference_columns(rows['PDR5'])
    gap_res = {r: gcol[rc[r]] for r in cov if r in rc}
    return {'colmap': colmap, 'na': na, 'gap': gap_res, 'n_rows': len(seqs), 'n_cols': len(gcol), 'n_na_cols': len(na)}

def analyse(a11, refs, cds, gene_records, ctx_primary, ctx_second, cov, classes_h5, truth, frozen_res, clusters, n_draws=10000, n_boot=2000, n_rep_boot=200, n_draws_boot=1000, seed=SEED, with_gene_runs=True):
    o = {g: P.o_per_position(cds[g], refs[g]) for g in refs}
    var_all = {g: P.gene_variation(g, refs[g], gene_records[g]) for g in refs}
    var_nostop = {g: P.gene_variation(g, refs[g], gene_records[g], mask_downstream_of_stop=True) for g in refs}
    strains_all = sorted(set().union(*[set(v.diffs) for v in var_all.values()]))
    major = {s for s, c in clusters.items() if c == 1}
    subsets = {'all': var_all, 'all_nostop': var_nostop}
    subsets['major'] = {g: P.subset_variation(v, [s for s in v.diffs if s in major]) for g, v in var_all.items()}
    subsets['notmajor'] = {g: P.subset_variation(v, [s for s in v.diffs if s not in major]) for g, v in var_all.items()}
    exclude_h5 = set(frozen_res)
    exclude_hub = set(frozen_res) | set(truth.get('H6_ALL', set()))
    tests = [(name, res, exclude_h5, 'F5') for name, res in sorted(classes_h5.items())]
    tests.append(('H6_PRIMARY', set(truth['H6_PRIMARY']), exclude_hub, 'F5b'))
    results, density_tables, info = ({}, {}, {'strains': len(strains_all), 'gene_records': {g: v.n_records for g, v in var_all.items()}, 'qc_flags': {g: v.qc_flags for g, v in var_all.items()}})
    for cname, cres, excl, fam in tests:
        universe = sorted((r for r in cov if r not in excl))
        results[cname] = {}
        for run, mode, genes, sub, cmode in RUN_TABLE:
            var = subsets[sub]
            gl = None if genes == 'all' else [genes]
            dens, diag = P.column_density(ctx_primary['colmap'], var, o, mode, gl, ctx_primary['na'])
            if run == 'PRIMARY':
                density_tables[cname] = {r: diag[r] for r in sorted(cres) if r in diag}
            r = P.h5_test(a11, sorted(cres), universe, cov, dens, ctx_primary['gap'], n_draws, seed, n_boot, comp_match=cmode == 'comp', exposed_only=cmode == 'exposed', diag=diag)
            r.update({'family': fam, 'class': cname, 'run': run, 'density_mode': mode, 'genes': genes, 'strains': sub})
            results[cname][run] = r
        if ctx_second is not None:
            dens, diag = P.column_density(ctx_second['colmap'], var_all, o, 'distinct', None, ctx_second['na'])
            r = P.h5_test(a11, sorted(cres), universe, cov, dens, ctx_second['gap'], n_draws, seed, n_boot, diag=diag)
            r.update({'family': fam, 'class': cname, 'run': 'V_muscle', 'density_mode': 'distinct', 'genes': 'all', 'strains': 'all'})
            results[cname]['V_muscle'] = r
        results[cname]['BOOT'] = P.cluster_bootstrap(a11, sorted(cres), universe, cov, ctx_primary['gap'], ctx_primary['colmap'], var_all, o, clusters, ctx_primary['na'], n_rep=n_rep_boot, n_draws=n_draws_boot, seed=seed)
        if with_gene_runs:
            for g in refs:
                dens, diag = P.column_density(ctx_primary['colmap'], var_all, o, 'distinct', [g], ctx_primary['na'])
                r = P.h5_test(a11, sorted(cres), universe, cov, dens, ctx_primary['gap'], min(n_draws, 2000), seed, 500, diag=diag)
                r.update({'family': 'descriptive', 'class': cname, 'run': 'D_gene_' + g, 'genes': g})
                results[cname]['D_gene_' + g] = r
    f5 = [c for c in sorted(classes_h5)]
    ps = [results[c]['PRIMARY'].get('p', 1.0) if results[c]['PRIMARY'].get('status') == 'OK' else 1.0 for c in f5]
    qs = P.bh(ps)
    qmap = {c: q for c, q in zip(f5, qs)}
    qmap['H6_PRIMARY'] = results['H6_PRIMARY']['PRIMARY'].get('p', 1.0) if results['H6_PRIMARY']['PRIMARY'].get('status') == 'OK' else 1.0
    labels = {}
    for cname in qmap:
        labels[cname] = P.label_h5(results[cname], qmap[cname])
        labels[cname]['q'] = qmap[cname]
        flags = results[cname]['PRIMARY'].get('flags', [])
        labels[cname]['power_flags'] = flags
    info['alignment_primary'] = {k: ctx_primary[k] for k in ('n_rows', 'n_cols', 'n_na_cols')}
    if ctx_second is not None:
        info['alignment_second'] = {k: ctx_second[k] for k in ('n_rows', 'n_cols', 'n_na_cols')}
    return (results, labels, density_tables, info)
COLS = ['class', 'family', 'run', 'density_mode', 'genes', 'strains', 'status', 'n_requested', 'n_eff', 'observed_mean', 'control_mean', 'ratio', 'ci95_ratio_lo', 'ci95_ratio_hi', 'delta', 'z', 'p', 'b', 'n_fallback_residues', 'weakly_matched', 'class_sum_m', 'flags', 'share_ratio_lt_1', 'n_replicates', 'ratio_median', 'dropped_NA']

def write_outputs(out_dir, results, labels, density_tables, info, manifest):
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, 'A14_H5_results.tsv'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\t'.join(COLS) + '\n')
        for cname, runs in results.items():
            for run, r in runs.items():
                row = dict(r)
                row['class'], row['run'] = (cname, run)
                ci = row.get('ci95_ratio', ('', ''))
                row['ci95_ratio_lo'], row['ci95_ratio_hi'] = ci
                row['flags'] = ','.join(row.get('flags', []))
                row['dropped_NA'] = ','.join((str(x) for x in row.get('dropped_NA', [])))
                fh.write('\t'.join((str(row.get(c, '')) for c in COLS)) + '\n')
    with open(os.path.join(out_dir, 'A14_H5_density_primary.tsv'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('class\tresidue\tsum_m\tsum_o\tn_genes\n')
        for cname, d in density_tables.items():
            for r, (m, oo, ng) in d.items():
                fh.write(f'{cname}\t{r}\t{m}\t{oo}\t{ng}\n')
    with open(os.path.join(out_dir, 'A14_H5_label.json'), 'w', encoding='utf-8') as fh:
        json.dump(labels, fh, indent=2)
    manifest['info'] = info
    manifest['end_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    with open(os.path.join(out_dir, 'A14_H5_manifest.json'), 'w', encoding='utf-8') as fh:
        json.dump(manifest, fh, indent=2, default=str)

def load_refs_and_cds(seeds_fasta, pangenome):
    refs = {}
    for hdr, seq in P.read_fasta_file(seeds_fasta):
        parts = hdr.split()[0].split('|')
        lab = parts[2].split('_')[0] if len(parts) >= 3 else hdr.split()[0]
        if lab in P.GENE_LABELS:
            refs[lab] = seq
    missing = [g for g in P.GENE_LABELS if g not in refs]
    if missing:
        raise RuntimeError(f'reference proteins missing: {missing}')
    cds = {}
    import re
    rx = re.compile('^\\d+-(Y[A-P][LR]\\d{3}[WC](?:-[A-Z])?)')
    sys2lab = {v: k for k, v in P.SYSTEMATIC.items()}
    for hdr, seq in P.read_fasta_file(pangenome):
        m = rx.match(hdr.split()[0])
        if m and m.group(1) in sys2lab and (sys2lab[m.group(1)] not in cds):
            cds[sys2lab[m.group(1)]] = seq
    missing = [g for g in P.GENE_LABELS if g not in cds]
    if missing:
        raise RuntimeError(f'reference ORFs missing: {missing}')
    return (refs, cds)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', required=True, help='project directory')
    ap.add_argument('--proteome-dir', required=True, help='folder with <SYSTEMATIC>.fasta of the 8 genes')
    ap.add_argument('--seeds-fasta', required=True, help='S288C UniProt proteins of the 8 genes')
    ap.add_argument('--pangenome', required=True, help='allORFs_pangenome.fasta.gz')
    ap.add_argument('--aln-mafft', required=True)
    ap.add_argument('--aln-muscle', required=True)
    ap.add_argument('--distance-matrix', required=True, help='1011DistanceMatrixBasedOnSNPs.tab.gz')
    ap.add_argument('--truth-sets', required=True, help='a11_truth_sets.tsv')
    ap.add_argument('--a11-lib', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--n-clusters', type=int, default=P.N_CLUSTERS_DEFAULT)
    ap.add_argument('--n-draws', type=int, default=10000)
    ap.add_argument('--check-only', action='store_true', help='integrity checks only; no density or test is computed')
    ap.add_argument('--toy', action='store_true', help='skip the real-data hash checks (synthetic inputs)')
    a = ap.parse_args()
    a11 = P.load_module(a.a11_lib, 'a11_lib', None if a.toy else A11_LIB_SHA256)
    R = a.root
    p_man = os.path.join(R, 'phase1_orthology', 'phase1b', 'outputs', 'set', 'analysis_set_manifest.tsv')
    p_cov = os.path.join(R, 'preregistration', 'outputs', 'structure_covariates_7P04.tsv')
    p_frz = os.path.join(R, 'preregistration', 'site_classes_FROZEN.tsv')
    if not a.toy:
        for p, want in ((p_man, MANIFEST_SHA256), (p_cov, COVARIATES_SHA256), (p_frz, FROZEN_SHA256), (a.aln_mafft, MAFFT_SHA256), (a.aln_muscle, MUSCLE_SHA256)):
            got = P.sha256_file(p)
            if got != want:
                sys.exit(f'HASH MISMATCH {p}: {got}')
    cov = a11.load_covariates(p_cov)
    refs, cds = load_refs_and_cds(a.seeds_fasta, a.pangenome)
    bad = a11.check_reference(refs['PDR5'], {r: v['aa'] for r, v in cov.items()})
    if bad:
        sys.exit('PDR5 reference differs from the structure covariates: ' + '; '.join(bad[:5]))
    frozen_rows = a11.read_tsv(p_frz)
    truth_rows = a11.read_tsv(a.truth_sets)
    classes_h5, truth = read_class_sets(frozen_rows, truth_rows)
    frozen_res = {int(r['residue']) for r in frozen_rows}
    man = [r['seq_id'] for r in a11.read_tsv(p_man) if r['role'] == 'INGROUP']
    if not a.toy and len(man) != 1253:
        sys.exit(f'unexpected INGROUP row count {len(man)} (expect 1253)')
    ctxs = []
    for path in (a.aln_mafft, a.aln_muscle):
        ids, seqs = a11.read_fasta(path)
        sids, sseq = a11.select_rows(ids, seqs, man)
        ctxs.append(build_alignment_context(a11, sids, sseq, refs, cov))
        print(path, 'rows', len(sids), 'columns', ctxs[-1]['n_cols'], 'NA columns', ctxs[-1]['n_na_cols'])
    gene_records = {}
    for g, sysn in P.GENES:
        gene_records[g] = P.read_fasta_file(os.path.join(a.proteome_dir, sysn + '.fasta'))
        print(g, 'records', len(gene_records[g]))
    ids, D = P.read_distance_matrix(a.distance_matrix)
    clusters = P.cluster_strains(ids, D, a.n_clusters)
    if a.check_only:
        print('check-only: inputs readable, S288C rows located, clusters:', len(set(clusters.values())))
        return
    manifest = {'python': sys.version, 'platform': platform.platform(), 'start_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'seed': SEED, 'n_draws': a.n_draws, 'n_clusters': a.n_clusters, 'inputs': {os.path.basename(p): P.sha256_file(p) for p in (a.aln_mafft, a.aln_muscle, p_man, p_cov, p_frz, a.truth_sets, a.a11_lib, a.distance_matrix, a.seeds_fasta, a.pangenome)}}
    res, labels, dens, info = analyse(a11, refs, cds, gene_records, ctxs[0], ctxs[1], cov, classes_h5, truth, frozen_res, clusters, n_draws=a.n_draws)
    write_outputs(a.out_dir, res, labels, dens, info, manifest)
    for c, l in labels.items():
        print(c, l['label'], 'q =', round(l['q'], 4), 'failed:', l['failed'], 'flags:', l['power_flags'])
if __name__ == '__main__':
    main()

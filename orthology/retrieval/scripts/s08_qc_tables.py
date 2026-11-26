import argparse, collections, hashlib, json, os, subprocess, tempfile
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from common import read_tsv, write_tsv, write_fasta, fasta_iter, log
ap = argparse.ArgumentParser()
ap.add_argument('--classification', required=True)
ap.add_argument('--retrans-log', required=True)
ap.add_argument('--final-faa', required=True)
ap.add_argument('--final-summary', required=True)
ap.add_argument('--species-table', required=True)
ap.add_argument('--outdir', required=True)
ap.add_argument('--threads', type=int, default=10)
ap.add_argument('--cdhit', default='cd-hit')
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
cls = read_tsv(a.classification)
rt = {r['seq_id']: r for r in read_tsv(a.retrans_log)}
fin = {r['seq_id']: r for r in read_tsv(a.final_summary)}
seqs = {sid: s for sid, _, s in fasta_iter(a.final_faa)}
st = {r['speciesCode']: r for r in read_tsv(a.species_table)}

def status_of(sid):
    r = fin.get(sid)
    return r['status'] if r else 'no_domain_hit'
rows = []
for r in cls:
    sid = r['seq_id']
    code = r['species_code']
    sp = st.get(code, {})
    t = rt.get(sid, {})
    s = seqs.get(sid, '')
    body = s
    flags = []
    n_stop = body.count('*')
    if n_stop:
        flags.append(f'INTERNAL_STOP({n_stop})')
    n_x = body.count('X')
    if n_x:
        flags.append('X_CONTENT(%.4f)' % (n_x / max(1, len(body))))
        if n_x / max(1, len(body)) > 0.02:
            flags.append('HIGH_X')
    if body and body[0] != 'M':
        flags.append('NO_START_M')
    L = len(body)
    ends = [int(x.split('-')[1]) for x in r['coords_PF01061'].split(';') if x] if r.get('coords_PF01061') else []
    tail = ''
    if ends:
        tail = L - max(ends)
        if tail < 100:
            flags.append(f'CTERM_TAIL_LT100({tail})')
    if r['retrieval_class'] == 'PASS':
        if 1100 <= L <= 1199:
            flags.append('LEN_NEAR_LOWER')
        if 1701 <= L <= 1800:
            flags.append('LEN_NEAR_UPPER')
    if int(r['n_PF00005'] or 0) >= 3 or int(r['n_PF01061'] or 0) >= 3:
        flags.append('FUSION_SUSPECT')
    for f in (r.get('rule_flags') or '').split(';'):
        if f:
            flags.append(f)
    ts = t.get('status', '')
    if t.get('cug_pending') == 'Y':
        flags.append('CUG_PENDING')
    if ts == 'MISMATCH_PEP_CDS':
        flags.append('RETRANS_MISMATCH')
    if ts == 'RETRANSLATED_MASKED_X':
        flags.append('MASKED_X_CTG')
    if sp.get('kind') == 'reference_annotation':
        flags.append('REFERENCE_ANNOTATION_NOT_A_GENOME')
    fs = status_of(sid)
    if fs != r['retrieval_class'] and (not (r['retrieval_class'].startswith('OG_FAIL') and fs in ('no_domain_hit', 'OTHER_ABC'))):
        flags.append('STATUS_CHANGED_AFTER_RETRANSLATION')
    nm = r['retrieval_class'].split(':', 1)[1] if ':' in r['retrieval_class'] else ''
    final_status = fs
    rows.append(dict(species_code=code, FINAL_Tip_ID=sp.get('FINAL_Tip_ID', ''), species_key=sp.get('species_key', ''), order=sp.get('order', ''), genome_kind=sp.get('kind', ''), pep_file=sp.get('pep_file', ''), gene_id=r['gene_id'], seq_id=sid, length_pep=t.get('len_pep', r['seq_len']), length_final=L, n_PF00005=r['n_PF00005'], n_PF01061=r['n_PF01061'], n_PF06422=r['n_PF06422'], n_PF14510=r['n_PF14510'], n_PF19055=r['n_PF19055'], coords_PF00005=r['coords_PF00005'] if 'coords_PF00005' in r else '', coords_PF01061=r.get('coords_PF01061', ''), coords_PF06422=r.get('coords_PF06422', ''), coords_PF14510=r.get('coords_PF14510', ''), coords_PF19055=r.get('coords_PF19055', ''), arch=r['arch'], PF14510_present=r['PF14510_present'], retrieval_class=r['retrieval_class'], final_class=final_status, passes_rule='Y' if r['retrieval_class'] == 'PASS' else 'N', near_miss_class=nm, near_miss_reasons=r['near_miss_reasons'], repair_candidate=r.get('repair_candidate', ''), passes_rule_final_seq='Y' if final_status == 'PASS' else 'N', in_OG0000012=r['in_OG0000012'], OG_id=r['OG_id'], ctg_call=sp.get('ctg_call', ''), CUG_code=sp.get('ctg_class', ''), retranslated=t.get('retranslated', 'N'), retrans_status=ts, retrans_table=t.get('table', ''), n_inframe_CTG=t.get('n_inframe_CTG', ''), n_ctg_positions_changed=t.get('n_changed', ''), tail_after_last_PF01061=tail, qc_flags=';'.join(flags), n_genomes_of_species=sp.get('n_genomes_of_species', ''), seq_sha256=hashlib.sha256(body.encode()).hexdigest()))
pass_rows = [r for r in rows if r['passes_rule_final_seq'] == 'Y' and r['genome_kind'] == 'y1000_genome']
log('PASS rows (final sequences, genomes)', len(pass_rows), 'of', len(rows), 'rows;', sum((1 for r in rows if r['retrieval_class'] == 'PASS' and r['passes_rule_final_seq'] != 'Y')), 'PASS (pep) rows no longer pass after re-translation;', sum((1 for r in rows if r['retrieval_class'] != 'PASS' and r['passes_rule_final_seq'] == 'Y')), 'non-PASS (pep) rows that pass after re-translation')
with open(os.path.join(a.outdir, 'retrieved_pass_final.faa'), 'w') as fh:
    for r in pass_rows:
        write_fasta(fh, r['seq_id'], seqs[r['seq_id']], 'len=%d' % len(seqs[r['seq_id']]))

def run_cdhit(fasta_in, prefix):
    out = prefix + '.nr'
    cmd = [a.cdhit, '-i', fasta_in, '-o', out, '-c', '0.99', '-n', '5', '-G', '1', '-d', '0', '-M', '0', '-T', str(a.threads)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit('cd-hit failed: ' + p.stderr[-500:])
    clusters, cur, rep = (collections.OrderedDict(), None, None)
    for line in open(out + '.clstr'):
        if line.startswith('>Cluster'):
            cur = int(line.split()[1])
            clusters[cur] = []
        else:
            sid = line.split('>', 1)[1].split('...')[0]
            clusters[cur].append(sid)
            if line.rstrip().endswith('*'):
                rep = sid
                clusters[cur].insert(0, clusters[cur].pop())
    return (clusters, out)
tmp = tempfile.mkdtemp(prefix='cdhit_', dir=a.outdir)
clusters_g, nrfile = run_cdhit(os.path.join(a.outdir, 'retrieved_pass_final.faa'), os.path.join(tmp, 'global'))
rep_g = {}
cl_g = {}
for cid, mem in clusters_g.items():
    for m in mem:
        rep_g[m] = mem[0]
        cl_g[m] = cid
with open(os.path.join(a.outdir, 'retrieved_pass_nr99.faa'), 'w') as fh:
    for sid, _, s in fasta_iter(nrfile):
        write_fasta(fh, sid, s)
by_sp = collections.defaultdict(list)
for r in pass_rows:
    by_sp[r['species_key']].append(r['seq_id'])
rep_s = {}
for sp, ids in by_sp.items():
    if len(ids) < 2:
        for i in ids:
            rep_s[i] = i
        continue
    fa = os.path.join(tmp, 'sp.faa')
    with open(fa, 'w') as fh:
        for i in ids:
            write_fasta(fh, i, seqs[i])
    cl, _ = run_cdhit(fa, os.path.join(tmp, 'sp'))
    for cid, mem in cl.items():
        for m in mem:
            rep_s[m] = mem[0]
code_of = {r['seq_id']: r for r in pass_rows}
clus_rows = []
rel = collections.Counter()
for sid in [r['seq_id'] for r in pass_rows]:
    rp = rep_g[sid]
    if rp == sid:
        relation = 'representative' if sum((1 for m in clusters_g[cl_g[sid]])) > 1 else 'singleton'
    else:
        a1, b1 = (code_of[sid], code_of[rp])
        if a1['species_code'] == b1['species_code']:
            relation = 'collapsed_same_genome'
        elif a1['species_key'] == b1['species_key']:
            relation = 'collapsed_same_species_other_strain'
        else:
            relation = 'collapsed_other_species'
    rel[relation] += 1
    clus_rows.append(dict(seq_id=sid, species_code=code_of[sid]['species_code'], species_key=code_of[sid]['species_key'], cluster_global=cl_g[sid], rep_global=rp, is_rep_global='Y' if rp == sid else 'N', relation=relation, rep_within_species=rep_s[sid], is_rep_within_species='Y' if rep_s[sid] == sid else 'N', length=code_of[sid]['length_final']))
write_tsv(os.path.join(a.outdir, 'cdhit99_clusters.tsv'), clus_rows)
cmap = {r['seq_id']: r for r in clus_rows}
for r in rows:
    c = cmap.get(r['seq_id'])
    r['cdhit99_rep_global'] = c['is_rep_global'] if c else 'NA'
    r['cdhit99_rep_id'] = c['rep_global'] if c else ''
    r['cdhit99_relation'] = c['relation'] if c else 'NA'
    r['cdhit99_rep_within_species'] = c['is_rep_within_species'] if c else 'NA'
order = ['species_code', 'FINAL_Tip_ID', 'species_key', 'order', 'genome_kind', 'pep_file', 'gene_id', 'seq_id', 'length_pep', 'length_final', 'n_PF00005', 'n_PF01061', 'n_PF06422', 'n_PF14510', 'n_PF19055', 'coords_PF00005', 'coords_PF01061', 'coords_PF06422', 'coords_PF14510', 'coords_PF19055', 'arch', 'PF14510_present', 'retrieval_class', 'final_class', 'passes_rule', 'passes_rule_final_seq', 'near_miss_class', 'near_miss_reasons', 'repair_candidate', 'in_OG0000012', 'OG_id', 'ctg_call', 'CUG_code', 'retranslated', 'retrans_status', 'retrans_table', 'n_inframe_CTG', 'n_ctg_positions_changed', 'tail_after_last_PF01061', 'qc_flags', 'cdhit99_rep_global', 'cdhit99_rep_id', 'cdhit99_relation', 'cdhit99_rep_within_species', 'n_genomes_of_species', 'seq_sha256']
rows.sort(key=lambda r: (r['order'], r['species_code'], r['gene_id']))
write_tsv(os.path.join(a.outdir, 'retrieval_table.tsv'), rows, order)
write_tsv(os.path.join(a.outdir, 'near_misses.tsv'), [r for r in rows if r['retrieval_class'] != 'PASS'], order)
sp_genomes = collections.defaultdict(set)
for r in st.values():
    if r['kind'] == 'y1000_genome':
        sp_genomes[r['species_key']].add(r['speciesCode'])
dup_rows = []
for sp, gs in sorted(sp_genomes.items()):
    if len(gs) < 2:
        continue
    per = {g: [r for r in pass_rows if r['species_code'] == g] for g in gs}
    for g in sorted(gs):
        mine = per[g]
        n_collapsed_other = sum((1 for r in mine if cmap[r['seq_id']]['relation'] in ('collapsed_same_species_other_strain',) or (cmap[r['seq_id']]['is_rep_global'] == 'Y' and any((m != r['seq_id'] and code_of[m]['species_code'] != g and (code_of[m]['species_key'] == sp) for m in clusters_g[cmap[r['seq_id']]['cluster_global']])))))
        dup_rows.append(dict(species_key=sp, speciesCode=g, genome_in_input='Y' if per[g] or g in st else 'N', n_genomes_of_species=len(gs), n_PASS=len(mine), n_PASS_sharing_99pct_cluster_with_other_strain=n_collapsed_other, fraction_shared=round(n_collapsed_other / len(mine), 3) if mine else ''))
write_tsv(os.path.join(a.outdir, 'strain_duplication.tsv'), dup_rows)
by_g = collections.defaultdict(lambda: collections.Counter())
for r in rows:
    k = r['species_code']
    by_g[k]['n_' + ('PASS' if r['retrieval_class'] == 'PASS' else 'NEAR_MISS' if r['retrieval_class'].startswith('NEAR_MISS') else 'OG_FAIL')] += 1
    if r['retrieval_class'] == 'PASS':
        by_g[k]['n_PASS_in_OG'] += r['in_OG0000012'] == 'Y'
        by_g[k]['n_PASS_nr99_rep'] += r['cdhit99_rep_global'] == 'Y'
        by_g[k]['n_PASS_with_flag'] += bool(r['qc_flags'])
write_tsv(os.path.join(a.outdir, 'qc_summary_by_genome.tsv'), [dict(species_code=k, FINAL_Tip_ID=st.get(k, {}).get('FINAL_Tip_ID', ''), order=st.get(k, {}).get('order', ''), CUG_code=st.get(k, {}).get('ctg_class', ''), **{kk: v for kk, v in c.items()}) for k, c in sorted(by_g.items())])
by_o = collections.defaultdict(collections.Counter)
genomes_o = collections.defaultdict(set)
for k, c in by_g.items():
    o = st.get(k, {}).get('order', '') or ('reference_annotation' if st.get(k, {}).get('kind') == 'reference_annotation' else '')
    genomes_o[o].add(k)
    by_o[o].update(c)
write_tsv(os.path.join(a.outdir, 'qc_summary_by_order.tsv'), [dict(order=o, n_genomes_with_hits=len(genomes_o[o]), **c) for o, c in sorted(by_o.items())])
flagc = collections.Counter((f.split('(')[0] for r in rows if r['retrieval_class'] == 'PASS' for f in r['qc_flags'].split(';') if f))
summ = dict(n_rows=len(rows), n_PASS_pep=sum((1 for r in rows if r['retrieval_class'] == 'PASS')), n_PASS_final_seq=len(pass_rows), n_NEAR_MISS=sum((1 for r in rows if r['retrieval_class'].startswith('NEAR_MISS'))), n_OG_FAIL=sum((1 for r in rows if r['retrieval_class'].startswith('OG_FAIL'))), cdhit99_relations=dict(rel), n_nr99_representatives=sum((1 for r in clus_rows if r['is_rep_global'] == 'Y')), n_within_species_representatives=sum((1 for r in clus_rows if r['is_rep_within_species'] == 'Y')), PASS_flag_counts=dict(flagc), retranslation_status=dict(collections.Counter((r['retrans_status'] for r in rows))))
json.dump(summ, open(os.path.join(a.outdir, 'qc_summary.json'), 'w'), indent=1)
log(json.dumps(summ))

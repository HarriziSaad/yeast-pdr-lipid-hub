import argparse, collections, csv, hashlib, json, math, os, re
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from common import read_tsv, write_tsv, species_key_of, sanitize, log, SEED
YGOB_POST = {'Vanderwaltozyma polyspora': 'Vanderwaltozyma_polyspora', 'Tetrapisispora phaffii': 'Tetrapisispora_phaffii', 'Tetrapisispora blattae': 'Tetrapisispora_blattae', 'Naumovozyma dairenensis': 'Naumovozyma_dairenensis', 'Naumovozyma castellii': 'Naumovozyma_castellii', 'Kazachstania naganishii': 'Kazachstania_naganishii', 'Kazachstania africana': 'Kazachstania_africana', 'Candida glabrata (= Nakaseomyces glabratus)': 'Nakaseomyces_glabratus', 'Saccharomyces bayanus var. uvarum (= S. uvarum)': 'Saccharomyces_uvarum', 'Saccharomyces kudriavzevii': 'Saccharomyces_kudriavzevii', 'Saccharomyces mikatae': 'Saccharomyces_mikatae', 'Saccharomyces cerevisiae': 'Saccharomyces_cerevisiae'}
YGOB_PRE = {'Zygosaccharomyces rouxii': 'Zygosaccharomyces_rouxii', 'Torulaspora delbrueckii': 'Torulaspora_delbrueckii', 'Kluyveromyces lactis': 'Kluyveromyces_lactis', 'Eremothecium gossypii': 'Eremothecium_gossypii', 'Eremothecium cymbalariae': 'Eremothecium_cymbalariae', 'Lachancea kluyveri': 'Lachancea_kluyveri', 'Lachancea thermotolerans': 'Lachancea_thermotolerans', 'Lachancea waltii': 'Lachancea_waltii'}
DIAS_SYN = {'Candida glabrata': 'Nakaseomyces_glabratus', 'Kluyveromyces marxianus var. marxianus': 'Kluyveromyces_marxianus', 'Lachancea kluyvery': 'Lachancea_kluyveri', 'Hansenula polymorpha': 'Ogataea_polymorpha', 'Dekkera bruxellensis': 'Brettanomyces_bruxellensis', 'Saccharomyces bayanus (HYBRID)': None, 'Zygosaccharomyces bailii (HYBRID)': None}
ANCHORS = {'Saccharomyces_cerevisiae': 'S288C (reference; the eight query paralogs)', 'Candida_albicans': 'SC5314 (Cdr1; H0c Shin 2026 anchor)', 'Candida_auris': 'CUG-Ser1 pathogen (Cdr1 orthologs)', 'Nakaseomyces_glabratus': 'C. glabrata (Cdr1 structures 9S2C etc.)', 'Komagataella_phaffii': 'Pichiales model', 'Yarrowia_lipolytica': 'Dias 2025: lacks Snq2/Pdr18', 'Lipomyces_starkeyi': 'Dias 2025: lacks Snq2/Pdr18', 'Nadsonia_fulvescens': 'Dias 2025: lacks Snq2/Pdr18', 'Hanseniaspora_valbyensis': 'Dias 2025: lacks Snq2/Pdr18', 'Pachysolen_tannophilus': 'CUG-Ala reference (NCBI table 26)', 'Ascoidea_rubescens': 'CUG-Ser2 / Dias 2025'}
ANCHOR_STRAIN = {'Saccharomyces_cerevisiae': 'Saccharomyces_cerevisiae', 'Candida_albicans': 'Candida_albicans_SC5314'}
ap = argparse.ArgumentParser()
ap.add_argument('--names', required=True)
ap.add_argument('--itol', required=True)
ap.add_argument('--tree', required=True)
ap.add_argument('--og-proxy', default='', help='OG0000012_per_species_counts.tsv (proxy for sequences per genome before retrieval is complete)')
ap.add_argument('--proxy-factor', type=float, default=1.0, help='multiplier applied to OG counts (PASS/OG ratio, e.g. 0.70 from the test prefix)')
ap.add_argument('--retrieval', default='', help='retrieval_table.tsv from the full run (PASS and NR-representative rows are counted)')
ap.add_argument('--dias', default='')
ap.add_argument('--outdir', required=True)
ap.add_argument('--seed', type=int, default=SEED)
ap.add_argument('--target-seqs', type=int, default=1250)
ap.add_argument('--min-seqs', type=int, default=1000)
ap.add_argument('--max-seqs', type=int, default=1500)
ap.add_argument('--min-per-order', type=int, default=3)
ap.add_argument('--wgd-mode', choices=['all', 'capped'], default='all')
ap.add_argument('--wgd-cap', type=int, default=8)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

def H(code):
    return int(hashlib.sha256(f'{a.seed}|{code}'.encode()).hexdigest(), 16)
tip2order = {}
lines = [l.rstrip('\n').split('\t') for l in open(a.itol, encoding='utf-8') if l.strip() and (not l.startswith('#'))]
cols = next((l for l in lines if l[0] == 'LEGEND_COLORS'))[1:]
labs = next((l for l in lines if l[0] == 'LEGEND_LABELS'))[1:]
c2o = {c.lower(): o for c, o in zip(cols, labs)}
for l in lines:
    if len(l) == 3 and l[1].startswith('#'):
        tip2order[l[0]] = c2o.get(l[1].lower(), '?')
genomes = []
for r in read_tsv(a.names):
    ft = r['FINAL_Tip_ID']
    k, strain = species_key_of(ft)
    genomes.append(dict(speciesCode=sanitize(ft), FINAL_Tip_ID=ft, species_key=k, strain=strain, genus=ft.split('_')[0], order=tip2order.get(ft, '?')))
by_sp = collections.defaultdict(list)
for g in genomes:
    by_sp[g['species_key']].append(g)
counts, count_src = ({}, '')
if a.retrieval:
    for r in read_tsv(a.retrieval):
        if r['retrieval_class'] == 'PASS' and r['genome_kind'] == 'y1000_genome' and (r['passes_rule_final_seq'] == 'Y') and (r['cdhit99_rep_global'] == 'Y'):
            counts[r['species_code']] = counts.get(r['species_code'], 0) + 1
    count_src = 'retrieved PASS sequences (CD-HIT 0.99 representatives) from ' + os.path.basename(a.retrieval)
elif a.og_proxy:
    for r in read_tsv(a.og_proxy):
        if r['column_kind'] == 'yeast_Y1000':
            counts[sanitize(r['FINAL_Tip_ID'])] = float(r['n_genes_in_OG']) * a.proxy_factor
    count_src = f'PROXY: OG0000012 genes per genome x {a.proxy_factor}'
for g in genomes:
    g['n_seq'] = counts.get(g['speciesCode'], 0)
rep = {}
for k, gs in by_sp.items():
    anc = ANCHOR_STRAIN.get(k)
    cand = [g for g in gs if g['speciesCode'] == anc] if anc else []
    rep[k] = cand[0] if cand else min(gs, key=lambda g: H(g['speciesCode']))
n_multi = sum((1 for gs in by_sp.values() if len(gs) > 1))
from ete3 import Tree
tree = Tree(a.tree, format=1)
leaf = {l.name: l for l in tree.iter_leaves()}
depth = {}
for n in tree.traverse('preorder'):
    depth[n] = 0.0 if n.up is None else depth[n.up] + (n.dist or 0.0)
anc_cache = {}

def ancestors(node):
    if node not in anc_cache:
        s, x = ([], node)
        while x is not None:
            s.append(x)
            x = x.up
        anc_cache[node] = s
    return anc_cache[node]

def pdist(a_, b_):
    A = ancestors(a_)
    sb = set(ancestors(b_))
    for x in A:
        if x in sb:
            return depth[a_] + depth[b_] - 2 * depth[x]
    return float('inf')
layer = collections.defaultdict(set)
post_tips = []
ygob_rows = []
for name, k in list(YGOB_POST.items()) + list(YGOB_PRE.items()):
    ok = k in rep
    ygob_rows.append(dict(ygob_species=name, group='post-WGD' if name in YGOB_POST else 'pre-WGD', y1000_species_key=k if ok else '', matched='Y' if ok else 'N', y1000_genome=rep[k]['speciesCode'] if ok else ''))
    if ok:
        layer[k].add('YGOB')
        if name in YGOB_POST:
            post_tips.append(leaf[rep[k]['FINAL_Tip_ID']])
mrca = tree.get_common_ancestor(post_tips)
wgd_species = sorted({species_key_of(l.name)[0] for l in mrca.iter_leaves()})
wgd_genomes_all = sum((1 for l in mrca.iter_leaves()))
dias_rows = []
if a.dias:
    seen = set()
    for r in read_tsv(a.dias):
        sp = r['species'].strip()
        if sp in seen:
            continue
        seen.add(sp)
        key = DIAS_SYN[sp] if sp in DIAS_SYN else sp.replace(' ', '_')
        ok = key is not None and key in rep
        dias_rows.append(dict(dias_species=sp, y1000_species_key=key or '', matched='Y' if ok else 'N', note='synonym table' if sp in DIAS_SYN and key else 'hybrid, not matched' if sp in DIAS_SYN else '', y1000_genome=rep[key]['speciesCode'] if ok else ''))
        if ok:
            layer[key].add('DIAS')
for k, why in ANCHORS.items():
    if k in rep:
        layer[k].add('ANCHOR')
write_tsv(os.path.join(a.outdir, 'ygob_match.tsv'), ygob_rows)
write_tsv(os.path.join(a.outdir, 'dias_match.tsv'), dias_rows)
protected = set(layer)
if a.wgd_mode == 'all':
    for k in wgd_species:
        layer[k].add('WGD')
else:
    by_genus = collections.defaultdict(list)
    for k in wgd_species:
        by_genus[k.split('_')[0]].append(k)
    for gen, ks in by_genus.items():
        keep = [k for k in ks if k in protected]
        rest = sorted([k for k in ks if k not in protected], key=lambda k: H(rep[k]['speciesCode']))
        chosen = list(keep)
        pool = list(rest)
        while len(chosen) < a.wgd_cap and pool:
            if not chosen:
                nxt = pool[0]
            else:
                nxt = max(pool, key=lambda k: (min((pdist(leaf[rep[k]['FINAL_Tip_ID']], leaf[rep[c]['FINAL_Tip_ID']]) for c in chosen)), -H(rep[k]['speciesCode'])))
            chosen.append(nxt)
            pool.remove(nxt)
        for k in chosen:
            layer[k].add('WGD')
mand = set(layer)
mand_genomes = {k: rep[k] for k in mand}
S_mand = sum((g['n_seq'] for g in mand_genomes.values()))
orders = sorted({g['order'] for g in genomes})
elig = {o: sorted([k for k in by_sp if rep[k]['order'] == o and k not in mand], key=lambda k: H(rep[k]['speciesCode'])) for o in orders}
mean_seq = {}
allg = [g for g in genomes if g['n_seq'] > 0]
glob_mean = sum((g['n_seq'] for g in allg)) / max(1, len(allg))
for o in orders:
    vs = [rep[k]['n_seq'] for k in by_sp if rep[k]['order'] == o and rep[k]['n_seq'] > 0]
    mean_seq[o] = sum(vs) / len(vs) if vs else glob_mean

def quotas(c):
    q = {}
    for o in orders:
        N = len(elig[o])
        if N == 0:
            q[o] = 0
        elif N <= a.min_per_order:
            q[o] = N
        else:
            q[o] = min(N, max(a.min_per_order, int(math.floor(c * math.sqrt(N)))))
    return q

def expected(c):
    q = quotas(c)
    return S_mand + sum((q[o] * mean_seq[o] for o in orders))
lo, hi = (0.0, 50.0)
for _ in range(60):
    mid = (lo + hi) / 2
    if expected(mid) <= a.target_seqs:
        lo = mid
    else:
        hi = mid
c_star = lo
quota = quotas(c_star)
chosen = {k: set(v) for k, v in layer.items()}
for o in orders:
    sel = [k for k in mand if rep[k]['order'] == o]
    pool = list(elig[o])
    cur = list(sel)
    take = quota[o]
    picked = []
    while len(picked) < take and pool:
        if not cur:
            nxt = pool[0]
        else:
            curtips = [leaf[rep[c]['FINAL_Tip_ID']] for c in cur]
            nxt = max(pool, key=lambda k: (min((pdist(leaf[rep[k]['FINAL_Tip_ID']], t) for t in curtips)), -H(rep[k]['speciesCode'])))
        picked.append(nxt)
        cur.append(nxt)
        pool.remove(nxt)
    for k in picked:
        chosen.setdefault(k, set()).add('STRATIFIED')
rows = []
for g in genomes:
    k = g['species_key']
    is_rep = rep[k]['speciesCode'] == g['speciesCode']
    ch = k in chosen and is_rep
    rows.append(dict(speciesCode=g['speciesCode'], FINAL_Tip_ID=g['FINAL_Tip_ID'], species_key=k, order=g['order'], n_genomes_of_species=len(by_sp[k]), representative_of_species='Y' if is_rep else 'N', layer='+'.join(sorted(chosen[k])) if ch else '', selected='Y' if ch else 'N', n_sequences='%.1f' % g['n_seq'] if isinstance(g['n_seq'], float) else g['n_seq'], hash=H(g['speciesCode']) % 10 ** 8))
write_tsv(os.path.join(a.outdir, 'sampling_plan.tsv'), rows)
sel_rows = [r for r in rows if r['selected'] == 'Y']
tot = sum((float(r['n_sequences'] or 0) for r in sel_rows))
by_order = {}
for o in orders:
    allo = [r for r in rows if r['order'] == o and r['representative_of_species'] == 'Y']
    s = [r for r in sel_rows if r['order'] == o]
    by_order[o] = dict(genomes=sum((1 for g in genomes if g['order'] == o)), species=len(allo), selected=len(s), selected_mandatory=sum((1 for r in s if r['layer'] != 'STRATIFIED')), selected_stratified=sum((1 for r in s if r['layer'] == 'STRATIFIED')), sequences=round(sum((float(r['n_sequences'] or 0) for r in s)), 1))
summ = dict(seed=a.seed, count_source=count_src, wgd_mode=a.wgd_mode, n_genomes=len(genomes), n_species_keys=len(by_sp), n_species_with_more_than_one_genome=n_multi, wgd_clade_genomes=wgd_genomes_all, wgd_clade_species=len(wgd_species), wgd_genera=dict(collections.Counter((k.split('_')[0] for k in wgd_species))), ygob_matched=f"{sum((1 for r in ygob_rows if r['matched'] == 'Y'))}/20", dias_species_in_table=len(dias_rows), dias_matched=sum((1 for r in dias_rows if r['matched'] == 'Y')), dias_unmatched=[r['dias_species'] for r in dias_rows if r['matched'] == 'N'], mandatory_species=len(mand), mandatory_by_layer=dict(collections.Counter((l for v in layer.values() for l in v))), mandatory_sequences=round(S_mand, 1), c_star=round(c_star, 4), quotas=quota, selected_genomes=len(sel_rows), selected_sequences=round(tot, 1), target=[a.min_seqs, a.target_seqs, a.max_seqs], in_target_band=bool(a.min_seqs <= tot <= a.max_seqs), by_order=by_order)
json.dump(summ, open(os.path.join(a.outdir, 'sampling_summary.json'), 'w'), indent=1)
log(json.dumps({k: v for k, v in summ.items() if k not in ('by_order', 'quotas', 'wgd_genera')}))
for o, v in by_order.items():
    log(o, v)

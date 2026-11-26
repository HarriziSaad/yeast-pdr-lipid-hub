import argparse, collections, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_tsv, write_tsv, log
from treeutil import Tree
ap = argparse.ArgumentParser()
ap.add_argument('--tree', required=True)
ap.add_argument('--manifest', default='/home/pdr/work/phase1b/set/analysis_set_manifest.tsv')
ap.add_argument('--ygob-map', default='/home/pdr/work/phase1b/set/ygob/ygob_map.tsv')
ap.add_argument('--outdir', required=True)
ap.add_argument('--sh-min', type=float, default=80.0)
ap.add_argument('--ufb-min', type=float, default=95.0)
ap.add_argument('--viol-max', type=float, default=0.05)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
S288C = {'PDR5': 'Saccharomyces_cerevisiae|g002312.m1', 'PDR15': 'Saccharomyces_cerevisiae|g004097.m1', 'PDR10': 'Saccharomyces_cerevisiae|g002467.m1', 'SNQ2': 'Saccharomyces_cerevisiae|g003736.m1', 'PDR18': 'Saccharomyces_cerevisiae|g002018.m1', 'PDR12': 'Saccharomyces_cerevisiae|g002739.m1', 'PDR11': 'Saccharomyces_cerevisiae|g005516.m1', 'AUS1': 'Saccharomyces_cerevisiae|g002189.m1'}
CRITICAL = ['PDR5', 'PDR15', 'SNQ2', 'PDR11', 'AUS1']
ALL8 = ['PDR5', 'PDR15', 'PDR10', 'SNQ2', 'PDR18', 'PDR12', 'PDR11', 'AUS1']
WGD_UNIVERSE = {'PDR5', 'PDR15', 'PDR11', 'AUS1'}
FAMILY = {'PDR5': 'PDR5-like', 'PDR15': 'PDR5-like', 'PDR10': 'PDR5-like', 'SNQ2': 'SNQ2-like', 'PDR18': 'SNQ2-like', 'PDR12': 'PDR12-like', 'PDR11': 'PDR11/AUS1-like', 'AUS1': 'PDR11/AUS1-like'}
EXPECT_NEAREST = {'PDR5': {'PDR15', 'PDR10'}, 'PDR15': {'PDR5', 'PDR10'}, 'PDR10': {'PDR5', 'PDR15'}, 'SNQ2': {'PDR18'}, 'PDR18': {'SNQ2'}, 'PDR11': {'AUS1'}, 'AUS1': {'PDR11'}, 'PDR12': None}
man = {r['seq_id']: r for r in read_tsv(a.manifest)}
tree = Tree(open(a.tree).read())
name2node = {tree.name[n]: n for n in tree.leaves()}
assert set(name2node) == set(man), len(set(name2node) ^ set(man))
has_ufb = any((e['ufb'] is not None for e in tree.edge.values()))
log('leaves', len(name2node), 'edges', len(tree.edge), 'UFBoot values present:', has_ufb)
og = {name2node[i] for i, r in man.items() if r['role'] == 'OUT_OG658'}
tmp = tree.root_on_edge(next(iter(tree.edge)))
total = len(name2node)
best = None
for n in tmp.order:
    ls = tmp.leafnodes(n)
    k = len(ls & og)
    for side_n, side_k in ((len(ls), k), (total - len(ls), len(og) - k)):
        f1 = 2 * side_k / (len(og) + side_n)
        e = tree.edge[tmp.stem[n]]
        cand = (f1, e['sh'] if e['sh'] is not None else -1, -side_n, tmp.stem[n], side_n, side_k)
        if best is None or cand[:3] > best[:3]:
            best = cand
f1, _, _, root_eid, side_n, side_k = best
R = tree.root_on_edge(root_eid)
root_info = dict(outgroup_n=len(og), outgroup_side_size=side_n, outgroup_in_side=side_k, F1=round(f1, 4), outgroup_monophyletic=f1 == 1.0, root_edge_sh=tree.edge[root_eid]['sh'], root_edge_ufb=tree.edge[root_eid]['ufb'])
sideA, sideB = R.top
root_og_side = sideA if len(R.leafnodes(sideA) & og) >= len(R.leafnodes(sideB) & og) else sideB
dik = {name2node[i] for i, r in man.items() if r['role'] == 'OUT_DIKARYA'}
root_info['dikarya_on_outgroup_side'] = len(R.leafnodes(root_og_side) & dik)
root_info['dikarya_total'] = len(dik)
log('root', root_info)
sp = {n: man[tree.name[n]]['species_code'] for n in tree.leaves()}
order_of = {n: man[tree.name[n]]['order'] for n in tree.leaves()}
role = {n: man[tree.name[n]]['role'] for n in tree.leaves()}
wgd = {n: 'WGD' in man[tree.name[n]]['layer'].split('+') and role[n] == 'INGROUP' for n in tree.leaves()}
best_hit = {n: man[tree.name[n]]['best_hit'] for n in tree.leaves()}
fam = {n: man[tree.name[n]]['family'] for n in tree.leaves()}
s288 = {x: name2node[i] for x, i in S288C.items()}
s288_set = set(s288.values())

def sup_ok(n):
    sh, ufb = R.support(n)
    sh_ok = sh is not None and sh >= a.sh_min
    ufb_ok = True if not has_ufb else ufb is not None and ufb >= a.ufb_min
    return (sh, ufb, sh_ok and ufb_ok)
clade = {}
for x in ALL8:
    n = s288[x]
    best_n = n
    cur = n
    while R.parent[cur] != 'ROOT':
        p = R.parent[cur]
        if len(R.leafnodes(p) & s288_set) > 1:
            break
        cur = p
        best_n = cur
    clade[x] = best_n
rows = []
viol_rows = []
res = {}
for x in ALL8:
    c = clade[x]
    ls = R.leafnodes(c)
    species = {sp[n] for n in ls}
    sh, ufb, ok_sup = sup_ok(c)
    univ = [n for n in tree.leaves() if role[n] == 'INGROUP' and (wgd[n] if x in WGD_UNIVERSE else True)]
    A = [n for n in univ if best_hit[n] == x]
    viol = [n for n in A if n not in ls]
    rate = len(viol) / len(A) if A else float('nan')
    in_fam = sum((1 for n in ls if role[n] == 'INGROUP' and fam[n] == FAMILY[x]))
    n_ing = sum((1 for n in ls if role[n] == 'INGROUP'))
    size_ok = len(ls) >= 10 and len(species) >= 5
    t_ok = ok_sup and len(ls) >= 2
    l_ok = t_ok and (len(A) == 0 or rate <= a.viol_max)
    where = collections.Counter()
    for n in viol:
        w = 'elsewhere'
        for y in ALL8:
            if y != x and n in R.leafnodes(clade[y]):
                w = 'in_C_' + y
                break
        where[w] += 1
    res[x] = dict(clade=c, n_seq=len(ls), n_species=len(species), n_ingroup=n_ing, n_wgd_species=len({sp[n] for n in ls if wgd[n]}), sh=sh, ufb=ufb, support_ok=ok_sup, scoring_group_3f=size_ok, assigned_best_hit=len(A), violators=len(viol), violator_rate=rate, pass_1T=t_ok, pass_1L=l_ok)
    rows.append(dict(clade=x, family=FAMILY[x], n_sequences=len(ls), n_ingroup=n_ing, n_species=len(species), n_wgd_species=res[x]['n_wgd_species'], stem_SH_aLRT=sh, stem_UFBoot=ufb, support_ok='Y' if ok_sup else 'N', scoring_group_3f='Y' if size_ok else 'N', family_purity=round(in_fam / n_ing, 3) if n_ing else '', assigned_best_hit_in_universe=len(A), violators_outside_clade=len(viol), violator_rate=round(rate, 4) if A else '', violators_where=json.dumps(dict(where)), pass_1T='Y' if t_ok else 'N', pass_1L='Y' if l_ok else 'N', orders=json.dumps(dict(collections.Counter((order_of[n] for n in ls if role[n] == 'INGROUP'))))))
    for n in viol:
        viol_rows.append(dict(clade=x, seq_id=tree.name[n], species=sp[n], best_hit=best_hit[n]))
write_tsv(os.path.join(a.outdir, 'g1_clades.tsv'), rows)
write_tsv(os.path.join(a.outdir, 'g1_literal_violators.tsv'), viol_rows)
pairs = {}
for x, y in (('PDR5', 'PDR15'), ('PDR11', 'AUS1')):
    cx, cy = (clade[x], clade[y])
    px, py = (R.parent[cx], R.parent[cy])
    sister = px == py and px != 'ROOT' and (set(R.children[px]) == {cx, cy})
    if sister:
        pn = px
    else:
        pn = R.mrca([s288[x], s288[y]])
    sh, ufb, ok_sup = sup_ok(pn) if pn is not None else (None, None, False)
    pls = R.leafnodes(pn) if pn is not None else set()
    non_sacch = sorted({order_of[n] for n in pls if role[n] == 'INGROUP' and order_of[n] != 'Saccharomycetales'})
    nested = len(non_sacch) == 0 and len(pls) > 0
    pairs[x, y] = dict(pair_node=pn, sister=sister, sh=sh, ufb=ufb, support_ok=ok_sup, n_seq=len(pls), n_species=len({sp[n] for n in pls}), orders_non_sacch=non_sacch, nested_in_Saccharomycetales=nested, n_wgd=sum((1 for n in pls if wgd[n])), n_nonwgd_ingroup=sum((1 for n in pls if role[n] == 'INGROUP' and (not wgd[n]))), other_s288c_in_pair=sorted((k for k, v in s288.items() if k not in (x, y) and v in pls)))
yg = read_tsv(a.ygob_map)
ycheck = []
by_genome = collections.defaultdict(list)
for r in yg:
    by_genome[r['pillar'], r['y1000_genome']].append(r)
for r in yg:
    pil, g = (r['pillar'], r['y1000_genome'])
    row = dict(pillar=pil, family=r['family'], ygob_species=r['ygob_species'], group=r['group'], position=r['position'], ygob_gene=r['ygob_gene'], mapped_seq=r['mapped_seq'], status=r['status'])
    if not r['mapped_seq']:
        row.update(expected='NA', observed='unmapped', ok='NA')
        ycheck.append(row)
        continue
    n = name2node[r['mapped_seq']]
    where = [x for x in ALL8 if n in R.leafnodes(clade[x])]
    row['observed'] = ','.join(('C_' + w for w in where)) if where else 'none'
    if r['family'] == 'PDR5/PDR15':
        pn = pairs['PDR5', 'PDR15']['pair_node']
        inpair = pn is not None and n in R.leafnodes(pn)
        if r['group'] == 'post-WGD':
            row['expected'] = 'inside pair clade (C_PDR5 or C_PDR15)'
            row['ok'] = 'Y' if inpair else 'N'
        else:
            row['expected'] = 'outside pair clade'
            row['ok'] = 'Y' if not inpair else 'N'
    elif r['family'] == 'PDR11/AUS1':
        pn = pairs['PDR11', 'AUS1']['pair_node']
        inpair = pn is not None and n in R.leafnodes(pn)
        if r['group'] == 'post-WGD':
            row['expected'] = 'inside pair clade (C_PDR11 or C_AUS1)'
            row['ok'] = 'Y' if inpair else 'N'
        else:
            row['expected'] = 'outside pair clade'
            row['ok'] = 'Y' if not inpair else 'N'
    elif r['family'] == 'SNQ2':
        row['expected'] = 'C_SNQ2'
        row['ok'] = 'Y' if n in R.leafnodes(clade['SNQ2']) else 'N'
    else:
        row['expected'] = 'C_PDR12'
        row['ok'] = 'Y' if n in R.leafnodes(clade['PDR12']) else 'N'
    ycheck.append(row)
two = []
for (pil, g), rs in by_genome.items():
    if pil not in ('Anc_5.500', 'Anc_7.123a'):
        continue
    rs = [r for r in rs if r['mapped_seq'] and r['group'] == 'post-WGD']
    if len(rs) >= 2:
        xs, ys = ('PDR5', 'PDR15') if pil == 'Anc_5.500' else ('PDR11', 'AUS1')
        obs = []
        for r in rs:
            n = name2node[r['mapped_seq']]
            obs.append('C_' + xs if n in R.leafnodes(clade[xs]) else 'C_' + ys if n in R.leafnodes(clade[ys]) else 'other')
        ok = sorted(obs) == sorted(['C_' + xs, 'C_' + ys]) if len(rs) == 2 else len(set(obs)) >= 2
        two.append(dict(pillar=pil, genome=g, members=';'.join((r['mapped_seq'] for r in rs)), observed=';'.join(obs), ok='Y' if ok else 'N'))
write_tsv(os.path.join(a.outdir, 'g1_ygob_check.tsv'), ycheck)
write_tsv(os.path.join(a.outdir, 'g1_ygob_two_copy_genomes.tsv'), two)
ycnt = collections.Counter((r['ok'] for r in ycheck))
log('YGOB check', dict(ycnt), 'two-copy genomes', collections.Counter((r['ok'] for r in two)))

def patristic(x, y):

    def up(n):
        out, d, cur = ([], 0.0, n)
        while True:
            out.append((cur, d))
            if R.parent[cur] == 'ROOT':
                return (out, d)
            d += tree.edge[R.stem[cur]]['length'] or 0.0
            cur = R.parent[cur]
    px, tx = up(s288[x])
    py, ty = up(s288[y])
    dy = dict(py)
    for n, d1 in px:
        if n in dy:
            return d1 + dy[n]
    return tx + ty + (tree.edge[R.root_edge]['length'] or 0.0)
cons = []
for x in ALL8:
    others = [(patristic(x, y), y) for y in ALL8 if y != x]
    others.sort()
    nearest = others[0][1]
    exp = EXPECT_NEAREST[x]
    cons.append(dict(paralog=x, nearest_other_S288C=nearest, nearest_distance=round(others[0][0], 4), second=others[1][1], expected=','.join(sorted(exp)) if exp else 'none (singleton family)', consistent='Y' if exp is None or nearest in exp else 'N'))
write_tsv(os.path.join(a.outdir, 'g1_s288c_family_consistency.tsv'), cons)
c1_T = all((res[x]['pass_1T'] for x in CRITICAL))
c1_L = all((res[x]['pass_1L'] for x in CRITICAL))
pair_ok = all((pairs[k]['sister'] and pairs[k]['support_ok'] and pairs[k]['nested_in_Saccharomycetales'] for k in pairs))
ygob_ok = all((r['ok'] in ('Y', 'NA') for r in ycheck)) and all((r['ok'] == 'Y' for r in two))
c2 = pair_ok and ygob_ok
n_T = sum((1 for x in ALL8 if res[x]['pass_1T']))
n_L = sum((1 for x in ALL8 if res[x]['pass_1L']))
c3_T, c3_L = (n_T >= 6, n_L >= 6)
summary = dict(tree=a.tree, ufb_present=has_ufb, thresholds=dict(sh=a.sh_min, ufb=a.ufb_min, violators_max=a.viol_max), rooting=root_info, criterion1_T_all_five=c1_T, criterion1_L_all_five=c1_L, criterion2_pairs=pair_ok, criterion2_ygob=ygob_ok, criterion2=c2, n_clades_pass_1T=n_T, n_clades_pass_1L=n_L, criterion3_T=c3_T, criterion3_L=c3_L, G1_PASS_T=bool(c1_T and c2 and c3_T), G1_PASS_L=bool(c1_L and c2 and c3_L), pairs={f'{k[0]}/{k[1]}': {kk: vv for kk, vv in v.items() if kk != 'pair_node'} for k, v in pairs.items()}, ygob_counts=dict(ycnt), ygob_failures=[(r['pillar'], r['ygob_gene'], r['mapped_seq'], r['observed']) for r in ycheck if r['ok'] == 'N'], two_copy_failures=[r for r in two if r['ok'] == 'N'], clades={x: {k: v for k, v in res[x].items() if k != 'clade'} for x in ALL8}, s288c_consistency={r['paralog']: (r['nearest_other_S288C'], r['consistent']) for r in cons})
json.dump(summary, open(os.path.join(a.outdir, 'g1_summary.json'), 'w'), indent=1, default=str)
print(json.dumps(summary, indent=1, default=str))
inv = {v: k for k, v in s288.items()}
res_bb = {}
for n in reversed(R.order):
    if n in tree.name:
        res_bb[n] = inv.get(n)
        continue
    kids = [res_bb[k] for k in R.children[n] if res_bb[k] is not None]
    if not kids:
        res_bb[n] = None
    elif len(kids) == 1:
        res_bb[n] = kids[0]
    else:
        sh_, ufb_ = R.support(n)
        res_bb[n] = '(' + ','.join(kids) + ')' + ('%s' % sh_ if sh_ is not None else '') + ('/%s' % ufb_ if ufb_ is not None else '')
bb_top = [res_bb[t] for t in R.top if res_bb[t] is not None]
open(os.path.join(a.outdir, 'g1_backbone.txt'), 'w').write('(' + ','.join(bb_top) + ');' + chr(10))
split = collections.Counter()
for r in ycheck:
    split['ohnolog pillars' if r['pillar'] in ('Anc_5.500', 'Anc_7.123a') else 'non-ohnolog pillars (SNQ2, PDR12)', r['ok']] += 1
sibl = {}
for (x, y), v in pairs.items():
    pn = v['pair_node']
    if pn is not None and R.parent[pn] != 'ROOT':
        sib = [k for k in R.children[R.parent[pn]] if k != pn]
        sibl['%s/%s' % (x, y)] = dict(sibling_n_seq=[len(R.leafnodes(k)) for k in sib], sibling_has_s288c=[sorted((inv[m] for m in R.leafnodes(k) if m in inv)) for k in sib], sibling_orders=[dict(collections.Counter((order_of[m] for m in R.leafnodes(k)))) for k in sib])
summary['ygob_split_by_pillar_type'] = {'%s|%s' % k: v for k, v in split.items()}
summary['pair_sibling_context'] = sibl
json.dump(summary, open(os.path.join(a.outdir, 'g1_summary.json'), 'w'), indent=1, default=str)
lab = []
for n in tree.leaves():
    ls_ = [x for x in ALL8 if n in R.leafnodes(clade[x])]
    lab.append(dict(seq_id=tree.name[n], role=role[n], species=sp[n], order=order_of[n], best_hit=best_hit[n], tree_clade=ls_[0] if ls_ else 'unassigned'))
write_tsv(os.path.join(a.outdir, 'g1_tree_clade_labels.tsv'), lab)

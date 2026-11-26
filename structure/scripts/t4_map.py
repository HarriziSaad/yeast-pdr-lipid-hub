import sys, os, json, csv, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
import p3lib as L
info = json.load(open(os.path.join(L.OUT, 'structure_info.json')))
classes = L.load_tsv(os.path.join(L.IN, 'site_classes_FROZEN.tsv'))
AF = {'Pdr5': 'P33302', 'Snq2': 'P32568', 'Pdr10': 'P51533', 'Pdr11': 'P40550', 'Pdr12': 'Q02785', 'Pdr15': 'Q04182', 'Pdr18': 'P53756', 'Aus1': 'Q08409'}
CONTEXT = {'YOL075C': 'Q08234', 'Adp1': 'P25371', 'CaCdr1_AF': 'Q5ANA3'}
EXP = ['22FR', '9VT1', '9VTO', '9IUK']
DOMS = ['D_NBD1', 'D_TMD1', 'D_NBD2', 'D_TMD2']
nums7p04 = set(info['7P04']['nums'])
list7p04 = sorted(nums7p04)
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score, aligner.extend_gap_score = (-11, -1)
aligner.mode = 'global'
p5seq = info['Pdr5']['seq']

def seqmap(modelseq, modelnums):
    aln = aligner.align(p5seq, modelseq)[0]
    out = {}
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for k in range(a1 - a0):
            out[a0 + k + 1] = modelnums[b0 + k]
    return out

def load(m, ref):
    fp = os.path.join(L.OUT, 'align', f'{m}__to__{ref}.json')
    return json.load(open(fp)) if os.path.exists(fp) else None

def load_dom(m):
    d = {}
    for dn in DOMS:
        fp = os.path.join(L.OUT, 'align_dom', f'{m}__{dn}.json')
        if not os.path.exists(fp):
            return None
        for p in json.load(open(fp))['pairs']:
            d[p[1]] = (p[0], p[2])
    return d

def plddt_of(name):
    st, _ = L.read_clean(os.path.join(L.OUT, 'pdb', f'{name}.pdb'))
    return {r['num']: r['b'] for r in L.residues(st) if r['ca'] is not None}

def mode4(primary, dom, best, sq):
    votes = collections.Counter((x for x in (primary, dom, best, sq) if x not in ('', None)))
    if not votes:
        return ''
    top = max(votes.values())
    for cand in (dom, primary, best, sq):
        if cand not in ('', None) and votes[cand] == top:
            return cand

def map_one(name, evidence, is_af):
    nums = info[name]['nums']
    seq = info[name]['seq']
    pl = plddt_of(name) if is_af else {}
    a4 = load(name, '7P04')
    tms = {ref: load(name, ref)['TM_norm_ref'] for ref in ('7P04', '7P03', '7P06') if load(name, ref)}
    best = max(tms, key=tms.get)
    ab = load(name, best)
    d4 = {p[1]: (p[0], p[2]) for p in a4['pairs']}
    db = {p[1]: (p[0], p[2]) for p in ab['pairs']}
    dd = load_dom(name) or {}
    sm = seqmap(seq, nums)
    aligned_ref = sorted(d4)
    first, last = (aligned_ref[0], aligned_ref[-1])
    out = []
    for c in classes:
        rn = int(c['residue'])
        flags = []
        pr, d, aa, plddt = ('', '', '', '')
        cons, nag = ('', '')
        r_best = db[rn][0] if rn in db else ''
        r_dom = dd[rn][0] if rn in dd else ''
        r_seq = sm.get(rn, '')
        note = ''
        if name == 'Pdr18' and rn >= 1473:
            raw = d4.get(rn)
            note = 'S288C Pdr18 lacks TMS12; raw 7P04 pairing ignored' + (' (raw partner %d, CA %.2f A)' % (raw[0], raw[1]) if raw else '')
            flags = ['NOT_ALIGNABLE', 'PDR18_NO_TMS12']
            out.append({'evidence': evidence, 'class': c['class'], 'P33302_residue': rn, 'P33302_aa': c['aa'], 'paralog': name, 'paralog_residue': '', 'aa': '', 'CA_distance_A': '', 'pLDDT': '', 'flags': ';'.join(flags), 'sidechain_ok': 'N', 'n_agree_of_3': '', 'consensus_residue': '', 'best_ref': best, 'residue_bestref': '', 'CA_dist_bestref_A': '', 'residue_domainwise': '', 'CA_dist_domainwise_A': '', 'residue_seqaln': '', 'TM_7P04': '%.3f' % a4['TM_norm_ref'], 'RMSD_7P04_A': '%.2f' % a4['RMSD'], 'note': note})
            continue
        if rn not in d4:
            idx = list7p04.index(rn)
            lo = idx
            while lo - 1 >= 0 and list7p04[lo - 1] not in d4:
                lo -= 1
            hi = idx
            while hi + 1 < len(list7p04) and list7p04[hi + 1] not in d4:
                hi += 1
            if rn < first or rn > last or hi - lo + 1 >= 5:
                flags.append('NOT_ALIGNABLE')
            else:
                flags.append('NO_ALIGNED_RESIDUE')
            cons = mode4('', r_dom, r_best, r_seq) if 'NOT_ALIGNABLE' not in flags else ''
        else:
            pr, d = d4[rn]
            aa = seq[nums.index(pr)]
            if is_af:
                plddt = pl[pr]
                if plddt < 70:
                    flags.append('PLDDT_LT70')
                elif plddt < 80:
                    flags.append('PLDDT_70_80')
            if d > 3.0:
                flags.append('MISALIGNED_GT3A')
            nag = sum((1 for x in (r_best, r_dom, r_seq) if x == pr))
            cons = mode4(pr, r_dom, r_best, r_seq)
            if nag == 0:
                flags.append('MAPPING_UNCERTAIN')
            if cons != pr:
                flags.append('MAPPING_DIFFERS')
        d_dom = dd[rn][1] if rn in dd else None
        if pr != '' and d_dom is not None and (d_dom <= 3.0) and (cons == pr) and (nag >= 2):
            if is_af:
                sc_ok = 'Y' if plddt >= 80 else 'N'
            else:
                sc_ok = 'NA(E)'
        else:
            sc_ok = 'N'
        out.append({'evidence': evidence, 'class': c['class'], 'P33302_residue': rn, 'P33302_aa': c['aa'], 'paralog': name, 'paralog_residue': pr, 'aa': aa, 'CA_distance_A': '' if d == '' else '%.2f' % d, 'pLDDT': '' if plddt == '' else '%.1f' % plddt, 'flags': ';'.join(flags) if flags else 'OK', 'sidechain_ok': sc_ok, 'n_agree_of_3': nag, 'consensus_residue': cons, 'best_ref': best, 'residue_bestref': r_best, 'CA_dist_bestref_A': '%.2f' % db[rn][1] if rn in db else '', 'residue_domainwise': r_dom, 'CA_dist_domainwise_A': '%.2f' % dd[rn][1] if rn in dd else '', 'residue_seqaln': r_seq, 'TM_7P04': '%.3f' % a4['TM_norm_ref'], 'RMSD_7P04_A': '%.2f' % a4['RMSD'], 'note': ''})
    return out
cols = ['evidence', 'class', 'P33302_residue', 'P33302_aa', 'paralog', 'paralog_residue', 'aa', 'CA_distance_A', 'pLDDT', 'flags', 'sidechain_ok', 'n_agree_of_3', 'consensus_residue', 'best_ref', 'residue_bestref', 'CA_dist_bestref_A', 'residue_domainwise', 'CA_dist_domainwise_A', 'residue_seqaln', 'TM_7P04', 'RMSD_7P04_A', 'note']

def write(path, rows):
    with open(path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join((str(r[c]) for c in cols)) + '\n')

def ready(name):
    return all((load(name, r) for r in ('7P04', '7P03', '7P06'))) and load_dom(name) is not None
main, ctx, exp = ([], [], [])
for name in AF:
    if ready(name):
        main += map_one(name, 'C-AlphaFold', True)
    else:
        print('SKIP (alignments not ready):', name)
write(os.path.join(L.OUT, 'site_class_map.tsv'), main)
for name in CONTEXT:
    if ready(name):
        ctx += map_one(name, 'C-AlphaFold', True)
    else:
        print('SKIP (alignments not ready):', name)
write(os.path.join(L.OUT, 'site_class_map_context_outgroups.tsv'), ctx)
for name in EXP:
    if ready(name):
        exp += map_one(name, 'E-cryoEM', False)
    else:
        print('SKIP (alignments not ready):', name)
write(os.path.join(L.OUT, 'site_class_map_Cdr1_experimental.tsv'), exp)

def summarise(rows, path):
    agg = collections.OrderedDict()
    for r in rows:
        k = (r['paralog'], r['class'])
        a = agg.setdefault(k, {'n': 0, 'mapped': 0, 'NA': 0, 'noal': 0, 'gt3': 0, 'lt70': 0, 'p7080': 0, 'ge80': 0, 'pl': [], 'sc': 0, 'unc': 0, 'diff': 0})
        a['n'] += 1
        fl = r['flags'].split(';')
        if r['paralog_residue'] != '':
            a['mapped'] += 1
        a['NA'] += 'NOT_ALIGNABLE' in fl
        a['noal'] += 'NO_ALIGNED_RESIDUE' in fl
        a['gt3'] += 'MISALIGNED_GT3A' in fl
        a['unc'] += 'MAPPING_UNCERTAIN' in fl
        a['diff'] += 'MAPPING_DIFFERS' in fl
        if r['pLDDT'] != '':
            v = float(r['pLDDT'])
            a['pl'].append(v)
            a['lt70' if v < 70 else 'p7080' if v < 80 else 'ge80'] += 1
        a['sc'] += r['sidechain_ok'] == 'Y'
    with open(path, 'w') as f:
        f.write('paralog\tclass\tevidence\tn_positions\tn_mapped\tn_not_alignable\tn_no_aligned_residue\tn_CA_gt3A\tn_mapping_uncertain\tn_mapping_differs\tn_pLDDT_lt70\tn_pLDDT_70_80\tn_pLDDT_ge80\tmean_pLDDT_mapped\tn_sidechain_ok\n')
        for (p, c), a in agg.items():
            mp = '%.1f' % np.mean(a['pl']) if a['pl'] else ''
            f.write(f"{p}\t{c}\t{rows[0]['evidence']}\t{a['n']}\t{a['mapped']}\t{a['NA']}\t{a['noal']}\t{a['gt3']}\t{a['unc']}\t{a['diff']}\t{a['lt70']}\t{a['p7080']}\t{a['ge80']}\t{mp}\t{a['sc']}\n")
if main:
    summarise(main, os.path.join(L.OUT, 'site_class_map_summary_by_class.tsv'))
if ctx:
    summarise(ctx, os.path.join(L.OUT, 'site_class_map_summary_context.tsv'))
if exp:
    summarise(exp, os.path.join(L.OUT, 'site_class_map_summary_Cdr1_experimental.tsv'))
allrows = main + ctx + exp
with open(os.path.join(L.OUT, 'site_class_map_coverage_per_paralog.tsv'), 'w') as f:
    f.write('model\tevidence\tdistinct_class_residues\tmapped\tmapped_CA_le3A\tnot_alignable\tno_aligned_residue\tflags_OK_only\tsidechain_ok_Y\tmapped_n_agree_ge2\tPdr5_AF_selftest_primary_eq_identity\tPdr5_AF_selftest_consensus_eq_identity\n')
    for name in list(AF) + list(CONTEXT) + EXP:
        rows = [r for r in allrows if r['paralog'] == name]
        if not rows:
            continue
        seen = {r['P33302_residue']: r for r in rows}
        n = len(seen)
        mapped = sum((1 for r in seen.values() if r['paralog_residue'] != ''))
        ok3 = sum((1 for r in seen.values() if r['paralog_residue'] != '' and float(r['CA_distance_A']) <= 3.0))
        na = sum((1 for r in seen.values() if 'NOT_ALIGNABLE' in r['flags']))
        noal = sum((1 for r in seen.values() if 'NO_ALIGNED_RESIDUE' in r['flags']))
        clean = sum((1 for r in seen.values() if r['flags'] == 'OK'))
        sc = sum((1 for r in seen.values() if r['sidechain_ok'] == 'Y'))
        ag2 = sum((1 for r in seen.values() if r['paralog_residue'] != '' and r['n_agree_of_3'] >= 2))
        s1 = s2 = ''
        if name == 'Pdr5':
            s1 = sum((1 for r in seen.values() if r['paralog_residue'] != '' and int(r['paralog_residue']) == int(r['P33302_residue'])))
            s2 = sum((1 for r in seen.values() if r['consensus_residue'] != '' and int(r['consensus_residue']) == int(r['P33302_residue'])))
        f.write(f"{name}\t{rows[0]['evidence']}\t{n}\t{mapped}\t{ok3}\t{na}\t{noal}\t{clean}\t{sc}\t{ag2}\t{s1}\t{s2}\n")
        print(name, dict(n=n, mapped=mapped, ok3=ok3, na=na, noal=noal, clean=clean, sc=sc, ag2=ag2, selftest=(s1, s2)))

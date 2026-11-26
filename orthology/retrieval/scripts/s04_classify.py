import argparse, collections
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from common import read_tsv, write_tsv, log
LEN_MIN, LEN_MAX = (1100, 1800)
FAM = ['PF00005', 'PF01061', 'PF06422', 'PF14510', 'PF19055']
LETTER = {'PF00005': 'A', 'PF01061': 'G'}

def classify(counts, length):
    a, g, c, n = (counts['PF00005'], counts['PF01061'], counts['PF06422'], counts['PF14510'])
    dom_ok = a >= 2 and g >= 2 and (c >= 1)
    len_ok = LEN_MIN <= length <= LEN_MAX
    if dom_ok and len_ok:
        return ('PASS', '')
    reasons = []
    if not dom_ok:
        if a < 2:
            reasons.append('PF00005<2')
        if g < 2:
            reasons.append('PF01061<2')
        if c < 1:
            reasons.append('PF06422=0')
    if not len_ok:
        reasons.append('LEN_SHORT' if length < LEN_MIN else 'LEN_LONG')
    if g + c + n == 0:
        return ('OTHER_ABC', ';'.join(reasons))
    if dom_ok and (not len_ok):
        cls = 'LEN_SHORT' if length < LEN_MIN else 'LEN_LONG'
    elif a >= 2 and g >= 2 and (c == 0) and len_ok:
        cls = 'NO_PF06422'
    elif a == 1 and g == 1:
        cls = 'HALF_SIZE_ABCG'
    else:
        cls = 'PARTIAL_DOMAINS'
    return ('NEAR_MISS:' + cls, ';'.join(reasons))

def merged_count(doms, gap=20):
    doms = sorted(doms)
    n, last_end = (0, None)
    for s, e in doms:
        if last_end is None or s - last_end > gap:
            n += 1
        last_end = e if last_end is None else max(last_end, e)
    return n
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--hits', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    by = collections.OrderedDict()
    for r in read_tsv(a.hits):
        d = by.setdefault(r['seq_id'], dict(species_code=r['species_code'], gene_id=r['gene_id'], seq_len=int(r['seq_len']), doms=[]))
        d['doms'].append((r['hmm_acc'], int(r['ali_from']), int(r['ali_to']), int(r['env_from']), int(r['env_to']), float(r['dom_score'])))
    rows = []
    for sid, d in by.items():
        cnt = {f: 0 for f in FAM}
        coords = {f: [] for f in FAM}
        for acc, af, at, ef, et, sc in sorted(d['doms'], key=lambda x: x[1]):
            cnt[acc] += 1
            coords[acc].append(f'{af}-{at}')
        order = [x[0] for x in sorted(d['doms'], key=lambda x: x[1])]
        arch = ','.join(order)
        ag = ''.join((LETTER[f] for f in order if f in LETTER))
        status, reasons = classify(cnt, d['seq_len'])
        flags = []
        if status == 'PASS':
            if ag != 'AGAG':
                flags.append('ARCH_NOT_AGAG:' + ag)
            envs = {f: [(x[3], x[4]) for x in d['doms'] if x[0] == f] for f in ('PF00005', 'PF01061')}
            if merged_count(envs['PF00005']) < 2 or merged_count(envs['PF01061']) < 2:
                flags.append('PASS_ONLY_VIA_SPLIT_DOMAIN')
        rep = ''
        if status.startswith('NEAR_MISS') and cnt['PF00005'] >= 1 and (cnt['PF01061'] >= 1) and (d['seq_len'] < LEN_MIN) and (cnt['PF06422'] >= 1 or cnt['PF00005'] >= 2 or cnt['PF01061'] >= 2):
            rep = 'Y'
        rows.append(dict(seq_id=sid, species_code=d['species_code'], gene_id=d['gene_id'], seq_len=d['seq_len'], n_PF00005=cnt['PF00005'], n_PF01061=cnt['PF01061'], n_PF06422=cnt['PF06422'], n_PF14510=cnt['PF14510'], n_PF19055=cnt['PF19055'], coords_PF00005=';'.join(coords['PF00005']), coords_PF01061=';'.join(coords['PF01061']), coords_PF06422=';'.join(coords['PF06422']), coords_PF14510=';'.join(coords['PF14510']), coords_PF19055=';'.join(coords['PF19055']), arch=arch, arch_AG=ag, PF14510_present='Y' if cnt['PF14510'] else 'N', status=status, near_miss_reasons=reasons, repair_candidate=rep, rule_flags=';'.join(flags)))
    write_tsv(a.out, rows)
    c = collections.Counter((r['status'] for r in rows))
    log('proteins with >=1 included domain:', len(rows), dict(c))

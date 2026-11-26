from __future__ import annotations
import csv
import hashlib
import math
import random
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
AA = 'ARNDCQEGHILKMFPSTWYV'
AA_INDEX = {a: i for i, a in enumerate(AA)}
_BG_RAW = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]
BG_BLOSUM62 = [x / sum(_BG_RAW) for x in _BG_RAW]
PSEUDOCOUNT = 1e-07
GAP_NA_THRESHOLD = 0.5
SEED_BASE = 20261118
N_DRAWS = 10000
N_BOOT = 2000
EXPOSED_RSA = 0.05
MIN_CANDIDATES = 10
MIN_CLASS_COLUMNS = 10
CALIPER_TIERS = [(0.1, 0.1, 3.0), (0.15, 0.15, 4.5), (0.2, 0.2, 6.0)]
FALLBACK_WEAK_FRACTION = 0.25
COMP_GROUP = {}
for _grp, _members in {'ALI': 'VLIM', 'ARO': 'FWY', 'SML': 'AGSTC', 'POL': 'NQHP', 'CHG': 'DEKR'}.items():
    for _a in _members:
        COMP_GROUP[_a] = _grp

def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def read_fasta(path: str) -> Tuple[List[str], List[str]]:
    ids: List[str] = []
    seqs: List[str] = []
    cur: List[str] = []
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if ids:
                    seqs.append(''.join(cur).upper())
                    cur = []
                ids.append(line[1:].split()[0])
            elif line:
                cur.append(line.strip())
    if ids:
        seqs.append(''.join(cur).upper())
    if len(ids) != len(seqs):
        raise ValueError('FASTA parse error: ids/seqs mismatch')
    if len(set(ids)) != len(ids):
        raise ValueError('duplicate ids in alignment')
    L = {len(s) for s in seqs}
    if len(L) > 1:
        raise ValueError(f'alignment rows have different lengths: {sorted(L)[:5]}')
    return (ids, seqs)

def read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, 'r', encoding='utf-8', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def select_rows(ids: List[str], seqs: List[str], wanted: Sequence[str]) -> Tuple[List[str], List[str]]:
    pos = {i: k for k, i in enumerate(ids)}
    missing = [w for w in wanted if w not in pos]
    if missing:
        raise KeyError(f'{len(missing)} requested rows are not in the alignment, e.g. {missing[:3]}')
    return (list(wanted), [seqs[pos[w]] for w in wanted])

def is_residue(ch: str) -> bool:
    return ch in AA_INDEX

def columns_of(seqs: Sequence[str]):
    return zip(*seqs)

def gap_fractions(seqs: Sequence[str]) -> List[float]:
    n = len(seqs)
    out = []
    for col in columns_of(seqs):
        miss = sum((1 for ch in col if ch not in AA_INDEX))
        out.append(miss / n)
    return out

def reference_columns(ref_seq: str) -> Dict[int, int]:
    out: Dict[int, int] = {}
    k = 0
    for c, ch in enumerate(ref_seq):
        if ch not in '-.':
            k += 1
            out[k] = c
    return out

def check_reference(ref_seq: str, residue_aa: Dict[int, str]) -> List[str]:
    rc = reference_columns(ref_seq)
    bad = []
    for r, a in residue_aa.items():
        if r not in rc:
            bad.append(f'{r}{a}: beyond reference length')
        elif ref_seq[rc[r]] != a:
            bad.append(f'{r}: expected {a}, alignment has {ref_seq[rc[r]]}')
    return bad

def henikoff_weights(seqs: Sequence[str], columns: Optional[Iterable[int]]=None) -> List[float]:
    n = len(seqs)
    w = [0.0] * n
    cols = range(len(seqs[0])) if columns is None else columns
    for c in cols:
        cnt: Dict[str, int] = {}
        for s in seqs:
            ch = s[c]
            if ch in AA_INDEX:
                cnt[ch] = cnt.get(ch, 0) + 1
        r = len(cnt)
        if r == 0:
            continue
        for j, s in enumerate(seqs):
            ch = s[c]
            if ch in AA_INDEX:
                w[j] += 1.0 / (r * cnt[ch])
    tot = sum(w)
    if tot <= 0:
        raise ValueError('all Henikoff weights are zero')
    return [x / tot for x in w]

def empirical_background(seqs: Sequence[str], weights: Sequence[float], columns: Iterable[int]) -> List[float]:
    f = [0.0] * 20
    for c in columns:
        for s, wj in zip(seqs, weights):
            ch = s[c]
            if ch in AA_INDEX:
                f[AA_INDEX[ch]] += wj
    t = sum(f)
    return [x / t for x in f]

def jsd_column(col: Sequence[str], weights: Sequence[float], bg: Sequence[float], gap_penalty: bool=True) -> float:
    fc = [PSEUDOCOUNT] * 20
    gap_w = 0.0
    for ch, wj in zip(col, weights):
        k = AA_INDEX.get(ch)
        if k is None:
            gap_w += wj
        else:
            fc[k] += wj
    tw = sum(weights)
    s = sum(fc)
    fc = [x / s for x in fc]
    d = 0.0
    for p, q in zip(fc, bg):
        m = 0.5 * (p + q)
        if m > 0.0:
            if p > 0.0:
                d += p * math.log2(p / m)
            if q > 0.0:
                d += q * math.log2(q / m)
    d *= 0.5
    if gap_penalty:
        d *= 1.0 - gap_w / tw
    return d

def js_conservation_columns(seqs: Sequence[str], columns: Iterable[int], weights: Sequence[float], bg: Sequence[float]=BG_BLOSUM62, gap_penalty: bool=True) -> Dict[int, float]:
    out: Dict[int, float] = {}
    for c in columns:
        col = [s[c] for s in seqs]
        out[c] = jsd_column(col, weights, bg, gap_penalty)
    return out

def window_scores(res_scores: Dict[int, float], half_width: int, lam: float=0.5) -> Dict[int, float]:
    if half_width <= 0:
        return dict(res_scores)
    keys = sorted(res_scores)
    lo, hi = (keys[0], keys[-1])
    out = dict(res_scores)
    for r in keys:
        if r < lo + half_width or r > hi - half_width:
            continue
        nb = [res_scores[j] for j in range(r - half_width, r + half_width + 1) if j != r and j in res_scores]
        if nb:
            out[r] = (1.0 - lam) * (sum(nb) / len(nb)) + lam * res_scores[r]
    return out

def load_covariates(path: str) -> Dict[int, Dict[str, object]]:
    cov: Dict[int, Dict[str, object]] = {}
    for r in read_tsv(path):
        k = int(r['residue'])
        cov[k] = {'aa': r['aa'], 'rsa': float(r['RSA']), 'z': float(r['z_CA_OPM_A']), 'helix': int(r['in_HELX_P']), 'domain': r['domain'], 'pool': r['pool_label']}
    return cov

def candidate_sets(class_res: Sequence[int], universe: Sequence[int], cov: Dict[int, Dict[str, object]], gap: Optional[Dict[int, float]]=None, comp_match: bool=False, exposed_only: bool=False) -> Tuple[Dict[int, List[int]], Dict[int, str]]:
    uni = [u for u in universe if not exposed_only or cov[u]['rsa'] >= EXPOSED_RSA]
    out: Dict[int, List[int]] = {}
    level: Dict[int, str] = {}
    for i in class_res:
        ci = cov[i]
        e_i = ci['rsa'] >= EXPOSED_RSA
        g_i = None if gap is None else gap.get(i)
        base = [u for u in uni if u != i and cov[u]['pool'] == ci['pool'] and (cov[u]['helix'] == ci['helix']) and (exposed_only or (cov[u]['rsa'] >= EXPOSED_RSA) == e_i)]
        tag_extra = ''
        if comp_match:
            grp = COMP_GROUP.get(ci['aa'])
            sub = [u for u in base if COMP_GROUP.get(cov[u]['aa']) == grp]
            if len(sub) >= MIN_CANDIDATES:
                base = sub
            else:
                tag_extra = '_nocomp'
        chosen: List[int] = []
        lvl = ''
        for t, (cr, cg, cz) in enumerate(CALIPER_TIERS):
            cand = []
            for u in base:
                cu = cov[u]
                if abs(cu['rsa'] - ci['rsa']) > cr:
                    continue
                if gap is not None and g_i is not None and (u in gap) and (abs(gap[u] - g_i) > cg):
                    continue
                if gap is not None and u not in gap:
                    continue
                if ci['pool'] in ('TM', 'TMDx'):
                    if abs(cu['z'] - ci['z']) > cz:
                        continue
                elif cu['domain'] != ci['domain']:
                    continue
                cand.append(u)
            if len(cand) >= MIN_CANDIDATES:
                chosen = cand
                lvl = f'tier{t}'
                break
        if not chosen:

            def dist(u: int) -> float:
                cu = cov[u]
                dz = cu['z'] - ci['z'] if ci['pool'] in ('TM', 'TMDx') else 0.0
                dg = 0.0 if gap is None or g_i is None or u not in gap else gap[u] - g_i
                return math.sqrt(((cu['rsa'] - ci['rsa']) / 0.1) ** 2 + (dz / 3.0) ** 2 + (dg / 0.1) ** 2)
            pool_same = [u for u in base if gap is None or u in gap]
            pool_same.sort(key=dist)
            chosen = pool_same[:MIN_CANDIDATES]
            lvl = 'fallback'
        out[i] = chosen
        level[i] = lvl + tag_extra
    return (out, level)

def draw_controls(order_rng: random.Random, class_res: Sequence[int], cands: Dict[int, List[int]]) -> Tuple[List[int], int]:
    order = list(class_res)
    order_rng.shuffle(order)
    used = set()
    chosen: List[int] = []
    forced = 0
    for i in order:
        avail = [u for u in cands[i] if u not in used]
        if not avail:
            avail = cands[i]
            forced += 1
        u = order_rng.choice(avail)
        used.add(u)
        chosen.append(u)
    return (chosen, forced)

def permutation_test(class_res: Sequence[int], universe: Sequence[int], cov: Dict[int, Dict[str, object]], score: Dict[int, float], gap: Optional[Dict[int, float]]=None, n_draws: int=N_DRAWS, seed: int=SEED_BASE, n_boot: int=N_BOOT, comp_match: bool=False, exposed_only: bool=False, tail: str='upper') -> Dict[str, object]:
    cls = [r for r in class_res if r in score]
    n_req, n_eff = (len(class_res), len(cls))
    res: Dict[str, object] = {'n_requested': n_req, 'n_eff': n_eff, 'dropped_NA': sorted(set(class_res) - set(cls))}
    if n_eff < MIN_CLASS_COLUMNS:
        res.update({'status': 'NOT_TESTABLE', 'reason': f'fewer than {MIN_CLASS_COLUMNS} scorable class columns'})
        return res
    uni = [u for u in universe if u in score]
    cands, level = candidate_sets(cls, uni, cov, gap, comp_match=comp_match, exposed_only=exposed_only)
    empty = [i for i in cls if not cands[i]]
    if empty:
        res.update({'status': 'NOT_TESTABLE', 'reason': f'no candidate controls for residues {empty}'})
        return res
    obs = sum((score[r] for r in cls)) / n_eff
    rng = random.Random(seed)
    cm: List[float] = []
    forced_total = 0
    for _ in range(n_draws):
        ctrl, f = draw_controls(rng, cls, cands)
        forced_total += f
        cm.append(sum((score[u] for u in ctrl)) / n_eff)
    mean_c = sum(cm) / n_draws
    sd_c = math.sqrt(sum(((x - mean_c) ** 2 for x in cm)) / (n_draws - 1))
    eps = 1e-12
    if tail == 'upper':
        b = sum((1 for x in cm if x >= obs - eps))
    else:
        b = sum((1 for x in cm if x <= obs + eps))
    p = (1 + b) / (1 + n_draws)
    delta = obs - mean_c
    z = delta / sd_c if sd_c > 0 else float('nan')
    brng = random.Random(seed + 1)
    vals = [score[r] for r in cls]
    bs = []
    for _ in range(n_boot):
        bs.append(sum((vals[brng.randrange(n_eff)] for _ in range(n_eff))) / n_eff - mean_c)
    bs.sort()
    lo, hi = (bs[int(0.025 * n_boot)], bs[int(0.975 * n_boot) - 1])
    n_fb = sum((1 for v in level.values() if v.startswith('fallback')))
    res.update({'status': 'OK', 'tail': tail, 'observed_mean': obs, 'control_mean': mean_c, 'control_sd': sd_c, 'delta': delta, 'z': z, 'ci95_delta': (lo, hi), 'b': b, 'p': p, 'n_draws': n_draws, 'seed': seed, 'n_fallback_residues': n_fb, 'weakly_matched': n_fb / n_eff > FALLBACK_WEAK_FRACTION, 'match_levels': {k: level[k] for k in cls}, 'n_forced_replacement_events': forced_total, 'mde_80pct_power_alpha05_onesided': (1.6449 + 0.8416) * sd_c})
    return res

def robust_label_h6(runs: Dict[str, Dict[str, object]]) -> Dict[str, object]:

    def ok(name: str) -> bool:
        r = runs.get(name)
        return bool(r) and r.get('status') == 'OK'

    def pos(name: str) -> bool:
        return ok(name) and runs[name]['delta'] > 0
    checks: Dict[str, bool] = {}
    checks['1_significance_R0_p<0.05'] = ok('R0') and runs['R0']['p'] < 0.05 and (runs['R0']['delta'] > 0)
    checks['2_both_aligners_delta>0_and_p<0.05'] = checks['1_significance_R0_p<0.05'] and ok('R1') and (runs['R1']['p'] < 0.05) and (runs['R1']['delta'] > 0)
    checks['3_subset_sign_both_aligners'] = pos('R2') and pos('R3')
    checks['4_score_variants_sign_>=3_of_4'] = sum((pos(n) for n in ('S_unweighted', 'S_nopenalty', 'S_window', 'S_empbg'))) >= 3
    checks['5_confound_controls_sign_comp_exposed_indep'] = pos('C_comp') and pos('C_exposed') and pos('T_indep')
    if not checks['1_significance_R0_p<0.05']:
        label = 'not supported'
    elif all(checks.values()):
        label = 'robust'
    else:
        label = 'supported, not robust'
    return {'label': label, 'checks': checks, 'failed': [k for k, v in checks.items() if not v]}

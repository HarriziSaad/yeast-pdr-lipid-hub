from __future__ import annotations
import csv
import gzip
import hashlib
import importlib.util
import math
import random
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
GENES = [('PDR5', 'YOR153W'), ('SNQ2', 'YDR011W'), ('PDR15', 'YDR406W'), ('PDR10', 'YOR328W'), ('PDR12', 'YPL058C'), ('PDR11', 'YIL013C'), ('PDR18', 'YNR070W'), ('AUS1', 'YOR011W')]
GENE_LABELS = [g for g, _ in GENES]
SYSTEMATIC = dict(GENES)
AA20 = 'ACDEFGHIKLMNPQRSTVWY'
NA_CHARS = set('XBZJUOxbzjuo?-.')
MIN_COVERAGE = 0.9
NW_FALLBACK_DELTA = 5
LOW_IDENTITY = 0.8
H5_MIN_CLASS_COLUMNS = 3
H5_UNDERPOWERED_ALLELES = 10
H5_UNDERPOWERED_COLUMNS = 10
N_CLUSTERS_DEFAULT = 20
BASES = 'TCAG'
_AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate(((x, y, z) for x in BASES for y in BASES for z in BASES))}

def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def load_module(path: str, name: str, expected_sha256: Optional[str]=None):
    if expected_sha256 is not None and sha256_file(path).lower() != expected_sha256.lower():
        raise RuntimeError(f'hash mismatch for {path}')
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def parse_fasta_text(text: str) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    hdr, buf = (None, [])
    for line in text.splitlines():
        if line.startswith('>'):
            if hdr is not None:
                out.append((hdr, ''.join(buf)))
            hdr, buf = (line[1:].rstrip(), [])
        elif line.strip():
            buf.append(line.strip())
    if hdr is not None:
        out.append((hdr, ''.join(buf)))
    return out

def read_fasta_file(path: str) -> List[Tuple[str, str]]:
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt', encoding='utf-8', errors='replace') as fh:
        return parse_fasta_text(fh.read())

def isolate_of(header: str, systematic: str) -> Optional[str]:
    first = header.split()[0]
    key = '_' + systematic + '_'
    i = first.find(key)
    return first[:i] if i > 0 else None

def translate(cds: str) -> str:
    cds = cds.upper()
    return ''.join((CODON.get(cds[i:i + 3], 'X') for i in range(0, len(cds) - len(cds) % 3, 3)))

def missense_possible(codon: str) -> int:
    codon = codon.upper()
    aa = CODON[codon]
    if aa == '*':
        return 0
    n = 0
    for pos in range(3):
        for b in 'ACGT':
            if b == codon[pos]:
                continue
            a2 = CODON[codon[:pos] + b + codon[pos + 1:]]
            if a2 != '*' and a2 != aa:
                n += 1
    return n

def o_per_position(cds: str, ref_prot: str) -> List[int]:
    cds = ''.join(cds.split()).upper()
    L = len(ref_prot)
    if len(cds) not in (3 * L, 3 * L + 3):
        raise ValueError(f'CDS length {len(cds)} does not fit protein length {L}')
    prot = translate(cds[:3 * L])
    if prot != ref_prot:
        bad = [i + 1 for i, (a, b) in enumerate(zip(prot, ref_prot)) if a != b][:5]
        raise ValueError(f'CDS translation differs from the reference protein at {bad}')
    return [missense_possible(cds[3 * i:3 * i + 3]) for i in range(L)]

def nw_align(a: str, b: str, match: int=2, mismatch: int=-1, gap_open: int=-5, gap_ext: int=-1) -> Tuple[str, str]:
    n, m = (len(a), len(b))
    NEG = -10 ** 9
    ptr = [[bytearray(m + 1) for _ in range(n + 1)] for _ in range(3)]
    Mp = [NEG] * (m + 1)
    Xp = [NEG] * (m + 1)
    Yp = [NEG] * (m + 1)
    Mp[0] = 0
    for j in range(1, m + 1):
        Yp[j] = gap_open + (j - 1) * gap_ext
        ptr[2][0][j] = 2 if j > 1 else 0
    for i in range(1, n + 1):
        Mc = [NEG] * (m + 1)
        Xc = [NEG] * (m + 1)
        Yc = [NEG] * (m + 1)
        Xc[0] = gap_open + (i - 1) * gap_ext
        ptr[1][i][0] = 1 if i > 1 else 0
        ai = a[i - 1]
        for j in range(1, m + 1):
            s = match if ai == b[j - 1] else mismatch
            pm, px, py = (Mp[j - 1], Xp[j - 1], Yp[j - 1])
            if pm >= px and pm >= py:
                Mc[j] = pm + s
                ptr[0][i][j] = 0
            elif px >= py:
                Mc[j] = px + s
                ptr[0][i][j] = 1
            else:
                Mc[j] = py + s
                ptr[0][i][j] = 2
            c0, c1, c2 = (Mp[j] + gap_open, Xp[j] + gap_ext, Yp[j] + gap_open)
            if c1 >= c0 and c1 >= c2:
                Xc[j] = c1
                ptr[1][i][j] = 1
            elif c0 >= c2:
                Xc[j] = c0
                ptr[1][i][j] = 0
            else:
                Xc[j] = c2
                ptr[1][i][j] = 2
            d0, d1, d2 = (Mc[j - 1] + gap_open, Yc[j - 1] + gap_ext, Xc[j - 1] + gap_open)
            if d1 >= d0 and d1 >= d2:
                Yc[j] = d1
                ptr[2][i][j] = 2
            elif d0 >= d2:
                Yc[j] = d0
                ptr[2][i][j] = 0
            else:
                Yc[j] = d2
                ptr[2][i][j] = 1
        Mp, Xp, Yp = (Mc, Xc, Yc)
    state = max(range(3), key=lambda k: (Mp[m], Xp[m], Yp[m])[k])
    i, j = (n, m)
    ra, rb = ([], [])
    while i > 0 or j > 0:
        nxt = ptr[state][i][j]
        if state == 0:
            ra.append(a[i - 1])
            rb.append(b[j - 1])
            i -= 1
            j -= 1
        elif state == 1:
            ra.append(a[i - 1])
            rb.append('-')
            i -= 1
        else:
            ra.append('-')
            rb.append(b[j - 1])
            j -= 1
        state = nxt
    return (''.join(reversed(ra)), ''.join(reversed(rb)))

def call_states(ref: str, rec: str) -> Tuple[str, Dict[str, object]]:
    L = len(ref)
    core = rec[:-1] if rec.endswith('*') else rec
    qc: Dict[str, object] = {'len_core': len(core), 'mode': 'direct', 'flag': ''}
    if abs(len(core) - L) > NW_FALLBACK_DELTA:
        ra, rb = nw_align(core, ref)
        out = []
        for x, y in zip(ra, rb):
            if y == '-':
                continue
            out.append(x)
        seq = ''.join(out)
        qc['mode'] = 'nw'
    else:
        seq = core[:L]
        if len(core) > L:
            qc['flag'] = 'EXTENSION'
        elif len(core) < L:
            seq = seq + '?' * (L - len(core))
            qc['flag'] = 'SHORT'
    diff = []
    same = cmp_ = 0
    for r, t in zip(ref, seq):
        if t == '*':
            diff.append('*')
        elif t in NA_CHARS or t not in AA20:
            diff.append('?')
        elif t == r:
            diff.append('.')
            same += 1
            cmp_ += 1
        else:
            diff.append(t)
            cmp_ += 1
    if cmp_ and qc['mode'] == 'direct' and (same / cmp_ < LOW_IDENTITY):
        qc['flag'] = 'LOWID'
        diff = ['?'] * L
    return (''.join(diff), qc)

class GeneVariation:

    def __init__(self, gene: str, ref: str, diffs: Dict[str, str], n_records: int, qc_flags: Dict[str, int]):
        self.gene = gene
        self.ref = ref
        self.L = len(ref)
        self.diffs = diffs
        self.n_records = n_records
        self.qc_flags = qc_flags
        self._tally()

    def _tally(self):
        L = self.L
        self.n_data = [0] * L
        self.alt: List[Dict[str, int]] = [dict() for _ in range(L)]
        self.n_stop = [0] * L
        for s, d in self.diffs.items():
            for j, ch in enumerate(d):
                if ch == '.':
                    self.n_data[j] += 1
                elif ch == '*':
                    self.n_stop[j] += 1
                elif ch != '?':
                    self.n_data[j] += 1
                    a = self.alt[j]
                    a[ch] = a.get(ch, 0) + 1

    def usable(self, j: int, min_cov: float=MIN_COVERAGE) -> bool:
        return self.n_records > 0 and self.n_data[j] / self.n_records >= min_cov

    def m(self, j: int, mode: str='distinct') -> float:
        a = self.alt[j]
        if mode == 'distinct':
            return float(len(a))
        if mode == 'nosingle':
            return float(sum((1 for c in a.values() if c >= 2)))
        if mode == 'af':
            return sum(a.values()) / self.n_data[j] if self.n_data[j] else 0.0
        raise ValueError(mode)

def gene_variation(gene: str, ref: str, records: Sequence[Tuple[str, str]], systematic: Optional[str]=None, mask_downstream_of_stop: bool=False) -> GeneVariation:
    systematic = systematic or SYSTEMATIC.get(gene, gene)
    diffs: Dict[str, str] = {}
    flags: Dict[str, int] = {}
    for hdr, seq in records:
        iso = isolate_of(hdr, systematic)
        if iso is None:
            flags['BAD_HEADER'] = flags.get('BAD_HEADER', 0) + 1
            continue
        if iso in diffs:
            raise ValueError(f'duplicate isolate code {iso} in gene {gene}')
        d, qc = call_states(ref, seq)
        if qc['flag']:
            flags[qc['flag']] = flags.get(qc['flag'], 0) + 1
        if qc['mode'] == 'nw':
            flags['NW'] = flags.get('NW', 0) + 1
        if mask_downstream_of_stop and '*' in d:
            k = d.index('*')
            d = d[:k + 1] + '?' * (len(d) - k - 1)
        diffs[iso] = d
    return GeneVariation(gene, ref, diffs, len(diffs), flags)

def subset_variation(gv: GeneVariation, strains: Iterable[str]) -> GeneVariation:
    keep = set(strains)
    diffs = {s: d for s, d in gv.diffs.items() if s in keep}
    return GeneVariation(gv.gene, gv.ref, diffs, len(diffs), dict(gv.qc_flags))

def build_column_map(aligned: Dict[str, str], pdr5_label: str='PDR5') -> Dict[int, Dict[str, int]]:
    lens = {len(v) for v in aligned.values()}
    if len(lens) != 1:
        raise ValueError('aligned rows differ in length')
    ncol = lens.pop()
    counters = {g: 0 for g in aligned}
    out: Dict[int, Dict[str, int]] = {}
    col_of: Dict[int, int] = {}
    for c in range(ncol):
        here = {}
        for g, row in aligned.items():
            ch = row[c]
            if ch not in '-.':
                counters[g] += 1
                here[g] = counters[g]
        if pdr5_label in here:
            out[here[pdr5_label]] = here
            col_of[here[pdr5_label]] = c
    out_cols = {i: col_of[i] for i in out}
    out['__col__'] = out_cols
    return out

def column_density(colmap: Dict[int, Dict[str, int]], var: Dict[str, GeneVariation], o: Dict[str, List[int]], mode: str='distinct', genes: Optional[Sequence[str]]=None, na_columns: Optional[set]=None, min_cov: float=MIN_COVERAGE) -> Tuple[Dict[int, float], Dict[int, Tuple[float, int, int]]]:
    genes = list(genes) if genes is not None else [g for g in GENE_LABELS if g in var]
    cols = colmap.get('__col__', {})
    dens: Dict[int, float] = {}
    diag: Dict[int, Tuple[float, int, int]] = {}
    for i, per in colmap.items():
        if i == '__col__':
            continue
        if na_columns is not None and cols.get(i) in na_columns:
            continue
        M = 0.0
        O = 0
        ng = 0
        for g in genes:
            j = per.get(g)
            if j is None or g not in var:
                continue
            gv = var[g]
            if not gv.usable(j - 1, min_cov):
                continue
            oo = o[g][j - 1]
            if oo <= 0:
                continue
            M += gv.m(j - 1, mode)
            O += oo
            ng += 1
        if O > 0:
            dens[i] = M / O
            diag[i] = (M, O, ng)
    return (dens, diag)

def h5_test(a11, class_res: Sequence[int], universe: Sequence[int], cov, score: Dict[int, float], gap: Optional[Dict[int, float]], n_draws: int, seed: int, n_boot: int=2000, comp_match: bool=False, exposed_only: bool=False, diag=None) -> Dict[str, object]:
    old = a11.MIN_CLASS_COLUMNS
    a11.MIN_CLASS_COLUMNS = H5_MIN_CLASS_COLUMNS
    try:
        res = a11.permutation_test(list(class_res), list(universe), cov, score, gap, n_draws=n_draws, seed=seed, n_boot=n_boot, comp_match=comp_match, exposed_only=exposed_only, tail='lower')
    finally:
        a11.MIN_CLASS_COLUMNS = old
    if res.get('status') == 'OK':
        mc = res['control_mean']
        res['ratio'] = res['observed_mean'] / mc if mc > 0 else float('nan')
        lo, hi = res['ci95_delta']
        res['ci95_ratio'] = ((lo + mc) / mc, (hi + mc) / mc) if mc > 0 else (float('nan'), float('nan'))
        flags = []
        if res['n_eff'] < H5_UNDERPOWERED_COLUMNS:
            flags.append('UNDERPOWERED_FEW_COLUMNS')
        if diag is not None:
            tot = sum((diag[r][0] for r in class_res if r in diag))
            res['class_sum_m'] = tot
            if tot < H5_UNDERPOWERED_ALLELES:
                flags.append('UNDERPOWERED_FEW_ALLELES')
        res['flags'] = flags
    return res

def bh(pvals: Sequence[float]) -> List[float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * m / rank)
        q[i] = prev
    return q

def read_distance_matrix(path: str) -> Tuple[List[str], 'object']:
    import numpy as np
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt', encoding='utf-8', newline='') as fh:
        rd = csv.reader(fh, delimiter='\t')
        header = next(rd)
        ids = header[1:]
        rows = []
        rid = []
        for r in rd:
            rid.append(r[0])
            rows.append([float(x) for x in r[1:]])
    if rid != ids:
        raise ValueError('distance matrix: row and column identifiers differ')
    D = np.array(rows, dtype=float)
    D = (D + D.T) / 2.0
    return (ids, D)

def cluster_strains(ids: Sequence[str], D, k: int=N_CLUSTERS_DEFAULT) -> Dict[str, int]:
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    import numpy as np
    Dm = np.array(D, dtype=float)
    np.fill_diagonal(Dm, 0.0)
    Z = linkage(squareform(Dm, checks=False), method='average')
    lab = fcluster(Z, t=k, criterion='maxclust')
    groups: Dict[int, List[str]] = {}
    for s, l in zip(ids, lab):
        groups.setdefault(int(l), []).append(s)
    ordered = sorted(groups.values(), key=lambda g: (-len(g), min(g)))
    out: Dict[str, int] = {}
    for n, g in enumerate(ordered, start=1):
        for s in g:
            out[s] = n
    return out

def cluster_bootstrap(a11, class_res: Sequence[int], universe: Sequence[int], cov, gap, colmap, var: Dict[str, GeneVariation], o: Dict[str, List[int]], clusters: Dict[str, int], na_columns: Optional[set], n_rep: int=200, n_draws: int=1000, seed: int=20261118, min_cov: float=MIN_COVERAGE) -> Dict[str, object]:
    full_d, _ = column_density(colmap, var, o, 'distinct', None, na_columns, min_cov)
    cls = [r for r in class_res if r in full_d]
    uni = [u for u in universe if u in full_d]
    if len(cls) < H5_MIN_CLASS_COLUMNS:
        return {'status': 'NOT_TESTABLE'}
    cands, _lvl = a11.candidate_sets(cls, uni, cov, gap)
    rng = random.Random(seed)
    draws = []
    for _ in range(n_draws):
        ctrl, _f = a11.draw_controls(rng, cls, cands)
        draws.append(ctrl)
    need = set(cls)
    for ctrl in draws:
        need.update(ctrl)
    cl_ids = sorted(set(clusters.values()))
    alle: Dict[str, Dict[int, Dict[int, set]]] = {}
    for g, gv in var.items():
        per = {c: {} for c in cl_ids}
        for s, d in gv.diffs.items():
            c = clusters.get(s)
            if c is None:
                continue
            store = per[c]
            for j, ch in enumerate(d):
                if ch not in '.?*':
                    store.setdefault(j, set()).add(ch)
        alle[g] = per
    cols = colmap.get('__col__', {})
    terms: Dict[int, List[Tuple[str, int, int]]] = {}
    for i in need:
        per = colmap.get(i)
        lst = []
        for g, gv in var.items():
            j = per.get(g) if per else None
            if j is None or not gv.usable(j - 1, min_cov) or o[g][j - 1] <= 0:
                continue
            lst.append((g, j - 1, o[g][j - 1]))
        terms[i] = lst
    brng = random.Random(seed + 7)
    ratios, neg = ([], 0)
    for _ in range(n_rep):
        pick = sorted(set((cl_ids[brng.randrange(len(cl_ids))] for _ in cl_ids)))
        d: Dict[int, float] = {}
        for i in need:
            M = 0
            O = 0
            for g, j0, oo in terms[i]:
                u: set = set()
                for c in pick:
                    s_ = alle[g][c].get(j0)
                    if s_:
                        u |= s_
                M += len(u)
                O += oo
            if O > 0:
                d[i] = M / O
        cls_r = [r for r in cls if r in d]
        if len(cls_r) < H5_MIN_CLASS_COLUMNS:
            continue
        obs = sum((d[r] for r in cls_r)) / len(cls_r)
        cm = []
        for ctrl in draws:
            vals = [d[u] for u in ctrl if u in d]
            if vals:
                cm.append(sum(vals) / len(vals))
        mc = sum(cm) / len(cm) if cm else float('nan')
        if mc and mc > 0:
            ratios.append(obs / mc)
            neg += 1 if obs < mc else 0
    ratios.sort()
    n = len(ratios)
    out = {'status': 'OK', 'n_replicates': n, 'share_ratio_lt_1': neg / n if n else float('nan')}
    if n:
        out['ratio_q025'] = ratios[int(0.025 * n)]
        out['ratio_median'] = ratios[n // 2]
        out['ratio_q975'] = ratios[min(n - 1, int(0.975 * n))]
    return out

def rint(values: Sequence[float]) -> List[float]:
    from statistics import NormalDist
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    nd = NormalDist()
    return [nd.inv_cdf((r - 0.5) / n) for r in ranks]

def kinship_from_distance(D):
    import numpy as np
    D = np.array(D, dtype=float)
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    K = -0.5 * J @ D ** 2 @ J
    K = (K + K.T) / 2.0
    w, U = np.linalg.eigh(K)
    w = np.clip(w, 0.0, None)
    tr = w.sum()
    if tr <= 0:
        raise ValueError('degenerate kinship')
    w = w * (n / tr)
    idx = np.argsort(-w)
    return (w[idx], U[:, idx])

def _reml_neg(logdelta: float, ys, Xs, S) -> float:
    import numpy as np
    delta = math.exp(logdelta)
    h = S + delta
    wi = 1.0 / h
    XtWX = (Xs * wi[:, None]).T @ Xs
    try:
        beta = np.linalg.solve(XtWX, (Xs * wi[:, None]).T @ ys)
    except np.linalg.LinAlgError:
        return 1e+18
    r = ys - Xs @ beta
    n, p = Xs.shape
    s2 = float(r * wi @ r) / (n - p)
    if s2 <= 0:
        return 1e+18
    sign, logdet = np.linalg.slogdet(XtWX)
    return 0.5 * ((n - p) * math.log(s2) + float(np.log(h).sum()) + logdet)

def lmm_p3d(y: Sequence[float], x_cols: Sequence[Sequence[float]], S, U, covars: Optional[Sequence[Sequence[float]]]=None) -> Dict[str, object]:
    import numpy as np
    from scipy import optimize, stats
    y = np.asarray(y, dtype=float)
    n = len(y)
    X0 = np.ones((n, 1))
    if covars:
        X0 = np.hstack([X0] + [np.asarray(c, dtype=float).reshape(n, 1) for c in covars])
    Ut = U.T
    ys = Ut @ y
    X0s = Ut @ X0
    res = optimize.minimize_scalar(_reml_neg, bounds=(-8.0, 8.0), args=(ys, X0s, S), method='bounded', options={'xatol': 0.0001})
    delta = math.exp(float(res.x))
    wi = 1.0 / (S + delta)
    out = []
    for xc in x_cols:
        x = np.asarray(xc, dtype=float)
        if np.ptp(x) == 0:
            out.append({'status': 'CONSTANT'})
            continue
        X = np.hstack([X0, x.reshape(n, 1)])
        Xs = Ut @ X
        XtW = (Xs * wi[:, None]).T
        A = XtW @ Xs
        try:
            Ainv = np.linalg.inv(A)
        except np.linalg.LinAlgError:
            out.append({'status': 'SINGULAR'})
            continue
        beta = Ainv @ (XtW @ ys)
        r = ys - Xs @ beta
        p_ = X.shape[1]
        s2 = float(r * wi @ r) / (n - p_)
        se = math.sqrt(s2 * Ainv[-1, -1])
        t = beta[-1] / se
        pval = float(2 * stats.t.sf(abs(t), n - p_))
        out.append({'status': 'OK', 'beta': float(beta[-1]), 'se': se, 't': float(t), 'p': pval, 'delta': delta, 'h2': 1.0 / (1.0 + delta), 'n': n})
    return {'delta': delta, 'tests': out}

def ols_test(y: Sequence[float], x: Sequence[float]) -> Dict[str, float]:
    import numpy as np
    from scipy import stats
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    n = len(y)
    X = np.column_stack([np.ones(n), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    s2 = float(r @ r) / (n - 2)
    se = math.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
    t = beta[1] / se
    return {'beta': float(beta[1]), 'se': se, 'p': float(2 * stats.t.sf(abs(t), n - 2))}

def class_burden(var: Dict[str, GeneVariation], colmap: Dict[int, Dict[str, int]], class_res: Sequence[int], genes: Sequence[str], strains: Sequence[str]) -> Dict[str, int]:
    out = {s: 0 for s in strains}
    for i in class_res:
        per = colmap.get(i)
        if not per:
            continue
        for s in strains:
            hit = False
            for g in genes:
                j = per.get(g)
                if j is None or g not in var:
                    continue
                d = var[g].diffs.get(s)
                if d is None:
                    continue
                ch = d[j - 1]
                if ch not in '.?*':
                    hit = True
                    break
            if hit:
                out[s] += 1
    return out

def stopgain_carriers(gv: GeneVariation) -> Dict[str, int]:
    return {s: 1 if '*' in d else 0 for s, d in gv.diffs.items()}

def label_h5(runs: Dict[str, Dict[str, object]], q_reference: Optional[float]) -> Dict[str, object]:
    fails: List[str] = []
    notes: List[str] = []
    prim = runs.get('PRIMARY', {})
    c1 = prim.get('status') == 'OK' and q_reference is not None and (q_reference < 0.05) and (prim.get('ratio', 9) < 1)
    if not c1:
        return {'label': 'not supported', 'failed': ['1_significance'], 'notes': notes}

    def ok_ratio(name):
        r = runs.get(name, {})
        return r.get('status') == 'OK' and r.get('ratio', 9) < 1
    v = runs.get('V_muscle', {})
    if not (v.get('status') == 'OK' and v.get('ratio', 9) < 1 and (v.get('p', 1) < 0.05)):
        fails.append('2_V_muscle')
    for nm in ('S_nosingle', 'S_af', 'S_nostop', 'S_major', 'S_notmajor', 'S_pdr5only'):
        r = runs.get(nm, {})
        if r.get('status') == 'NOT_TESTABLE' and nm in ('S_major', 'S_notmajor', 'S_pdr5only'):
            notes.append(f'{nm} not evaluable')
            continue
        if not ok_ratio(nm):
            fails.append('3_' + nm)
    for nm in ('C_comp', 'C_exposed'):
        if not ok_ratio(nm):
            fails.append('4_' + nm)
    b = runs.get('BOOT', {})
    if not (b.get('status') == 'OK' and b.get('share_ratio_lt_1', 0) >= 0.9):
        fails.append('5_cluster_bootstrap')
    return {'label': 'robust' if not fails else 'supported, not robust', 'failed': fails, 'notes': notes}

import csv, os, sys

def rd(p):
    with open(p, encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def f(x):
    try:
        return float(x)
    except Exception:
        return float('nan')

def bh(ps):
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    q = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, ps[i] * n / rank)
        q[i] = prev
    return q

def main(run):
    out = []
    mp = os.path.join(run, 'M', 'a12_M_results.tsv')
    if os.path.exists(mp):
        M = rd(mp)
        prim = [r for r in M if r['variant'] == 'PRIMARY']
        fm1 = [r for r in prim if r['test'] in ('M1a', 'M1b')]
        qs = bh([f(r['p']) for r in fm1])
        out.append('FM1 BH check: ' + '; '.join((f"{r['test']} p={f(r['p']):.4f} q(recomputed)={q:.4f} q(driver)={f(r['q']):.4f}" for r, q in zip(fm1, qs))))
        m2 = [r for r in prim if r['test'] == 'M2'][0]
        out.append(f"FM2 (single test): M2 p={f(m2['p']):.4f} q={f(m2['q']):.4f}")
    bp = os.path.join(run, 'B1', 'a12_B1_results.tsv')
    if os.path.exists(bp):
        B = rd(bp)
        prim = [r for r in B if r['variant'] == 'PRIMARY' and r['structure'] == '7P04']
        qs = bh([f(r['p']) for r in prim])
        out.append('FB1 BH check: ' + '; '.join((f"{r['test']} p={f(r['p']):.4f} q(recomputed)={q:.4f} q(driver)={f(r['q']):.4f}" for r, q in zip(prim, qs))))
        for r, q in zip(prim, qs):
            t = r['test']
            checks = []
            for struct, var in (('7P03', 'cutoff 15, modes all'), ('7P06', 'cutoff 15, modes all'), ('7P04', 'cutoff 12, modes all'), ('7P04', 'cutoff 18, modes all')):
                rr = [x for x in B if x['test'] == t and x['structure'] == struct and (x['variant'] == var)]
                checks.append((struct, var, f(rr[0]['delta']) if rr else float('nan')))
            sign0 = f(r['delta']) > 0
            agree = all(((c[2] > 0) == sign0 for c in checks))
            lab = 'consistent' if q < 0.1 and agree else 'no evidence (q >= 0.1 at the primary)' if q >= 0.1 else 'not consistent (q < 0.1 but sign differs)'
            out.append(f"B1 label {t}: primary delta={f(r['delta']):.4g} p={f(r['p']):.4f} q={q:.4f}; sign checks " + ', '.join((f'{s}/{v}: delta={d:.4g}' for s, v, d in checks)) + f' -> {lab}')
    b2 = os.path.join(run, 'B2', 'a12_B2_results.tsv')
    if os.path.exists(b2):
        B = rd(b2)
        prim = [r for r in B if r['variant'] == 'PRIMARY']
        k = sum((1 for r in prim if f(r['p']) < 0.1))
        neg = all((f(r['delta']) < 0 for r in prim))
        out.append(f'B2 consistency: {len(prim)} primary entries, p<0.1 in {k}, delta<0 in all: {neg} -> ' + ('consistent' if len(prim) >= 4 and k >= 3 and neg else 'not consistent'))
    txt = '\n'.join(out)
    open(os.path.join(run, 'analysis', 'a12_label_check.txt'), 'w', encoding='utf-8').write(txt + '\n')
    print(txt)
if __name__ == '__main__':
    main(sys.argv[1])

import os
import sys
if __name__ == '__main__':
    _T = sys.argv[sys.argv.index('--threads') + 1] if '--threads' in sys.argv else '8'
    for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[_v] = _T
SCRIPTS = 'C:\\Users\\MSI\\Documents\\Advisor\\project\\phase5_validation\\A12\\scripts'
sys.path.insert(0, SCRIPTS)
import hashlib
import time
import numpy as np
import a12_b as B
import a12_enm as E
import a12_lib as L
import a12_run_b as RB
import a12_run_b2_fast as D
import a12_b2_fast as F
TOL = 1e-08

def rel(a, b):
    a, b = (float(a), float(b))
    den = max(abs(a), abs(b))
    return 0.0 if den == 0 else abs(a - b) / den

def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()

def main(argv):
    root, outdir = (argv[1], argv[2])
    npar = int(argv[argv.index('--npar') + 1]) if '--npar' in argv else 30
    workers = int(argv[argv.index('--workers') + 1]) if '--workers' in argv else 4
    scratch = argv[argv.index('--scratch') + 1] if '--scratch' in argv else os.environ.get('TEMP', '.')
    os.makedirs(outdir, exist_ok=True)
    out_lines = []

    def log(m=''):
        line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {m}"
        print(line, flush=True)
        out_lines.append(line)
    log(f"B2 PARITY (registered vs fast) npar={npar} workers={workers} main-process BLAS threads {(sys.argv[sys.argv.index('--threads') + 1] if '--threads' in sys.argv else '8')} (set before numpy import)")
    log('registered scripts (SHA-256 against A12 item 12):')
    bad = D.check_registered_scripts(log)
    if bad:
        log('STOP: registered script hash mismatch ' + ','.join(bad))
        return 1
    here = os.path.join(SCRIPTS)
    log('a12_b2_fast.py ' + sha(os.path.join(here, 'a12_b2_fast.py')))
    log('a12_run_b2_fast.py ' + sha(os.path.join(here, 'a12_run_b2_fast.py')))
    L.guard_amendment_appended(RB.P(root, 'preregistration', 'AMENDMENTS.md'))
    cov, sets, structs, hits = RB.load_b_inputs(root)
    tmp_nodes = os.path.join(scratch, 'a12_b2_parity_sterols')
    qc, ster = RB.build_sterols(root, tmp_nodes)
    b2dir = RB.P(root, 'phase5_validation', 'A12', 'run', 'B2')
    for e in ster:
        f = f'a12_sterol_nodes_{e}.tsv'
        same = os.path.exists(os.path.join(b2dir, f)) and sha(os.path.join(b2dir, f)) == sha(os.path.join(tmp_nodes, f))
        log(f"  sterol node table {f}: regenerated file identical to run/B2 copy: {same}; n_mol {qc[e]['n_mol']}, n_atoms {qc[e]['n_atoms']}")
    s = structs['7P04']
    P, res = (s['X'], s['res'])
    t = time.time()
    G = E.covariance(*E.modes(E.anm(P, 15.0)))
    log(f'G (registered functions) computed in {time.time() - t:.1f} s; N = {len(P)} nodes')
    pos = {r: i for i, r in enumerate(res)}
    gate_all = [pos[r] for r in sorted(sets['GATE_ALL']) if r in pos]
    gate_core = [pos[r] for r in sorted(sets['GATE_CORE']) if r in pos]
    anchors, outward = E.lipid_anchors(s['ca'], cov)
    log(f'gate nodes: GATE_ALL {len(gate_all)}, GATE_CORE {len(gate_core)}; anchors {len(anchors)}')
    reg = {}
    t_reg = 0.0
    n_eval = 0
    for e, d in ster.items():
        xyz, mol = (d['xyz'], d['mol'])
        templates = [xyz[mol == k] - xyz[mol == k].mean(0) for k in np.unique(mol)]
        for variant, key, cj in D.variants_for(e):
            nodes = gate_core if key == 'core' else gate_all
            t = time.time()
            o = B.b2_test(P, G, nodes, xyz, templates, anchors, outward, s['heavy'], n_draws=npar, seed_offset=sum(map(ord, e)) % 1000, cut_ps=D.CUTS[cj], cut_ss=10.0)
            t_reg += time.time() - t
            n_eval += 1 + npar
            reg[e, variant] = o
        log(f'  registered code done for {e} ({t_reg:.0f} s so far)')
    log(f'registered code: {t_reg:.1f} s for {n_eval} evaluations (20 observed + {n_eval - 20} null draws, each null draw = one placement + one condensation + Woodbury; single process, {t_reg / n_eval:.3f} s per evaluation)')
    t = time.time()
    rows_reg, rows_all, detail, info = D.run_b2_fast(P, res, cov, s['ca'], sets, s['heavy'], ster, nd=npar, nd_registered=npar, workers=workers, G=G, scratch=scratch, chunk=5, log=log)
    t_fast = time.time() - t
    log(f"fast code: {t_fast:.1f} s wall for the same {len(ster) * 4 * (1 + npar)} (entry, variant, draw) values using {workers} worker processes (phases: {info['phases']})")
    table = []
    worst5 = worst_all = worst_obs = 0.0
    worst_nc = 0
    log('per entry and variant: max relative difference over observed + first 5 null draws | over observed + all null draws | contact-node counts equal')
    for (e, variant), o in reg.items():
        dd = detail[e, variant]
        r_obs = rel(o['obs'], dd['obs'])
        nulls = [rel(a, b) for a, b in zip(o['null_dist'], dd['null'])]
        nc_ok = int(o['n_contact_nodes']) == int(dd['n_contact_nodes'])
        worst_obs = max(worst_obs, r_obs)
        m5 = max([r_obs] + nulls[:5])
        mall = max([r_obs] + nulls)
        worst5, worst_all = (max(worst5, m5), max(worst_all, mall))
        worst_nc += 0 if nc_ok else 1
        log(f'  {e} {variant:<28s} {m5:.2e} | {mall:.2e} | {nc_ok}')
        table.append((e, variant, 'observed', o['obs'], dd['obs'], r_obs))
        for i, (a, b, r) in enumerate(zip(o['null_dist'], dd['null'], nulls)):
            table.append((e, variant, f'null_{i}', a, b, r))
    with open(os.path.join(outdir, 'B2_parity_table.tsv'), 'w', encoding='utf-8', newline='') as fh:
        fh.write('structure\tvariant\tquantity\tregistered_code\tfast_code\trelative_difference\n')
        for r in table:
            fh.write(f'{r[0]}\t{r[1]}\t{r[2]}\t{r[3]!r}\t{r[4]!r}\t{r[5]:.3e}\n')
    n_pl = 50
    tree = F.make_tree(s['heavy'])
    all_eq = True
    for e, d in ster.items():
        xyz, mol = (d['xyz'], d['mol'])
        templates = [xyz[mol == k] - xyz[mol == k].mean(0) for k in np.unique(mol)]
        r1, r2 = (L.make_rng(sum(map(ord, e)) % 1000), L.make_rng(sum(map(ord, e)) % 1000))
        eq = True
        for _ in range(n_pl):
            a, _m = E.place_templates(templates, anchors, outward, s['heavy'], r1, n_mol=len(templates))
            b, _m2 = F.place_templates_fast(templates, anchors, outward, tree, r2, n_mol=len(templates))
            eq &= bool(np.array_equal(a, b))
        eq &= r1.bit_generator.state == r2.bit_generator.state
        all_eq &= eq
        log(f'  placement coordinates bitwise identical over the first {n_pl} draws and RNG state equal: {e}: {eq}')
    log('row-level check: registered a12_run_b.run_b2(nd=5) against the fast driver (nd=5)')
    rows_ref = RB.run_b2(P, res, cov, s['ca'], sets, s['heavy'], ster, nd=5)
    rows5, _, _, _ = D.run_b2_fast(P, res, cov, s['ca'], sets, s['heavy'], ster, nd=5, nd_registered=5, workers=workers, G=G, scratch=scratch, chunk=5, log=lambda m: None)
    worst_rows = 0.0
    str_ok = True
    for a, b in zip(rows_ref, rows5):
        for k in RB.FIELDS:
            if isinstance(a[k], str) or a[k] == '':
                str_ok &= a[k] == b[k]
            else:
                x, y = (float(a[k]), float(b[k]))
                if x != x and y != y:
                    continue
                worst_rows = max(worst_rows, rel(x, y))
    log(f'  {len(rows_ref)} rows; text columns identical: {str_ok}; max relative difference over all numeric columns (obs, null_mean, delta, z, p, mde80): {worst_rows:.2e}')
    ok = worst5 <= TOL and worst_all <= TOL and (worst_nc == 0) and all_eq and str_ok and (worst_rows <= TOL)
    log('SUMMARY')
    log(f'  max relative difference, observed configurations (5 entries x 4 variants): {worst_obs:.3e}')
    log(f'  max relative difference, observed + first 5 null draws (required parity set): {worst5:.3e}')
    log(f'  max relative difference, observed + first {npar} null draws: {worst_all:.3e}')
    log(f'  contact-node count mismatches: {worst_nc}; placements bitwise equal: {all_eq}; row-level max relative difference {worst_rows:.3e}')
    log(f'  speed on this machine (single registered process vs {workers} workers): {t_reg:.1f} s vs {t_fast:.1f} s')
    log(f"VERDICT: {('PARITY PASS (<= 1e-8)' if ok else 'PARITY FAIL -> STOP')}")
    with open(os.path.join(outdir, 'B2_parity_log.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(out_lines) + '\n')
    return 0 if ok else 1
if __name__ == '__main__':
    sys.exit(main(sys.argv))

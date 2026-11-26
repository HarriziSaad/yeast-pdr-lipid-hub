import os
import sys
if __name__ == '__main__':
    _T = '8'
    if '--threads' in sys.argv:
        _T = sys.argv[sys.argv.index('--threads') + 1]
    for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[_v] = _T
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import concurrent.futures as cf
import hashlib
import json
import multiprocessing as mp
import platform
import time
import numpy as np
import a12_b2_fast as F
import a12_enm as E
import a12_lib as L
import a12_run_b as RB
CUTS = (10.0, 8.0, 12.0)
CUT_SS = 10.0
REGISTERED_SHA256 = {'a12_b.py': 'e21870fef239ff4de85c86e16ed4f7f248d69b1e0fe71674bdd671154a4bd704', 'a12_enm.py': '62ed5b821446608b5f26a6812dfc5b3fb2d5e007ce1406ddc28704a61ed7f878', 'a12_inputs.py': 'e7c5656bdd87834e1e3c59e96fa3566d715bb82eca5fa2c0d98fa93f3ff6e72b', 'a12_lib.py': '327cadc6587bfa453ac3e6181dc14e9bc39ea4d4f00d3c9458a6d44604fd430b', 'a12_m.py': '426097d6b9aacba4db76ee2368c4b94719c31555f5c2037dc2a8b1115f578283', 'a12_run_b.py': '53183f34486caf623740089cf578c27c4d6b0c28b14073440bf11f1ae5f59ecd', 'a12_run_m.py': 'b6fac5a724b0c7d6f2bbee62ed8a58666bb2c7d80cde69a940afdc3841385286', 'a12_select.py': '24f5fba32395cc17e71b98e351f5939843a3cc2ba257fac81306d968bd644891', 'a12_sets.py': 'f6338edcee4a06b23579485052d3cd38ebe0375fed66f8cc8cc655c018c969e7'}

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for blk in iter(lambda: fh.read(1 << 20), b''):
            h.update(blk)
    return h.hexdigest()

def check_registered_scripts(log=print):
    here = os.path.dirname(os.path.abspath(__file__))
    bad = []
    for name, want in REGISTERED_SHA256.items():
        got = sha256_file(os.path.join(here, name))
        log(f"  {name} {got} {('OK' if got == want else 'MISMATCH (registered ' + want + ')')}")
        if got != want:
            bad.append(name)
    return bad

def variants_for(entry):
    return [('PRIMARY' if entry in RB.INWARD_ENTRIES else 'outward-facing entry', 'all', 0), ('gate-core only', 'core', 0), ('sterol-protein cutoff 8', 'all', 1), ('sterol-protein cutoff 12', 'all', 2)]

def compute_G(P, cutoff=15.0, log=print):
    t = time.time()
    G = E.covariance(*E.modes(E.anm(P, cutoff)))
    log(f'  G = H^+ (protein ANM, cutoff {cutoff:g}) computed in {time.time() - t:.1f} s')
    return G

def run_b2_fast(P, res, cov, ca_opm, sets, heavy_xyz, sterol_sets, nd=1000, nd_registered=None, workers=1, cutoff=15.0, G=None, scratch=None, chunk=10, log=print, entries=None):
    nd_registered = nd if nd_registered is None else nd_registered
    t_start = time.time()
    info = dict(nd=nd, nd_registered=nd_registered, workers=workers, phases={})
    if G is None:
        G = compute_G(P, cutoff, log)
    info['phases']['G_s'] = time.time() - t_start
    pos = {r: i for i, r in enumerate(res)}
    gate_all = [pos[r] for r in sorted(sets['GATE_ALL']) if r in pos]
    gate_core = [pos[r] for r in sorted(sets['GATE_CORE']) if r in pos]
    union = np.array(sorted(set(gate_all) | set(gate_core)), int)
    col = {n: i for i, n in enumerate(union)}
    cols = dict(all=[col[n] for n in gate_all], core=[col[n] for n in gate_core])
    msf0 = F.msf_diag(G)
    m0u = msf0[union]
    anchors, outward = E.lipid_anchors(ca_opm, cov)
    entries = list(sterol_sets) if entries is None else list(entries)
    P = np.asarray(P, float)
    templates, seed_off, n_mol = ({}, {}, {})
    for e in entries:
        xyz, mol = (sterol_sets[e]['xyz'], sterol_sets[e]['mol'])
        templates[e] = [xyz[mol == k] - xyz[mol == k].mean(0) for k in np.unique(mol)]
        seed_off[e] = sum(map(ord, e)) % 1000
        n_mol[e] = len(templates[e])
    t0 = time.time()
    obs_ln, obs_nm = ({}, {})
    for e in entries:
        obs_ln[e], obs_nm[e] = F.lnratio_nodes(G, m0u, np.asarray(sterol_sets[e]['xyz'], float), P, list(CUTS), union, cut_ss=CUT_SS)
    info['phases']['observed_s'] = time.time() - t0
    LN = {e: np.empty((nd, len(CUTS), len(union))) for e in entries}
    NM = {e: np.empty((nd, len(CUTS)), int) for e in entries}
    stats_tot = {}
    t0 = time.time()
    if workers <= 1:
        F.init_worker(G, P, union, m0u, CUTS, CUT_SS, anchors, outward, heavy_xyz)
        for e in entries:
            _, Sl = F.task_place((e, seed_off[e], templates[e], nd, n_mol[e]))
            for a in range(0, nd, chunk):
                _, st, l, n, s = F.task_eval((e, a, Sl[a:a + chunk]))
                LN[e][st:st + len(l)], NM[e][st:st + len(l)] = (l, n)
                for k, v in s.items():
                    stats_tot[k] = stats_tot.get(k, 0) + v
        log(f'  inline run done in {time.time() - t0:.1f} s')
    else:
        os.makedirs(scratch, exist_ok=True) if scratch else None
        gpath = os.path.join(scratch or os.environ.get('TEMP', '.'), f'a12_b2_G_{os.getpid()}.npy')
        np.save(gpath, G)
        for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            os.environ[_v] = '1'
        ctx = mp.get_context('spawn')
        n_done = {e: 0 for e in entries}
        try:
            with cf.ProcessPoolExecutor(max_workers=workers, mp_context=ctx, initializer=F.init_worker, initargs=(gpath, P, union, m0u, CUTS, CUT_SS, anchors, outward, heavy_xyz)) as ex:
                pend = {ex.submit(F.task_place, (e, seed_off[e], templates[e], nd, n_mol[e])): ('place', e) for e in entries}
                t_place_done = {}
                while pend:
                    done, _ = cf.wait(list(pend), return_when=cf.FIRST_COMPLETED)
                    for f in done:
                        kind, e = pend.pop(f)
                        r = f.result()
                        if kind == 'place':
                            Sl = r[1]
                            t_place_done[e] = time.time() - t0
                            log(f'  [{time.time() - t0:7.1f} s] placements done: {e} ({nd} draws)')
                            for a in range(0, nd, chunk):
                                pend[ex.submit(F.task_eval, (e, a, Sl[a:a + chunk]))] = ('eval', e)
                        else:
                            _, st, l, n, s = r
                            LN[e][st:st + len(l)], NM[e][st:st + len(l)] = (l, n)
                            n_done[e] += len(l)
                            for k, v in s.items():
                                stats_tot[k] = stats_tot.get(k, 0) + v
                            if n_done[e] == nd:
                                log(f'  [{time.time() - t0:7.1f} s] evaluations done: {e} ({nd} draws x {len(CUTS)} cutoffs)')
        finally:
            try:
                os.remove(gpath)
            except OSError:
                pass
        info['placement_done_s'] = t_place_done
    info['phases']['null_s'] = time.time() - t0
    info['stats'] = stats_tot
    info['total_s'] = time.time() - t_start
    rows_reg, rows_all, detail = ([], [], {})
    for e in entries:
        for variant, key, cj in variants_for(e):
            c = cols[key]
            obs = float(np.mean(obs_ln[e][cj][c]))
            null = LN[e][:, cj, :][:, c].mean(axis=1)
            detail[e, variant] = dict(obs=obs, null=null, n_contact_nodes=int(obs_nm[e][cj]), n_contact_null=NM[e][:, cj])
            for nd_use, rows in ((nd_registered, rows_reg), (nd, rows_all)):
                nl = null[:nd_use]
                o = L.effect_summary(obs, nl)
                o.update(p=L.perm_p(obs, nl, 'lower'), n_contact_nodes=int(obs_nm[e][cj]))
                rows.append(RB._row(test='B2', structure=e, variant=variant, n_eff=len(templates[e]), obs=o['obs'], null_mean=o['null_mean'], delta=o['delta'], z=o['z'], p=o['p'], mde80=o['mde80'], note=f"contact nodes {o['n_contact_nodes']}; statistic = mean ln(MSF'/MSF) over gate nodes (lower tail)"))
    return (rows_reg, rows_all, detail, info)

def write_null_draws(detail, path):
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        fh.write('structure\tvariant\tdraw\tmean_ln_msf_ratio\tn_contact_nodes\n')
        for (e, variant), d in detail.items():
            for i, (v, n) in enumerate(zip(d['null'], d['n_contact_null'])):
                fh.write(f'{e}\t{variant}\t{i}\t{v!r}\t{int(n)}\n')

def main(argv):
    k = argv.index('--real')
    root, outdir = (argv[k + 1], argv[k + 2])

    def opt(name, default):
        return type(default)(argv[argv.index(name) + 1]) if name in argv else default
    nd_reg, nd, workers, threads = (opt('--nd-registered', 1000), opt('--nd', 2000), opt('--workers', 8), opt('--threads', 8))
    scratch = opt('--scratch', os.environ.get('TEMP', '.'))
    nd = max(nd, nd_reg)
    log_lines = []

    def log(msg=''):
        line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
        print(line, flush=True)
        log_lines.append(line)
    log(f'B2 fast run: python {platform.python_version()}, numpy {np.__version__}; main-process BLAS threads {threads}, workers {workers} (1 BLAS thread each); nd_registered {nd_reg}, nd {nd}')
    log('command line: ' + ' '.join(argv))
    log('registered scripts (SHA-256 against A12 item 12):')
    bad = check_registered_scripts(log)
    if bad:
        raise SystemExit('registered script hash mismatch: ' + ','.join(bad))
    log('fast module: a12_b2_fast.py ' + sha256_file(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a12_b2_fast.py')))
    log('driver: a12_run_b2_fast.py ' + sha256_file(os.path.abspath(__file__)))
    L.guard_amendment_appended(RB.P(root, 'preregistration', 'AMENDMENTS.md'))
    cov, sets, structs, hits = RB.load_b_inputs(root)
    qc, ster = RB.build_sterols(root, None)
    for e, v in qc.items():
        log(f"  sterol QC {e}: n_mol {v['n_mol']}, n_atoms {v['n_atoms']}, {v['D_TMD1']}, {v['D_TMD2']}")
    s = structs['7P04']
    t0 = time.time()
    rows_reg, rows_all, detail, info = run_b2_fast(s['X'], s['res'], cov, s['ca'], sets, s['heavy'], ster, nd=nd, nd_registered=nd_reg, workers=workers, scratch=scratch, log=log)
    os.makedirs(outdir, exist_ok=True)
    RB.write_rows(rows_reg, RB.P(outdir, 'a12_B2_results.tsv'))
    verdict = RB.b2_consistency(rows_reg)
    log(f'B2 written (registered null size B = {nd_reg}): {len(rows_reg)} rows; FB2: {verdict}')
    if nd > nd_reg:
        RB.write_rows(rows_all, RB.P(outdir, f'a12_B2_results_nd{nd}_SUPPLEMENT.tsv'))
        verdict2 = RB.b2_consistency(rows_all)
        log(f'SUPPLEMENT (B = {nd}, same stream, first {nd_reg} draws = registered): FB2: {verdict2}')
    write_null_draws(detail, RB.P(outdir, f'a12_B2_null_draws_nd{nd}.tsv'))
    info['verdict_registered'] = verdict
    info['wall_total_s'] = time.time() - t0
    info['threads_main'], info['workers'] = (threads, workers)
    info['versions'] = dict(python=platform.python_version(), numpy=np.__version__)
    with open(RB.P(outdir, 'a12_B2_run_info.json'), 'w', encoding='utf-8') as fh:
        json.dump(info, fh, indent=1, default=str)
    with open(RB.P(outdir, 'a12_B2_run_log.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(log_lines) + '\n')
    log('done')
if __name__ == '__main__':
    if '--real' in sys.argv:
        main(sys.argv)
    else:
        print(__doc__)

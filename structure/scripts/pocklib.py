import os, re, subprocess, shutil, csv, json
import numpy as np
import gemmi
import p3lib as L
P2R_HOME = '/home/pdr/opt/p2rank'

def kabsch(P, Q):
    P, Q = (np.asarray(P, float), np.asarray(Q, float))
    Pc, Qc = (P.mean(0), Q.mean(0))
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return (R, Qc - R @ Pc)

def ca_by_num(pdb):
    st, _ = L.read_clean(pdb)
    return {r['num']: r['ca'] for r in L.residues(st) if r['ca'] is not None}

def opm_ca():
    st = gemmi.read_structure('/home/pdr/work/phase3/in/opm_7p04.pdb')
    out = {}
    for r in st[0]['A']:
        a = r.find_atom('CA', '*')
        if a:
            out[r.seqid.num] = np.array(a.pos.tolist())
    return out

def write_local_launcher(path, xmx='1400m'):
    with open(path, 'w') as f:
        f.write(f'#!/usr/bin/env bash\nexport JAVA_OPTS="-Xmx{xmx}"\nINSTALL_DIR={P2R_HOME}\nCLASSPATH="${{INSTALL_DIR}}/bin/p2rank.jar:${{INSTALL_DIR}}/bin/lib/*"\nexec java $JAVA_OPTS -cp "${{CLASSPATH}}" cz.siret.prank.program.Main "$@"\n')
    os.chmod(path, 493)

def run_p2rank(pdb, outdir, config=None, launcher='/home/pdr/work/phase3/out/prank_local.sh', log=None):
    os.makedirs(outdir, exist_ok=True)
    cmd = [launcher, 'predict', '-f', pdb, '-o', outdir, '-threads', '1', '-visualizations', '0', '-log_to_console', '0']
    if config:
        cmd = [launcher, 'predict', '-c', config, '-f', pdb, '-o', outdir, '-threads', '1', '-visualizations', '0', '-log_to_console', '0']
    p = subprocess.run(cmd, capture_output=True, text=True)
    if log:
        open(log, 'w').write('CMD: ' + ' '.join(cmd) + '\n' + p.stdout + '\n' + p.stderr)
    if p.returncode != 0:
        raise RuntimeError('P2Rank failed: ' + p.stderr[-800:])
    base = os.path.basename(pdb)
    return os.path.join(outdir, base + '_predictions.csv')

def parse_p2rank(csvpath):
    pockets = []
    with open(csvpath) as f:
        rd = csv.reader(f, skipinitialspace=True)
        hdr = [h.strip() for h in next(rd)]
        for row in rd:
            d = dict(zip(hdr, [x.strip() for x in row]))
            res = []
            for tok in d['residue_ids'].split():
                ch, num = (tok.split('_')[0], tok.split('_')[1])
                res.append(int(re.sub('\\D', '', num)))
            pockets.append({'method': 'P2Rank', 'rank': int(d['rank']), 'score': float(d['score']), 'prob': float(d['probability']), 'center': np.array([float(d['center_x']), float(d['center_y']), float(d['center_z'])]), 'residues': sorted(set(res)), 'n_sas': int(d['sas_points']), 'extra': {}})
    return pockets

def run_fpocket(pdb, workdir, log=None):
    os.makedirs(workdir, exist_ok=True)
    base = os.path.basename(pdb)
    local = os.path.join(workdir, base)
    shutil.copyfile(pdb, local)
    p = subprocess.run(['fpocket', '-f', local], capture_output=True, text=True, cwd=workdir)
    if log:
        open(log, 'w').write(p.stdout + '\n' + p.stderr)
    if p.returncode != 0:
        raise RuntimeError('fpocket failed: ' + p.stderr[-800:])
    return os.path.join(workdir, base[:-4] + '_out')

def parse_fpocket(outdir):
    stem = os.path.basename(outdir)[:-4]
    info = open(os.path.join(outdir, stem + '_info.txt')).read()
    blocks = re.split('\\n(?=Pocket\\s+\\d+\\s*:)', info)
    pockets = []
    for b in blocks:
        m = re.match('Pocket\\s+(\\d+)\\s*:', b)
        if not m:
            continue
        n = int(m.group(1))
        kv = {}
        for line in b.splitlines()[1:]:
            if ':' in line:
                k, v = line.split(':', 1)
                try:
                    kv[k.strip()] = float(v.strip())
                except ValueError:
                    pass
        vert = os.path.join(outdir, 'pockets', f'pocket{n}_vert.pqr')
        atm = os.path.join(outdir, 'pockets', f'pocket{n}_atm.pdb')
        pts = []
        for line in open(vert):
            if line.startswith(('ATOM', 'HETATM')):
                f = line.split()
                pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        res = set()
        for line in open(atm):
            if line.startswith(('ATOM', 'HETATM')):
                res.add(int(line[22:26]))
        pockets.append({'method': 'fpocket', 'rank': n, 'score': kv.get('Score', float('nan')), 'prob': float('nan'), 'center': np.mean(pts, axis=0), 'residues': sorted(res), 'n_sas': len(pts), 'extra': kv, 'alpha_spheres': np.array(pts)})
    pockets.sort(key=lambda p: p['rank'])
    return pockets

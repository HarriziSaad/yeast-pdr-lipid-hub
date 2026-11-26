import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import p3lib as L
cov = {int(r['residue']): r for r in L.load_tsv(os.path.join(L.IN, 'structure_covariates_7P04.tsv'))}

def coords(name):
    st, _ = L.read_clean(os.path.join(L.OUT, 'pdb', f'{name}.pdb'))
    return {r['num']: r['ca'] for r in L.residues(st) if r['ca'] is not None}
refs = {r: coords(r) for r in ('7P03', '7P04', '7P05', '7P06')}
common = set.intersection(*[set(x) for x in refs.values()])

def band(tmd, lo, hi):
    out = []
    for n, r in cov.items():
        if r['TM_helix'] == '1' and r['domain'] in (('TMD1', 'ECD1') if tmd == 1 else ('TMD2', 'ECD2')) and (lo <= float(r['z_CA_OPM_A']) <= hi) and (n in common):
            out.append(n)
    return out

def cen(x, ns):
    return np.mean([x[n] for n in ns], axis=0)
print('bands (z lo hi): d between TMD1 and TMD2 centroids')
for lo, hi in [(-15.2, -10), (-10, -5), (-5, 0), (0, 5), (5, 10), (10, 15.2), (-15.2, -5), (5, 15.2), (10, 20.2)]:
    b1, b2 = (band(1, lo, hi), band(2, lo, hi))
    print(lo, hi, len(b1), len(b2), {r: round(float(np.linalg.norm(cen(refs[r], b1) - cen(refs[r], b2))), 2) for r in refs})

def d(a, b, r):
    return round(float(np.linalg.norm(refs[r][a] - refs[r][b])), 2)
for a, b in [(683, 1373), (682, 1371), (683, 1371), (682, 1373), (680, 1369), (685, 1372), (684, 1373)]:
    print('CA-CA', a, b, {r: d(a, b, r) for r in refs})
seq = ''.join((L.THREE2ONE[r.name] for r in L.read_clean(os.path.join(L.OUT, 'pdb', '7P04.pdb'))[0][0]['A']))
for n in (911,):
    pass

def kabsch(P, Q):
    Pc, Qc = (P.mean(0), Q.mean(0))
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return (R, Qc - R @ Pc)
for dom in ('NBD1', 'NBD2', 'TMD1', 'TMD2'):
    ns = [n for n, r in cov.items() if r['domain'] == dom and n in common]
    P = np.array([refs['7P06'][n] for n in ns])
    Q = np.array([refs['7P04'][n] for n in ns])
    R, t = kabsch(P, Q)
    rm = np.sqrt(np.mean(np.sum((P @ R.T + t - Q) ** 2, axis=1)))
    print('7P06 vs 7P04 superposed on', dom, 'RMSD', round(float(rm), 2), 'n', len(ns))

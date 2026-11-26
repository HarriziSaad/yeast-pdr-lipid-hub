import csv
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def read_classes(path):
    cl = {}
    with open(path, encoding='utf-8', newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            cl.setdefault(r['class'], set()).add(int(r['residue']))
    return cl

def read_truth(path):
    s = {}
    with open(path, encoding='utf-8', newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            s.setdefault(r['set'], set()).add(int(r['residue']))
    return s

def build_sets(classes, truth):
    c = classes
    hub = c['LH-MIPC-exit'] | c['LH-ergC'] | c['LH-ergH'] | truth['H6_PRIMARY']
    sets = {'R_HUB': hub, 'R_HUBGATE': hub | c['GATE-core'] | c['GATE-region'], 'R_LHCONF': c['LH-CONF'] | truth['H6_PRIMARY'], 'R_H6': set(truth['H6_PRIMARY']), 'CAV': set(c['CAV-R6G']), 'GATE_CORE': set(c['GATE-core']), 'GATE_ALL': c['GATE-core'] | c['GATE-region'], 'NBS': set(c['NBS-MOTIF'])}
    return sets

def write_sets(sets, path):
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['set', 'residue'])
        for k in sorted(sets):
            for r in sorted(sets[k]):
                w.writerow([k, r])

def size_table(sets, modelled=None):
    lines = []
    for k in sorted(sets):
        n = len(sets[k])
        m = len(sets[k] & modelled) if modelled is not None else n
        lines.append(f'{k}: n={n} modelled_in_7P04={m}')
    lines.append(f"R_HUB & CAV overlap: n={len(sets['R_HUB'] & sets['CAV'])}")
    lines.append(f"R_HUB & GATE_ALL overlap: n={len(sets['R_HUB'] & sets['GATE_ALL'])}")
    lines.append(f"R_HUBGATE & CAV overlap: n={len(sets['R_HUBGATE'] & sets['CAV'])}")
    lines.append(f"R_LHCONF & CAV overlap: n={len(sets['R_LHCONF'] & sets['CAV'])}")
    return '\n'.join(lines)
if __name__ == '__main__':
    classes_path, truth_path, cov_path, out_path = sys.argv[1:5]
    from a12_lib import read_covariates
    cov = read_covariates(cov_path)
    sets = build_sets(read_classes(classes_path), read_truth(truth_path))
    write_sets(sets, out_path)
    print(size_table(sets, set(cov)))

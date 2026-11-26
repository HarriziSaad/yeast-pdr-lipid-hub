import re, sys, html, json, itertools
txt = open(sys.argv[1], encoding='utf-8', errors='replace').read()
plain = html.unescape(re.sub('<[^>]+>', '\n', txt))
lines = [l.strip() for l in plain.split('\n') if l.strip()]
tabs = {}
cur = None
for i, l in enumerate(lines):
    m = re.match('^(\\d+)\\.\\s+(.*)$', l)
    if m and ('Code' in l or 'transl_table' in l or 'Mitochondrial' in l or ('Nuclear' in l) or ('Plastid' in l) or ('Genetic' in l)):
        cur = int(m.group(1))
        tabs[cur] = dict(title=m.group(2))
    if cur is not None and l.startswith('AAs'):
        tabs[cur]['aa'] = l.split('=', 1)[1].strip()
    if cur is not None and l.startswith('Starts'):
        tabs[cur]['starts'] = l.split('=', 1)[1].strip()
    if cur is not None and l.startswith('Base1'):
        tabs[cur]['b1'] = l.split('=', 1)[1].strip()
    if cur is not None and l.startswith('Base2'):
        tabs[cur]['b2'] = l.split('=', 1)[1].strip()
    if cur is not None and l.startswith('Base3'):
        tabs[cur]['b3'] = l.split('=', 1)[1].strip()
print('tables parsed:', len(tabs), 'ids', sorted(tabs))
std = tabs[1]

def codon(t, i):
    return t['b1'][i] + t['b2'][i] + t['b3'][i]
out = {}
for k in sorted(tabs):
    t = tabs[k]
    if 'aa' not in t:
        print(k, 'no AAs line')
        continue
    assert all((codon(t, i) == codon(std, i) for i in range(64))), k
    d = [(codon(t, i), std['aa'][i], t['aa'][i]) for i in range(64) if t['aa'][i] != std['aa'][i]]
    ctg = t['aa'][[codon(t, i) for i in range(64)].index('CTG')]
    print(k, '|', t['title'], '| CTG =', ctg, '| starts-with-CTG:', t.get('starts', '')[[codon(t, i) for i in range(64)].index('CTG')] if t.get('starts') else '', '| n diffs vs table 1:', len(d), d if ctg != 'L' else '')
    out[k] = dict(title=t['title'], aa=t['aa'], starts=t.get('starts', ''), b1=t['b1'], b2=t['b2'], b3=t['b3'])
if '--json' in sys.argv:
    json.dump(out, open(sys.argv[sys.argv.index('--json') + 1], 'w'), indent=1)

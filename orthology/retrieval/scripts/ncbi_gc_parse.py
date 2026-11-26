import re, sys, itertools, html
BASES = 'TCAG'
CODONS = [''.join(c) for c in itertools.product(BASES, repeat=3)]

def parse_prt(path):
    txt = open(path, encoding='utf-8', errors='replace').read()
    out = {}
    for blk in re.split('\\n\\s*\\{\\s*\\n\\s*name', txt)[1:]:
        names = re.findall('"([^"]+)"', blk.split('id')[0])
        m_id = re.search('\\bid\\s+(\\d+)\\s*,', blk)
        m_aa = re.search('ncbieaa\\s+"([^"]+)"', blk)
        m_st = re.search('sncbieaa\\s+"([^"]+)"', blk)
        if m_id and m_aa:
            out[int(m_id.group(1))] = dict(names=names, aa=m_aa.group(1), starts=m_st.group(1) if m_st else '')
    return out

def parse_html(path):
    txt = open(path, encoding='utf-8', errors='replace').read()
    return txt
if __name__ == '__main__':
    prt = parse_prt(sys.argv[1])
    print('gc.prt tables:', len(prt))
    std = prt[1]['aa']
    for k in sorted(prt):
        aa = prt[k]['aa']
        diffs = [(CODONS[i], std[i], aa[i]) for i in range(64) if aa[i] != std[i]]
        ctg = aa[CODONS.index('CTG')]
        flag = '  <-- CTG != Leu' if ctg != 'L' else ''
        print(k, '|'.join(prt[k]['names']), '| CTG=', ctg, '| diffs vs table 1:', diffs, flag)

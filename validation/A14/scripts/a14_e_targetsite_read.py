import argparse
import csv
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a14_p_lib as P
TARGETS = {'eL42': ('eL42_RPL42A', 56, 53, 'QTKPVFHK'), 'uL15': ('uL15_RPL28', 38, 36, 'GGQHH')}
MIN_PIDENT, MIN_QCOV, LEN_LO, LEN_HI, COPY_FRAC, MAX_COPIES, FLANK_MIN_ID = (50.0, 0.8, 0.8, 1.3, 0.7, 4, 3)

def stream_fasta(path):
    hdr, buf = (None, [])
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if hdr is not None:
                    yield (hdr, ''.join(buf))
                hdr, buf = (line[1:].split()[0], [])
            elif line:
                buf.append(line)
    if hdr is not None:
        yield (hdr, ''.join(buf))

def prefilter(merged_fasta, out_fasta, min_len=80, max_len=200):
    n_in = n_out = 0
    with open(out_fasta, 'w', encoding='utf-8') as fo:
        for h, s in stream_fasta(merged_fasta):
            n_in += 1
            if min_len <= len(s.rstrip('*')) <= max_len:
                fo.write(f">{h}\n{s.rstrip('*')}\n")
                n_out += 1
    return (n_in, n_out)

def check_reference(refs):
    for name, (label, pos, start, ctx) in TARGETS.items():
        seq = refs[label]
        if seq[start - 1:start - 1 + len(ctx)] != ctx or seq[pos - 1] != {'eL42': 'P', 'uL15': 'Q'}[name]:
            raise ValueError(f"reference {label}: residues {start}-{start + len(ctx) - 1} are '{seq[start - 1:start - 1 + len(ctx)]}', expected '{ctx}'")

def read_position(ref, hit, pos):
    ra, rb = P.nw_align(ref, hit)
    k = 0
    cols = {}
    for x, y in zip(ra, rb):
        if x != '-':
            k += 1
            cols[k] = y
    res = cols.get(pos)
    ctx = ''.join((cols.get(p, '-') for p in range(pos - 3, pos + 4)))
    ident = sum((1 for p in range(pos - 3, pos + 4) if p != pos and cols.get(p) == ref[p - 1]))
    return (None if res in (None, '-') else res, ctx, ident)

def build_table(hits_tsv, small_fasta, refs, ctg_calls, species_list=None):
    check_reference(refs)
    hits = {}
    with open(hits_tsv, encoding='utf-8', newline='') as fh:
        for r in csv.reader(fh, delimiter='\t'):
            q, s, pid, aln, qlen, slen, qs, qe, ev, bit = (r[0], r[1], float(r[2]), int(r[3]), int(r[4]), int(r[5]), int(r[6]), int(r[7]), float(r[8]), float(r[9]))
            if pid < MIN_PIDENT or (qe - qs + 1) / qlen < MIN_QCOV or (not LEN_LO * qlen <= slen <= LEN_HI * qlen):
                continue
            hits.setdefault((q, s.split('|')[0]), []).append((bit, s))
    wanted = {s for v in hits.values() for _, s in v}
    seqs = {h: s for h, s in stream_fasta(small_fasta) if h in wanted}
    species = set(species_list or [])
    species.update((sp for _, sp in hits))
    rows = []
    for sp in sorted(species):
        row = {'species': sp}
        for name, (label, pos, _st, _ctx) in TARGETS.items():
            cand = sorted(hits.get((label, sp), []), key=lambda t: -t[0])
            res_list, ctx_list, flags = ([], [], [])
            best = cand[0][0] if cand else None
            for bit, sid in cand:
                if len(res_list) >= MAX_COPIES or bit < COPY_FRAC * best:
                    continue
                res, ctx, ident = read_position(refs[label], seqs[sid], pos)
                if ident < FLANK_MIN_ID:
                    flags.append('FLANK_MISMATCH')
                    continue
                if res is None:
                    flags.append('GAP_AT_POSITION')
                    continue
                if res == 'L' and ctg_calls.get(sp, 'L') != 'L':
                    flags.append('CTG_AMBIGUOUS')
                    res = 'NA'
                res_list.append(res)
                ctx_list.append(ctx)
            if not res_list or all((x == 'NA' for x in res_list)):
                flags.append('NO_HIT' if not cand else 'NO_USABLE_COPY')
                row[f'{name}_{pos}'] = 'NA'
            else:
                row[f'{name}_{pos}'] = ';'.join((x for x in res_list if x != 'NA'))
            row[f'{name}_n_copies'] = len(res_list)
            row[f'{name}_context'] = ';'.join(ctx_list)
            row[f'{name}_flags'] = ','.join(sorted(set(flags)))
        rows.append(row)
    return rows

def read_ctg_calls(index_tsv):
    out = {}
    with open(index_tsv, encoding='utf-8', newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            out[r['speciesCode']] = r.get('codetta_CTG', 'L') or 'L'
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('prefilter')
    a.add_argument('--merged-fasta', required=True)
    a.add_argument('--out', required=True)
    b = sub.add_parser('table')
    b.add_argument('--hits', required=True)
    b.add_argument('--small-fasta', required=True)
    b.add_argument('--refs', required=True)
    b.add_argument('--index', required=True)
    b.add_argument('--out', required=True)
    args = ap.parse_args()
    if args.cmd == 'prefilter':
        n_in, n_out = prefilter(args.merged_fasta, args.out)
        print(f'proteins read {n_in}; kept {n_out} (80-200 aa)')
        return
    refs = {h: s for h, s in stream_fasta(args.refs)}
    ctg = read_ctg_calls(args.index)
    rows = build_table(args.hits, args.small_fasta, refs, ctg, species_list=sorted(ctg))
    cols = ['species', 'eL42_56', 'eL42_n_copies', 'eL42_context', 'eL42_flags', 'uL15_38', 'uL15_n_copies', 'uL15_context', 'uL15_flags']
    with open(args.out, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join((str(r.get(c, '')) for c in cols)) + '\n')
    n_na = sum((1 for r in rows if r['eL42_56'] == 'NA' or r['uL15_38'] == 'NA'))
    print(f'genomes {len(rows)}; with an NA genotype {n_na}')
if __name__ == '__main__':
    main()

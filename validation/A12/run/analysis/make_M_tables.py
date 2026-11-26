import csv, sys, os
run = sys.argv[1]
rows = list(csv.DictReader(open(os.path.join(run, 'M', 'a12_M_results.tsv'), encoding='utf-8'), delimiter='\t'))

def g(x, n=4):
    try:
        v = float(x)
        return 'nan' if v != v else ('%.' + str(n) + 'g') % v
    except Exception:
        return x if x else '-'
out = []
out.append('| test | variant | region | hit set | hit res | n_eff | obs | null mean (M1b: cavity excess) | Delta | z | p | q | 95% CI Delta | MDE80 | fallback | weak | balance flag |')
out.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in rows:
    ci = '' if r['ci_lo'] in ('', 'nan') else '[%s, %s]' % (g(r['ci_lo'], 3), g(r['ci_hi'], 3))
    out.append('| ' + ' | '.join([r['test'], r['variant'], r['region'], r['hits'], r['n_hit_residues'], str(r['n_eff']), g(r['obs']), g(r['null_mean']), g(r['delta']), g(r['z'], 3), g(r['p']), g(r['q']), ci or '-', g(r['mde80'], 3), g(r['frac_fallback'], 3), r['weakly_matched'] or '-', r['balance_flag'] or '-']) + ' |')
open(os.path.join(run, 'analysis', 'M_results_table.md'), 'w', encoding='utf-8').write('\n'.join(out) + '\n')
print('ok', len(rows))

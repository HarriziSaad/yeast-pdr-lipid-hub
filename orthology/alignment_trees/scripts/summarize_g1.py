import json, os, sys, csv
d = sys.argv[1]
s = json.load(open(os.path.join(d, 'g1_summary.json')))
rd = lambda f: list(csv.DictReader(open(os.path.join(d, f)), delimiter='\t'))
print('Tree:', s['tree'], '| UFBoot present:', s['ufb_present'])
print('Rooting:', s['rooting'])
print()
print('| clade | seqs | species | WGD sp. | stem SH-aLRT | stem UFB | support ok | scoring group (3f) | best-hit assigned | violators | rate | 1T | 1L |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in rd('g1_clades.tsv'):
    print('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (r['clade'], r['n_sequences'], r['n_species'], r['n_wgd_species'], r['stem_SH_aLRT'], r['stem_UFBoot'], r['support_ok'], r['scoring_group_3f'], r['assigned_best_hit_in_universe'], r['violators_outside_clade'], r['violator_rate'], r['pass_1T'], r['pass_1L']))
print()
print('Pairs:')
for k, v in s['pairs'].items():
    print('-', k, v)
print()
print('YGOB counts:', s['ygob_counts'], '| failures:', s['ygob_failures'], '| two-copy failures:', s['two_copy_failures'])
print('S288C consistency:', s['s288c_consistency'])
print()
for k in ('criterion1_T_all_five', 'criterion1_L_all_five', 'criterion2_pairs', 'criterion2_ygob', 'n_clades_pass_1T', 'n_clades_pass_1L', 'criterion3_T', 'criterion3_L', 'G1_PASS_T', 'G1_PASS_L'):
    print(k, '=', s[k])

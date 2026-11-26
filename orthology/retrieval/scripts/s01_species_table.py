import argparse, collections
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from common import read_tsv, write_tsv, species_key_of, log
CLASS = {'L': 'standard', 'S': 'CUG-Ser', 'A': 'CUG-Ala', '?': 'CUG-undetermined', '': 'unknown'}
ap = argparse.ArgumentParser()
ap.add_argument('--map', required=True)
ap.add_argument('--out', required=True)
a = ap.parse_args()
rows = read_tsv(a.map)
out = []
for r in rows:
    st = r['name_match_status']
    if st == 'NAMES_ROW_WITHOUT_PEP_FILE':
        continue
    if st == 'EXTRA_NOT_FINAL_PEP':
        stem = r['pep_stem']
        code = {'candida_albicans.cgd': 'candida_albicans_CGDref', 'saccharomyces_cerevisiae.sgd': 'saccharomyces_cerevisiae_SGDref'}.get(stem, stem)
        out.append(dict(speciesCode=code, kind='reference_annotation', pep_stem=stem, pep_file=r['pep_file'], assembly_fullID_updated='', name_match_status=st, FINAL_Tip_ID='', species_key='Candida_albicans' if 'albicans' in stem else 'Saccharomyces_cerevisiae', strain='SC5314 (CGD)' if 'albicans' in stem else 'S288C (SGD)', genus='Candida' if 'albicans' in stem else 'Saccharomyces', order='', family_2018='', ctg_call='S' if 'albicans' in stem else 'L', ctg_class='CUG-Ser' if 'albicans' in stem else 'standard'))
        continue
    ftid = r['FINAL_Tip_ID']
    sk, strain = species_key_of(ftid) if ftid else (r['speciesCode'], '')
    call = r['codetta_CTG']
    out.append(dict(speciesCode=r['speciesCode'], kind='y1000_genome', pep_stem=r['pep_stem'], pep_file=r['pep_file'], assembly_fullID_updated=r['assembly_fullID_updated'], name_match_status=st, FINAL_Tip_ID=ftid, species_key=sk, strain=strain, genus=ftid.split('_')[0] if ftid else '', order=r['order_2024_itol'], family_2018=r['family_2018'], ctg_call=call, ctg_class=CLASS.get(call, 'unknown:' + call), codetta_diffs_vs_table1=r['codetta_diffs_vs_table1']))
cnt = collections.Counter((r['species_key'] for r in out if r['kind'] == 'y1000_genome'))
for r in out:
    r['n_genomes_of_species'] = cnt.get(r['species_key'], '') if r['kind'] == 'y1000_genome' else ''
fields = ['speciesCode', 'kind', 'pep_stem', 'pep_file', 'assembly_fullID_updated', 'name_match_status', 'FINAL_Tip_ID', 'species_key', 'strain', 'genus', 'order', 'family_2018', 'ctg_call', 'ctg_class', 'codetta_diffs_vs_table1', 'n_genomes_of_species']
write_tsv(a.out, out, fields)
g = [r for r in out if r['kind'] == 'y1000_genome']
log('genomes', len(g), 'distinct species_key', len(cnt), 'species with >1 genome', sum((1 for v in cnt.values() if v > 1)), 'ctg classes', dict(collections.Counter((r['ctg_class'] for r in g))), 'orders', dict(collections.Counter((r['order'] for r in g))))

#!/bin/bash
n=$1
W=/home/pdr/work/phase1b
S=/mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/scripts
O=/mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/outputs
until grep -q 'Total wall-clock time used' $W/full/${n}_T1/full_${n}_T1.iqtree 2>/dev/null; do sleep 60; done
mkdir -p $W/g1/full_${n}_T1 $O/g1/full_${n}_T1 $O/trees
python -I $S/b05_g1_eval.py --tree $W/full/${n}_T1/full_${n}_T1.treefile --outdir $W/g1/full_${n}_T1 > $W/logs/g1_full_${n}_T1.log 2>&1
python -I - <<PY > $W/g1/full_${n}_T1/PAIR_NODE_SUPPORT.txt
import json, csv
d = json.load(open('$W/g1/full_${n}_T1/g1_summary.json'))
print('Full T1 tree of the $n alignment: node supports from IQ-TREE (SH-aLRT / UFBoot). Descriptive; NOT part of the A4 gate (A4 verdicts are those of the fast trees).')
for k, v in d['pairs'].items():
    print('sister-pair node', k, ': sisters =', v['sister'], '| SH-aLRT =', v['sh'], '| UFBoot =', v['ufb'], '| sequences', v['n_seq'], '| species', v['n_species'])
for r in csv.DictReader(open('$W/g1/full_${n}_T1/g1_clades.tsv'), delimiter='\t'):
    print('clade stem', r['clade'], ': sequences', r['n_sequences'], '| SH-aLRT', r['stem_SH_aLRT'], '| UFBoot', r['stem_UFBoot'])
print('backbone:', open('$W/g1/full_${n}_T1/g1_backbone.txt').read().strip())
print('(A4 verdict flags computed on this tree, for information only) G1_PASS_T =', d['G1_PASS_T'], ', G1_PASS_L =', d['G1_PASS_L'])
PY
cp $W/g1/full_${n}_T1/* $O/g1/full_${n}_T1/
cp $W/full/${n}_T1/full_${n}_T1.treefile $W/full/${n}_T1/full_${n}_T1.iqtree $W/full/${n}_T1/full_${n}_T1.contree $O/trees/ 2>/dev/null
echo "$(date -u +%FT%TZ) post-processing $n done"

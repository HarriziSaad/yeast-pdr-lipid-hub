#!/bin/bash
set -euo pipefail
P=/mnt/c/Users/MSI/Documents/Advisor/project
for A in mafft muscle; do
  python3 -I -X utf8 $P/phase1_orthology/phase1b/scripts/b05_g1_eval.py \
    --tree /home/pdr/work/phase1c/$A/fast_sacch_$A.treefile \
    --manifest $P/phase1_orthology/phase1c/outputs/subset_manifest.tsv \
    --ygob-map $P/phase1_orthology/phase1b/outputs/set/ygob/ygob_map.tsv \
    --outdir $P/phase1_orthology/phase1c/outputs/g1_sacch/fast_$A
  echo "== done $A"
done

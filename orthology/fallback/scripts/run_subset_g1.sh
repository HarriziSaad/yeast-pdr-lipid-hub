#!/bin/bash
set -euo pipefail
A=$1; T=${2:-4}
W=/home/pdr/work/phase1c/$A; mkdir -p "$W"; cd "$W"
P1B=/home/pdr/work/phase1b
IDS=/mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1c/outputs/subset_ids.txt
SCR=/mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1c/scripts
python3 "$SCR/c02_extract_rows.py" --aln "$P1B/aln/$A.trim50.faa" --ids "$IDS" --out "$A.sacch.trim50.faa" --report "$A.sacch.trim50.json" \
        --expect-cols $([ "$A" = mafft ] && echo 1444 || echo 1477)
iqtree2 -s "$A.sacch.trim50.faa" -m Q.yeast+G4 -fast -alrt 1000 -T "$T" --mem 4G --seed 20261118 -pre "fast_sacch_$A"

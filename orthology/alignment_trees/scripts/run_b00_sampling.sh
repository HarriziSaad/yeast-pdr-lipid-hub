#!/bin/bash
set -euo pipefail
P=/mnt/c/Users/MSI/Documents/Advisor/project
S=$P/phase1_orthology/phase1b/scripts
M=$P/phase1_orthology/data_staging/metadata
W=/home/pdr/work/phase1b/set
mkdir -p $W/sampling
python -I $S/b00_sampling_plan.py --names $M/ALL_NAMES_FILE.tsv --itol $M/itol_color_strip.txt --tree $M/1154yeasts_1403OGs_ml_timetree.tree \
  --dias $P/phase1_orthology/phase1a/inputs/dias2025_tableA1_strains.tsv --retrieval /home/pdr/work/phase1a/run_full/qc/retrieval_table.tsv \
  --rep-col cdhit99_rep_within_species --wgd-mode all --outdir $W/sampling

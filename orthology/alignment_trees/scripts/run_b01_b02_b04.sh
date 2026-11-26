#!/bin/bash
set -e
cd /home/pdr/work/phase1b
S=/mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/scripts
python -I $S/b01_build_set.py --plan set/sampling/sampling_plan.tsv --outdir set
python -I $S/b02_labels.py --setdir set --seeds /mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/data_staging/seeds_scerevisiae/scerevisiae_S288C_PDR8_uniprot.fasta
python -I $S/b04_ygob_map.py

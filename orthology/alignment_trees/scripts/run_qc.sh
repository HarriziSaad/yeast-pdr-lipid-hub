#!/bin/bash
cd /home/pdr/work/phase1b
python -I /mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/scripts/b03_aln_qc.py --aln $2 --name $1 --outdir aln/qc_$1
mv aln/qc_$1/$1.trim50.faa aln/$1.trim50.faa; mv aln/qc_$1/$1.kept_columns.txt aln/$1.kept_columns.txt

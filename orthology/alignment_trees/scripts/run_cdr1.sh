#!/bin/bash
cd /home/pdr/work/phase1b
python -I /mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/scripts/b06_cdr1_check.py --aln aln/$2 --name $1 --usalign struct/us_9VT1_7P04.txt --cdr1-cif struct/9VT1.cif --pdr5-cif struct/7P04.cif --outdir aln/qc_$1

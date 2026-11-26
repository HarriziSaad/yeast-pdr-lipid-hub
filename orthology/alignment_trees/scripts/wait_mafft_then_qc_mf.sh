#!/bin/bash
log=/home/pdr/work/phase1b/logs/mafft.stderr.log
until grep -q '^exit=' $log 2>/dev/null; do sleep 15; done
S=/mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/scripts
bash $S/run_qc.sh mafft /home/pdr/work/phase1b/aln/mafft.aln.faa > /home/pdr/work/phase1b/logs/qc_mafft.log 2>&1
bash $S/run_cdr1.sh mafft mafft.aln.faa > /home/pdr/work/phase1b/logs/cdr1_mafft.log 2>&1
bash $S/run_mf.sh mafft 8 > /home/pdr/work/phase1b/logs/mf_mafft.log 2>&1

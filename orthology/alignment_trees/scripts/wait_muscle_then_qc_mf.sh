#!/bin/bash
log=/home/pdr/work/phase1b/logs/muscle.stderr.log
until grep -q '^exit=' $log 2>/dev/null; do sleep 20; done
S=/mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/scripts
bash $S/run_qc.sh muscle /home/pdr/work/phase1b/aln/muscle.super5.afa > /home/pdr/work/phase1b/logs/qc_muscle.log 2>&1
bash $S/run_cdr1.sh muscle muscle.super5.afa > /home/pdr/work/phase1b/logs/cdr1_muscle.log 2>&1
bash $S/run_mf.sh muscle 6 > /home/pdr/work/phase1b/logs/mf_muscle.log 2>&1

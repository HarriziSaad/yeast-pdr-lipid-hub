#!/bin/bash
cd /mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/scripts
sha256sum b03_aln_qc.py b04_ygob_map.py b05_g1_eval.py b06_cdr1_check.py treeutil.py run_align.sh run_mf.sh run_fast.sh run_full.sh
sha256sum /home/pdr/work/phase1b/aln/mafft.aln.faa /home/pdr/work/phase1b/aln/mafft.trim50.faa /home/pdr/work/phase1b/set/ygob/ygob_map.tsv
sha256sum /mnt/c/Users/MSI/Documents/Advisor/project/preregistration/AMENDMENTS.md /mnt/c/Users/MSI/Documents/Advisor/project/preregistration/FREEZE.txt
date -u +%Y-%m-%dT%H:%M:%SZ

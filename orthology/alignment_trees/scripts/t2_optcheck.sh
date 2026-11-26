#!/bin/bash
W=/home/pdr/work/phase1b; mkdir -p $W/tmp/t2opt; cd $W/tmp/t2opt
timeout 120 iqtree2 -s $W/aln/mafft.trim50.faa -m LG+C20+F+G -ft $W/trees/mafft/fast_mafft.treefile -B 1000 -alrt 1000 --wbtl -T 1 --mem 100M --seed 20261118 -pre t2opt -redo 2>&1 | grep -v -E '^[0-9]+ +[A-Za-z]' | tail -12 | cut -c1-200

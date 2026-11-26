#!/bin/bash
W=/home/pdr/work/phase1b
n=$1; T=${2:-6}
mkdir -p $W/trees/$n; cd $W/trees/$n
iqtree2 -s $W/aln/$n.trim50.faa -m MF -mset LG,WAG,JTT,Q.yeast -mfreq FU,F -mrate E,G -T $T --mem 3G --seed 20261118 -pre mf_$n

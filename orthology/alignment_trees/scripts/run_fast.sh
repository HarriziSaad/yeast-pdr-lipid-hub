#!/bin/bash
W=/home/pdr/work/phase1b
n=$1; M=$2; T=${3:-6}
mkdir -p $W/trees/$n; cd $W/trees/$n
iqtree2 -s $W/aln/$n.trim50.faa -m "$M" -fast -alrt 1000 -T $T --mem 4G --seed 20261118 -pre fast_$n

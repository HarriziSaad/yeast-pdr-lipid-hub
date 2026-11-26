#!/bin/bash
W=/home/pdr/work/phase1b
n=$1; tag=$2; M=$3; T=${4:-6}; G=${5:-}
mkdir -p $W/full/${n}_$tag; cd $W/full/${n}_$tag
if [ -n "$G" ]; then
  iqtree2 -s $W/aln/$n.trim50.faa -m "$M" -ft "$G" -B 1000 -alrt 1000 -T $T --mem 5G --seed 20261118 -pre full_${n}_$tag
else
  iqtree2 -s $W/aln/$n.trim50.faa -m "$M" -B 1000 -alrt 1000 -T $T --mem 5G --seed 20261118 -pre full_${n}_$tag
fi

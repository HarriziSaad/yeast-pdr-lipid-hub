#!/bin/bash
n=$1; T=${2:-6}
W=/home/pdr/work/phase1b
M=$(awk '{print $NF}' $W/logs/model_$n.txt | head -1)
echo "$(date -u +%FT%TZ) lane $n model $M threads $T"
for tag in ${STAGES:-T1}; do
  d=$W/full/${n}_$tag; mkdir -p $d; cd $d
  if [ -s full_${n}_$tag.iqtree ] && grep -q 'Total wall-clock time used' full_${n}_$tag.iqtree; then echo "$tag already finished"; continue; fi
  echo "$(date -u +%FT%TZ) start $tag"
  if [ $tag = T1 ]; then
    iqtree2 -s $W/aln/$n.trim50.faa -m "$M" -B 1000 -alrt 1000 --wbtl -T $T --mem ${MEM:-3G} --seed 20261118 -pre full_${n}_T1
  else
    iqtree2 -s $W/aln/$n.trim50.faa -m LG+C20+F+G -ft $W/trees/$n/fast_$n.treefile -B 1000 -alrt 1000 --wbtl -T $T --mem ${MEM:-3G} --seed 20261118 -pre full_${n}_T2
  fi
  echo "$(date -u +%FT%TZ) end $tag exit=$?"
done

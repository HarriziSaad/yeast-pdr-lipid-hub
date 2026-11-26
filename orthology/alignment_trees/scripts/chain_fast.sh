#!/bin/bash
n=$1; T=${2:-6}
W=/home/pdr/work/phase1b; S=/mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/scripts
f=$W/trees/$n/mf_$n.iqtree
until [ -s $f ] && grep -q 'Best-fit model according to BIC' $f; do sleep 30; done
M=$(grep 'Best-fit model according to BIC' $f | sed 's/.*BIC: *//' | awk '{print $1}')
echo "$(date -u +%FT%TZ) model for $n: $M" | tee $W/logs/model_$n.txt
bash $S/run_fast.sh $n "$M" $T > $W/logs/fast_$n.log 2>&1
mkdir -p $W/g1/fast_$n
python -I $S/b05_g1_eval.py --tree $W/trees/$n/fast_$n.treefile --outdir $W/g1/fast_$n > $W/logs/g1_fast_$n.log 2>&1
echo "$(date -u +%FT%TZ) done fast+G1 for $n" >> $W/logs/model_$n.txt

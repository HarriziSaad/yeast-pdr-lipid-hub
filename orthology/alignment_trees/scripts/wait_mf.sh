#!/bin/bash
n=$1; max=${2:-160}
f=$W/trees/$n/mf_$n.log; f=/home/pdr/work/phase1b/trees/$n/mf_$n.log
t=0
while [ $t -lt $max ]; do
  grep -q 'Best-fit model' /home/pdr/work/phase1b/trees/$n/mf_$n.iqtree 2>/dev/null && break
  sleep 5; t=$((t+5))
done
echo "$(date -u +%T) models done: $(grep -c -E '^ *[0-9]+ +[A-Za-z]' $f)"
grep -E '^ *[0-9]+ +[A-Za-z]' $f | tail -2 | cut -c1-110
grep 'Best-fit' /home/pdr/work/phase1b/trees/$n/mf_$n.iqtree 2>/dev/null

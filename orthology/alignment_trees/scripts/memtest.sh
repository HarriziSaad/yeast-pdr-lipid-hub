#!/bin/bash
W=/home/pdr/work/phase1b; mkdir -p $W/tmp; cd $W/tmp
for m in "LG+C20+F+G" "LG+C10+F+G" "LG+F+G"; do
  echo "== $m"; iqtree2 -s $W/aln/$1.trim50.faa -m "$m" --mem 50M -T 1 -n 0 -pre memtest_$m -redo 2>&1 | grep -i -E "RAM|memory|MB|GB" | head -5
done

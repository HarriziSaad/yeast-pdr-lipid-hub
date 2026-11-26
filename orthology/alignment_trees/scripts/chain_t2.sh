#!/bin/bash
W=/home/pdr/work/phase1b
S=/mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/scripts
fin() { grep -q 'Total wall-clock time used' $W/full/$1_T1/full_$1_T1.iqtree 2>/dev/null; }
until fin mafft && fin muscle; do sleep 60; done
for n in mafft muscle; do
  until [ "$(awk '/MemAvailable/ {print int($2/1024)}' /proc/meminfo)" -ge 7000 ] && ! pgrep -x iqtree2 >/dev/null && ! pgrep -x iqtree3 >/dev/null; do sleep 60; done
  echo "$(date -u +%FT%TZ) starting T2 $n"
  STAGES=T2 MEM=6G bash $S/full_lane.sh $n 10
  echo "$(date -u +%FT%TZ) T2 $n finished/stopped"
done

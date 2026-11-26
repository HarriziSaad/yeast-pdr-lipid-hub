#!/bin/bash
echo "DEPRECATED by the coordinator (2026-11-19): T2 must run alone with 6G; use chain_t2.sh" >&2; exit 1
n=$1
W=/home/pdr/work/phase1b
S=/mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/scripts
until grep -q 'Total wall-clock time used' $W/full/${n}_T1/full_${n}_T1.iqtree 2>/dev/null; do sleep 60; done
echo "$(date -u +%FT%TZ) T1 $n finished; starting T2 $n"
STAGES=T2 MEM=3G bash $S/full_lane.sh $n 6
echo "$(date -u +%FT%TZ) T2 $n ended"

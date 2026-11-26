#!/bin/bash
set -u
W=/home/pdr/work/phase1b
S=$W/set/analysis_set.faa
mkdir -p $W/aln $W/logs
cd $W/aln
case "$1" in
  mafft)
    /usr/bin/time -v mafft --retree 2 --maxiterate 2 --thread ${2:-8} $S > mafft.aln.faa 2> $W/logs/mafft.stderr.log
    echo "exit=$?" >> $W/logs/mafft.stderr.log ;;
  muscle)
    /usr/bin/time -v muscle -super5 $S -output muscle.super5.afa -threads ${2:-4} > $W/logs/muscle.stdout.log 2> $W/logs/muscle.stderr.log
    echo "exit=$?" >> $W/logs/muscle.stderr.log ;;
esac
